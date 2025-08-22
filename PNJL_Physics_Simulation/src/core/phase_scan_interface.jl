module PhaseScanInterface

"""
参数轴结构：表示一个参数的采样点
字段：name, values, unit（可选）
"""
struct ParameterAxis
    name::String
    values::Vector{Float64}
    unit::Union{String, Nothing}
    ParameterAxis(name::String, values::Vector{Float64}; unit::Union{String, Nothing}=nothing) = new(name, values, unit)
end

"""
生成网格迭代器，返回每个格点的 NamedTuple（参数名 => 值）。
"""
function generate_grid(axes::Vector{ParameterAxis})
    names_vec = [a.name for a in axes]
    names = Tuple(Symbol.(names_vec))
    arrays = [a.values for a in axes]
    sizes = map(length, arrays)
    inds = CartesianIndices(Tuple(sizes))
    return (NamedTuple{names}(Tuple(arrays[i][I[i]] for i=1:length(arrays))) for I in inds)
end

"""
将项目的物理计算函数适配为 scanner 使用的统一接口：接受 NamedTuple 并返回标量值。
"""
function wrap_physics_fn(fn::Function)
    return params -> try
        fn(params)
    catch e
        try
            if haskey(params, :T) && haskey(params, :mu)
                return fn(params.T, params.mu)
            else
                return fn(params...)
            end
        catch
            rethrow(e)
        end
    end
end

"""
扫描参数空间并将结果写为 CSV 文件（首选）。

参数：
 - f: 接受 NamedTuple 的函数，返回标量 Omega
 - axes: ParameterAxis 数组
可选参数：outfile（CSV 路径），batchsize（保留，当前实现一次性写入），返回值为 outfile
"""
function scan_phase_space(f::Function, axes::Vector{ParameterAxis};
    mode::Symbol = :stream,
    outfile::AbstractString = "output/phase_scan.csv",
    outputs::Vector{Symbol} = [:Omega],
    solver::Union{Nothing, Function} = nothing,
    initial_guess::Union{Nothing, Function} = nothing,
    postprocess::Union{Nothing, Function} = nothing,
    solution_names::Union{Nothing, Vector{Symbol}} = nothing,
    header_probe::Symbol = :first_point, # :first_point or :first_batch
    probe_batchsize::Int = 16,
    solver_options::Dict = Dict(),
    batchsize::Int = 256,
    max_points::Int = 1_000_000,
    fill_failed_with_last::Bool = true)

    # 检查点数上限
    total_points = prod(map(a->length(a.values), axes))
    if total_points > max_points
        error("total_points ($(total_points)) > max_points ($(max_points)); 减小网格或使用自适应采样")
    end

    # 输出目录准备
    outdir = dirname(outfile)
    if !isdir(outdir)
        try mkpath(outdir) catch end
    end

    # 内部：写 CSV 行头（stream 模式）
    function _open_stream_writer(outfile, header_cols::Vector{String})
        io = open(outfile, "w")
        println(io, join(header_cols, ","))
        return io
    end

    # 内部：写入 meta.json（尽量使用 JSON.jl，否则手写）
    function _write_meta(outfile, meta::Dict)
        metafile = string(outfile, ".meta.json")
        try
            @eval using JSON
            open(metafile, "w") do io
                print(io, JSON.json(meta))
            end
        catch
            # 简单序列化为最小 JSON（仅包含简单键值）
            open(metafile, "w") do io
                print(io, "{")
                first = true
                for (k,v) in meta
                    if !first; print(io, ", ") end
                    print(io, "\"", string(k), "\": ")
                    if isa(v, String)
                        print(io, "\"", v, "\"")
                    else
                        print(io, v)
                    end
                    first = false
                end
                print(io, "}")
            end
        end
    end

    # 确保能使用 Dates.now() 获取时间字符串
    try
        @eval using Dates
    catch
    end

    # 生成列头：参数列按 axes 顺序，然后 outputs 展开，最后 solver meta
    param_names = [a.name for a in axes]
    param_syms = Symbol.(param_names)
    header_cols = String[]
    append!(header_cols, param_names)
    # outputs columns will be determined on-the-fly if vectors present
    # 记录向量输出长度（首次出现时固定长度）
    vector_output_length = Dict{Symbol, Int}()
    append!(header_cols, string.(outputs))
    push!(header_cols, "solver_converged")
    push!(header_cols, "n_iter")
    push!(header_cols, "residual")

    # 打开文件和写头（stream 模式）
    io = nothing

    # 用于统计与填充
    prev_solution_outvals = nothing    # 用于填充 CSV 的上一个输出字典
    prev_solution_vector = nothing     # 用于 warm-start 的上一个解向量（具体格式由 solver 定义）
    prev_success = false
    global_min = (val=Inf, params=nothing)
    failed_points = Vector{Any}()

    # helper: build params NamedTuple from index vector
    function params_from_idx(idx)
        values = [axes[i].values[idx[i]] for i in 1:length(axes)]
        names = Tuple(Symbol.(param_names))
        return NamedTuple{names}(Tuple(values))
    end

    # helper: run solver/initial_guess/postprocess for a params and return (outvals::Dict, solver_meta, success::Bool, solution_vector)
    function evaluate_point(params, prev_vec)
        g0 = nothing
        if initial_guess !== nothing
            try
                g0 = initial_guess(params)
            catch
                g0 = nothing
            end
        end
        if g0 === nothing && prev_vec !== nothing
            g0 = prev_vec
        end

        outvals = Dict{Symbol, Any}()
        solver_meta = (converged=true, niter=0, residual=0.0)
        success = false
        solvec = nothing
        try
            if solver !== nothing
                res = solver(params, g0; solver_options...)
                if isa(res, NamedTuple)
                    for k in keys(res)
                        outvals[Symbol(k)] = getfield(res, k)
                    end
                    # if there's a :zero field, treat as solution vector
                    if haskey(outvals, :zero)
                        solvec = outvals[:zero]
                    else
                        # if solution_names provided and fields match, assemble vector
                        if solution_names !== nothing
                            ok = all(n -> haskey(outvals, n), solution_names)
                            if ok
                                solvec = [outvals[n] for n in solution_names]
                            end
                        end
                    end
                elseif isa(res, Dict)
                    for (k,v) in res
                        outvals[Symbol(k)] = v
                    end
                    if haskey(outvals, :zero)
                        solvec = outvals[:zero]
                    end
                else
                    outvals[:Omega] = res
                end

                solver_meta = (converged=get(outvals, :converged, true), niter=get(outvals, :niter, 0), residual=get(outvals, :residual, 0.0))
                success = Bool(solver_meta.converged)
            else
                val = try f(params) catch e; NaN end
                outvals[:Omega] = val
                success = !(isnan(outvals[:Omega]) || isinf(outvals[:Omega]))
            end
        catch e
            success = false
            solver_meta = (converged=false, niter=0, residual=Inf)
        end

        # postprocess hook
        if success && postprocess !== nothing
            try
                pp = postprocess(solvec === nothing ? outvals : solvec, params)
                if isa(pp, NamedTuple)
                    for k in keys(pp)
                        outvals[Symbol(k)] = getfield(pp, k)
                    end
                elseif isa(pp, Dict)
                    for (k,v) in pp
                        outvals[Symbol(k)] = v
                    end
                end
            catch
                # postprocess failed -> mark failure
                success = false
                solver_meta = (converged=false, niter=get(solver_meta, :niter, 0), residual=get(solver_meta, :residual, Inf))
            end
        end

        return outvals, solver_meta, success, solvec
    end

    # 递归嵌套循环以实现按 axes 顺序的行优先遍历（axes[1] 为行外层）
    N = length(axes)
    idx = fill(1, N)

    # header probing (determine full header and possibly write first rows)
    first_params = NamedTuple{Tuple(Symbol.(param_names))}(Tuple(a.values[1] for a in axes))
    header_determined = false
    first_rows_buffer = Vector{Vector{Any}}()
    written_first_params = false
    if mode == :stream
        if header_probe == :first_point
            outvals_p, meta_p, succ_p, solvec_p = evaluate_point(first_params, nothing)
            if !succ_p && fill_failed_with_last
                # nothing to fill with; mark failure
                push!(failed_points, (params=first_params, reason="solver_or_postprocess_failed"))
            end
            # build header based on outvals_p and outputs/solution_names
            # start with parameter names
            header_cols = String[]
            append!(header_cols, param_names)
            # determine output columns from outputs and outvals_p
            for key in outputs
                v = get(outvals_p, key, nothing)
                if v === nothing
                    push!(header_cols, string(key))
                elseif isa(v, AbstractVector)
                    if solution_names !== nothing && key == :zero
                        for n in solution_names
                            push!(header_cols, string(n))
                        end
                    else
                        for i in 1:length(v)
                            push!(header_cols, string(key, "_", i))
                        end
                    end
                else
                    push!(header_cols, string(key))
                end
            end
            push!(header_cols, "solver_converged")
            push!(header_cols, "n_iter")
            push!(header_cols, "residual")

            # open and write header
            io = _open_stream_writer(outfile, header_cols)
            # write first row if successful or per fill strategy
            row = Any[]
            for i in 1:length(param_names)
                push!(row, getfield(first_params, param_syms[i]))
            end
            for key in outputs
                v = get(outvals_p, key, NaN)
                if isa(v, AbstractVector)
                    for comp in v
                        push!(row, comp)
                    end
                else
                    push!(row, v)
                end
            end
            push!(row, meta_p.converged)
            push!(row, meta_p.niter)
            push!(row, meta_p.residual)
            println(io, join(row, ","))
            written_first_params = true
            if succ_p
                prev_solution_outvals = outvals_p
                prev_solution_vector = solvec_p
                prev_success = true
            end
            header_determined = true
        elseif header_probe == :first_batch
            # cache first probe_batchsize points
            cached = Vector{Tuple{NamedTuple, Dict{Symbol,Any}, Any, Any}}()
            # simple indexed iteration for first batch
            it = 0
            function inc_idx!(idx)
                N = length(idx)
                for i in N:-1:1
                    idx[i] += 1
                    if idx[i] <= length(axes[i].values)
                        return true
                    else
                        idx[i] = 1
                    end
                end
                return false
            end
            tmp_idx = fill(1, N)
            while it < probe_batchsize
                p = params_from_idx(tmp_idx)
                outvals_p, meta_p, succ_p, solvec_p = evaluate_point(p, nothing)
                push!(cached, (p, outvals_p, meta_p, solvec_p))
                it += 1
                cont = inc_idx!(tmp_idx)
                if !cont; break end
            end
            # build header from cached results (use first non-nothing values)
            header_cols = String[]
            append!(header_cols, param_names)
            for key in outputs
                found = false
                for (_p, outvals_p, _m, _s) in cached
                    v = get(outvals_p, key, nothing)
                    if v !== nothing
                        if isa(v, AbstractVector)
                            for i in 1:length(v)
                                push!(header_cols, string(key, "_", i))
                            end
                        else
                            push!(header_cols, string(key))
                        end
                        found = true
                        break
                    end
                end
                if !found
                    push!(header_cols, string(key))
                end
            end
            push!(header_cols, "solver_converged")
            push!(header_cols, "n_iter")
            push!(header_cols, "residual")

            io = _open_stream_writer(outfile, header_cols)
            # write cached rows
            for (p, outvals_p, meta_p, solvec_p) in cached
                row = Any[]
                for i in 1:length(param_names)
                    push!(row, getfield(p, param_syms[i]))
                end
                for key in outputs
                    v = get(outvals_p, key, NaN)
                    if isa(v, AbstractVector)
                        for comp in v
                            push!(row, comp)
                        end
                    else
                        push!(row, v)
                    end
                end
                push!(row, meta_p.converged)
                push!(row, meta_p.niter)
                push!(row, meta_p.residual)
                println(io, join(row, ","))
                if meta_p.converged
                    prev_solution_outvals = outvals_p
                    prev_solution_vector = solvec_p
                    prev_success = true
                end
            end
            header_determined = true
        else
            error("未知的 header_probe: $(header_probe)")
        end
    else
        # memory mode: we'll collect and write later; determine header when writing
    end

    function recurse(level)
        if level > N
            # 构建 params NamedTuple
            values = [axes[i].values[idx[i]] for i in 1:N]
            names = Tuple(Symbol.(param_names))
            params = NamedTuple{names}(Tuple(values))

            # 计算 initial guess
            g0 = nothing
            if initial_guess !== nothing
                try
                    g0 = initial_guess(params)
                catch
                    g0 = nothing
                end
            end
            # warm-start：若没有提供 initial_guess，且启用填充，则使用 prev_solution（仅在同一行内）
            # 我们重置 prev_solution 当 axes[1] 索引变化（见外层调用）
            if g0 === nothing && prev_solution_vector !== nothing
                g0 = prev_solution_vector
            end

            # evaluate current point using helper
            outvals, solver_meta, success, solvec = evaluate_point(params, prev_solution_vector)

            # 填充失败策略
            if !success
                push!(failed_points, (params=params, reason="solver_or_postprocess_failed"))
                if fill_failed_with_last && prev_success && prev_solution_outvals !== nothing
                    for (k,v) in prev_solution_outvals
                        outvals[Symbol(k)] = v
                    end
                else
                    for key in outputs
                        outvals[key] = NaN
                    end
                end
            else
                # 成功时记录 prev_solution
                prev_solution_outvals = outvals
                prev_solution_vector = solvec
                prev_success = true
            end

            # 更新全局最小 Omega
            if haskey(outvals, :Omega) && !isnan(outvals[:Omega]) && outvals[:Omega] < global_min.val
                global_min = (val = outvals[:Omega], params = params)
            end

            # 写入 CSV（流式）
            row = Any[]
            for i in 1:length(param_names)
                push!(row, values[i])
            end
            # outputs 列
            for key in outputs
                v = get(outvals, key, NaN)
                if isa(v, AbstractVector)
                    # 首次出现记录长度
                    if !haskey(vector_output_length, key)
                        vector_output_length[key] = length(v)
                        # 插入多列 header（简单处理：当前实现不动态调整已写入 header）
                        # 首版限制：如果遇到向量输出，视为定长；后续实现可更改
                    end
                    expected = vector_output_length[key]
                    if length(v) != expected
                        error("可变长度向量输出不被首版支持：输出 $(key) 在不同点的长度不一致")
                    end
                    for comp in v
                        push!(row, comp)
                    end
                else
                    push!(row, v)
                end
            end
            push!(row, solver_meta.converged)
            push!(row, solver_meta.niter)
            push!(row, solver_meta.residual)

            if mode == :stream
                # skip if we already wrote the first_params when probing
                if written_first_params && params == first_params
                    # this is the initial probed point; skip
                else
                    println(io, join(row, ","))
                end
            else
                # memory: TODO collect
            end

            return
        else
            for i in 1:length(axes[level].values)
                idx[level] = i
                # 在每次进入新行（level==1）时重置 prev_solution
                if level == 1
                    prev_solution = nothing
                    prev_success = false
                end
                recurse(level+1)
            end
        end
    end

    # 启动递归遍历
    recurse(1)

    # 关闭流式写入
    if mode == :stream && io !== nothing
        close(io)
        println("已将扫描结果写入 CSV: $(outfile) (stream mode)")
    end

    # 写 meta.json
    meta = Dict(
        "axes" => string([Dict("name"=>a.name, "n"=>length(a.values)) for a in axes]),
        "outputs" => string(outputs),
        "total_points" => total_points,
        "created_at" => string(Dates.now()),
        "tool_version" => "phase_scan_interface_v1"
    )
    try
        _write_meta(outfile, meta)
    catch
    end

    println("全局最小 Omega = $(global_min.val) at $(global_min.params)")
    return outfile
end

export ParameterAxis, generate_grid, wrap_physics_fn, scan_phase_space

end # module
