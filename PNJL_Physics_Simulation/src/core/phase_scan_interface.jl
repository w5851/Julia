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
function scan_phase_space(f::Function, axes::Vector{ParameterAxis}; outfile::AbstractString="output/phase_scan.csv", batchsize::Int=256)
    # 评估并收集结果（小规模实现：收集到内存后写 CSV）
    results = Vector{NamedTuple{(:T, :mu, :Omega),Tuple{Float64,Float64,Float64}}}()
    npoints = 0
    for params in generate_grid(axes)
        T = haskey(params, :T) ? params.T : (haskey(params, :temperature) ? params.temperature : NaN)
        mu = haskey(params, :mu) ? params.mu : (haskey(params, :chemical_potential) ? params.chemical_potential : NaN)
        val = try f(params) catch e; NaN end
        push!(results, (T, mu, float(val)))
        npoints += 1
    end

    # 写入 CSV（依赖 CSV.jl, DataFrames.jl）
    # 尝试使用 CSV + DataFrames 写入；若不可用则使用轻量级后备实现
    function _write_csv_fallback(outfile, results)
        outdir = dirname(outfile)
        if !isdir(outdir)
            try
                mkpath(outdir)
            catch
            end
        end
        open(outfile, "w") do io
            println(io, "T,mu,Omega")
            for r in results
                println(io, "$(r.T),$(r.mu),$(r.Omega)")
            end
        end
        println("已将扫描结果写入 CSV (fallback): $(outfile)")
    end

    try
        @eval begin
            using CSV, DataFrames
        end
        df = DataFrame(T = [r.T for r in results], mu = [r.mu for r in results], Omega = [r.Omega for r in results])
        # 确保输出目录存在
        outdir = dirname(outfile)
        if !isdir(outdir)
            try
                mkpath(outdir)
            catch
            end
        end
        CSV.write(outfile, df)
        println("已将扫描结果写入 CSV: $(outfile)")
    catch e
        try
            _write_csv_fallback(outfile, results)
        catch e2
            println("写入 CSV 失败：$(e) — 备份写入也失败: $(e2)")
        end
    end
    return outfile
end

export ParameterAxis, generate_grid, wrap_physics_fn, scan_phase_space

end # module
