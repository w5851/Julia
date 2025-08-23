module EquationFactory

export create_equation_system_inplace, bind_params

"""
    create_equation_system_inplace(funcs_with_vars...; global_vars=nothing, param_names=nothing)

工厂函数：接受若干子方程 `(func, vars)`，其中 `vars` 是一组未知量名（Symbol 或 String）或数字索引。
可以通过 `param_names` 指定一组全局参数名（例如 `(:T,:mu)`）。
返回五元组 `(F_param!, names, param_names, pack, unpack)`：
- `F_param!(res, X, params)`：in-place 原型残差函数，接受额外参数对象 `params`（通常为 NamedTuple）。
- `names`：未知量顺序（Tuple{Symbol,...}）。
- `param_names`：参数名顺序（Tuple{Symbol,...}）。
- `pack(X)`：把向量转成 NamedTuple（字段顺序与 `names` 一致）。
- `unpack(nt)`：把 NamedTuple 转回向量（返回 Vector）。

示例：
```julia
F_param!, names, param_names, pack, unpack = create_equation_system_inplace((f, (:x,:y)), (g, (:y,:z)), (h, (:x,:z)); param_names=(:T,:mu))
```
"""
function create_equation_system_inplace(funcs_with_vars...; global_vars=nothing, param_names=nothing)
    # 验证输入形状
    if length(funcs_with_vars) == 0
        error("需要至少一个 (func, vars) 对")
    end

    # 1) 确定全局变量顺序
    if global_vars === nothing
        seen = Symbol[]
        for pair in funcs_with_vars
            if !(length(pair) >= 2)
                error("每个输入项必须是 (func, vars) 的元组")
            end
            vars = pair[2]
            for v in vars
                push!(seen, Symbol(v))
            end
        end
        names = unique(seen)
    else
        names = [Symbol(v) for v in global_vars]
    end

    name_to_idx = Dict{Symbol,Int}()
    for (i,n) in enumerate(names)
        name_to_idx[n] = i
    end

    # 参数名解析（要求每个子方程必须显式声明局部参数）
    # 每个输入项必须为 (func, vars, local_param_names)，local_param_names 可以为空 ()
    raw_parsed = Vector{Tuple{Any, Vector{Int}, Vector{Symbol}}}(undef, length(funcs_with_vars))
    for (i, pair) in enumerate(funcs_with_vars)
        if length(pair) != 3
            error("每个输入项必须是三元组 (func, vars, local_param_names)，例如 (f, (:x,:y), (:T,:mu))")
        end
        func = pair[1]
        vars = pair[2]
        idxs = Int[]
        for v in vars
            s = Symbol(v)
            if !haskey(name_to_idx, s)
                error("变量名 $(s) 未出现在全局变量列表中")
            end
            push!(idxs, name_to_idx[s])
        end
        # 局部参数名（允许为空元组）
        local_pnames = [Symbol(p) for p in pair[3]]
        raw_parsed[i] = (func, idxs, local_pnames)
    end

    # 构建全局 param_names（若用户未提供），否则使用用户提供的并校验局部声明
    if param_names === nothing
        # 按子方程首次出现顺序收集局部参数名
        seen_params = Symbol[]
        for (_func, _idxs, local_pnames) in raw_parsed
            for p in local_pnames
                if !(p in seen_params)
                    push!(seen_params, p)
                end
            end
        end
        param_names = seen_params
    else
        param_names = [Symbol(p) for p in param_names]
    end
    param_to_idx = Dict{Symbol,Int}(p => i for (i,p) in enumerate(param_names))

    # 为每个子方程构建局部参数索引（相对于全局 param_names）
    parsed = Vector{Tuple{Any, Vector{Int}, Vector{Int}}}(undef, length(raw_parsed))
    for (i, (func, idxs, local_pnames)) in enumerate(raw_parsed)
        local_idxs = Int[]
        for p in local_pnames
            if !haskey(param_to_idx, p)
                error("局部参数名 $(p) 未在全局 param_names 中找到；请在任一子方程中声明或在 param_names 参数中提供它们。")
            end
            push!(local_idxs, param_to_idx[p])
        end
        parsed[i] = (func, idxs, local_idxs)
    end

    # 3) in-place 原型残差：接受额外的 params（NamedTuple 或类似结构）
    # 简化：仅按局部参数展开调用（调用签名由用户通过三元组显式声明）
    function F_param!(res::AbstractVector, X::AbstractVector, params)
        @inbounds begin
            for (i,(func, idxs, local_pidxs)) in enumerate(parsed)
                # 将 X 中的未知量提取为位置参数
                # 简化：用户保证每个子方程有且只有 3 个未知量，因此直接展开三元组
                args_x = (X[idxs[1]], X[idxs[2]], X[idxs[3]])

                # 提取局部参数的值（按 global param_names 顺序的索引 local_pidxs）
                par_vals = isempty(local_pidxs) ? () : Tuple(getfield(params, param_names[k]) for k in local_pidxs)
                res[i] = func(args_x..., par_vals...)
            end
        end
        return res
    end

    # 为兼容求解器返回普通 F!(res, X) 的简单绑定器可以由调用方生成（见 bind_params）
    pack = X -> NamedTuple{Tuple(names)}(Tuple(X))
    unpack = nt -> collect(values(nt))

    return (F_param!, Tuple(names), Tuple(param_names), pack, unpack)
end

""" bind_params(F_param!, params) -> F!(res, X)

返回一个小闭包 `F!(res, X)`，在内部把当前 `params` 绑定并调用 `F_param!`。
通常 scanner 在每个点上构造一个短生命周期的 `params`（NamedTuple）并传给该闭包。
"""
function bind_params(F_param!, params)
    return (res, X) -> F_param!(res, X, params)
end
end # module
