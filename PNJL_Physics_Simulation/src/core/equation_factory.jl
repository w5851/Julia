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

    # 参数名解析（全局参数）
    if param_names === nothing
        param_names = Symbol[]
    else
        param_names = [Symbol(p) for p in param_names]
    end
    param_to_idx = Dict{Symbol,Int}(p => i for (i,p) in enumerate(param_names))

    # 2) 解析每个子方程的索引列表（未知量索引）
    parsed = Vector{Tuple{Any, Vector{Int}}}(undef, length(funcs_with_vars))
    for (i, pair) in enumerate(funcs_with_vars)
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
        parsed[i] = (func, idxs)
    end

    # 3) in-place 原型残差：接受额外的 params（NamedTuple 或类似结构）
    function F_param!(res::AbstractVector, X::AbstractVector, params)
        @inbounds begin
            for (i,(func, idxs)) in enumerate(parsed)
                # 将 X 中的未知量提取为位置参数
                # 以常见的小参数长度展开减少分配
                if length(idxs) == 1
                    args_x = (X[idxs[1]],)
                elseif length(idxs) == 2
                    args_x = (X[idxs[1]], X[idxs[2]])
                elseif length(idxs) == 3
                    args_x = (X[idxs[1]], X[idxs[2]], X[idxs[3]])
                else
                    args_x = Tuple(map(j -> X[j], idxs))
                end

                if length(param_names) == 0
                    res[i] = func(args_x...)
                else
                    # 将 params 按 param_names 顺序展开为位置参数
                    par_vals = Tuple(getfield(params, n) for n in param_names)
                    res[i] = func(args_x..., par_vals...)
                end
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
