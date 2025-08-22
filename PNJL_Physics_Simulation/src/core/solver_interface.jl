module SolverInterface
using NLsolve
export solve_equilibrium_equations

"""
    solve_equilibrium_equations(equation_system, initial_guess::Vector{T}, config; pack=nothing)

说明：
- 只依赖 NLsolve；接受 in-place 残差函数 `equation_system(res, x, config)` 并返回解（默认返回向量）。
- 可选关键字参数 `pack`：若提供，应为一个函数 `pack(::AbstractVector) -> NamedTuple`，求解完成后会调用 `pack(result.zero)` 并返回包含命名字段的 NamedTuple（并保留 `:zero` 字段）。
- 未收敛时抛出错误（调用方可通过 try/catch 处理）。
"""
function solve_equilibrium_equations(equation_system, initial_guess::Vector{T}, config::Any; pack=nothing) where T<:Real

    # residual_func! expects equation_system to be in-place: equation_system(res, x, config)
    function residual_func!(F, x)
        equation_system(F, x, config)
        return F
    end

    result = nlsolve(residual_func!, initial_guess; autodiff=:forward)
    if !converged(result)
        error("求解未收敛: residual_norm=$(result.residual_norm), iterations=$(result.iterations)")
    end
    if pack !== nothing
        nt = pack(result.zero)
        return merge(nt, (zero = result.zero,))
    end

    return result.zero
end

end # module SolverInterface

# 兼容绑定：便于旧代码直接调用函数名
const solve_equilibrium_equations = SolverInterface.solve_equilibrium_equations
