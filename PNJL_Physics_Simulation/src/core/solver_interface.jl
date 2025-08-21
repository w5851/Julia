module SolverInterface
using NLsolve
export solve_equilibrium_equations

"""
    solve_equilibrium_equations(equation_system, initial_guess::Vector{T}, config)

说明：
- 只依赖 NLsolve；接受残差函数 `equation_system(x, config)`，在收敛时返回解向量 `Vector{T}`。
- 未收敛时抛出错误（调用方可通过 try/catch 处理）。
"""
function solve_equilibrium_equations(equation_system, initial_guess::Vector{T}, config::Any) where T<:Real
    residual_func!(F, x) = (F .= equation_system(x, config))
    result = nlsolve(residual_func!, initial_guess; autodiff=:forward)
    if !converged(result)
        error("求解未收敛: residual_norm=$(result.residual_norm), iterations=$(result.iterations)")
    end
    return result.zero
end

end # module SolverInterface

# 兼容绑定：便于旧代码直接调用函数名
const solve_equilibrium_equations = SolverInterface.solve_equilibrium_equations
