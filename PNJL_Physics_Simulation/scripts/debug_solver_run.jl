# 调试脚本：直接调用 EquationFactory 与 SolverInterface，打印 nlsolve 详细信息
include(joinpath(@__DIR__, "..", "src", "core", "equation_factory.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "solver_interface.jl"))
using .EquationFactory
using .SolverInterface

f(x, y, T, mu) = x - T*0.1
g(y, x, T, mu) = y - mu*0.01

F_param!, names, param_names, pack, unpack = create_equation_system_inplace(
    (f, (:x, :y), (:T,)),
    (g, (:y, :x), (:mu,))
)

println("names=", names, " param_names=", param_names)
params = (T=100.0, mu=10.0)
X0 = [0.0, 0.0]

println("Calling solve_equilibrium_equations with initial guess ", X0)
try
    result = SolverInterface.solve_equilibrium_equations(F_param!, X0, params; pack=pack)
    println("solver returned: ", result)
catch e
    println("solver raised exception: ", e)
    # try to run nlsolve directly to inspect result
    try
        function residual_func!(F, x)
            F_param!(F, x, params)
            return F
        end
        using NLsolve
        res = nlsolve(residual_func!, X0; autodiff=:forward)
        println("nlsolve converged=", converged(res))
        println("residual_norm=", res.residual_norm)
        println("iterations=", res.iterations)
        println("zero type=", typeof(res.zero), " values=", res.zero)
    catch e2
        println("direct nlsolve raised: ", e2)
    end
end
