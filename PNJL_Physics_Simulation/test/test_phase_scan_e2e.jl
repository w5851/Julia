# 集成测试：使用 EquationFactory -> SolverInterface 的 in-place API + pack
using Test

# Load factory and solver (relative paths)
include(joinpath(@__DIR__, "..", "src", "core", "equation_factory.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "solver_interface.jl"))
using .EquationFactory
using .SolverInterface

# 定义一个简单线性方程组（有解析解），依赖未知量 x, y, z 和参数 T, mu
# 方程：
#   f: x + y - T == 0
#   g: y - z + mu == 0
#   h: x - z == 0
# 解析解： z = (T + mu)/2, x = z, y = z - mu

f(x,y,T,mu) = x + y - T
g(y,z,T,mu) = y - z + mu
h(x,z,T,mu) = x - z

F_param!, names, param_names, pack, unpack = EquationFactory.create_equation_system_inplace((f, (:x,:y)), (g, (:y,:z)), (h, (:x,:z)); param_names=(:T,:mu))
@info "factory names" names param_names

params_list = ( (T=1.0, mu=0.0), (T=0.36, mu=0.2) )

# initial guess
x0 = [0.5, 0.5, 0.5]

for p in params_list
    # call solver with factory in-place function and pack
    result = SolverInterface.solve_equilibrium_equations(F_param!, x0, p; pack=pack)
    @test isa(result, NamedTuple)
    # expected analytic solution
    z_expected = (p.T + p.mu)/2
    x_expected = z_expected
    y_expected = z_expected - p.mu
    @test isapprox(result.x, x_expected; atol=1e-8, rtol=1e-8)
    @test isapprox(result.y, y_expected; atol=1e-8, rtol=1e-8)
    @test isapprox(result.z, z_expected; atol=1e-8, rtol=1e-8)
    @test haskey(result, :zero)
    @test length(result.zero) == 3
end

println("test_phase_scan_e2e: factory+solver pack integration OK")
