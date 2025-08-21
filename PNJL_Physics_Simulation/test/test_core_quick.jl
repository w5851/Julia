# 轻量 core 单元测试（快速）
println("=== test_core_quick.jl ===")

# 加载项目环境
using Pkg
Pkg.instantiate()

# 引入 core 模块
include(joinpath(@__DIR__, "..", "src", "core", "math_utils.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "integration_interface.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "autodiff_interface.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "solver_interface.jl"))

using .MathUtils
using .IntegrationInterface
using .AutodiffInterface
using .SolverInterface

# 1. MathUtils 基本测试
println("1) MathUtils tests")
@assert MathUtils.safe_log(2.0) ≈ log(2.0)
@assert MathUtils.safe_log(-1.0) ≈ log(1e-16)
@assert MathUtils.safe_exp(1.0) ≈ exp(1.0)
@assert MathUtils.safe_exp(1000.0) ≈ exp(700.0)
@assert MathUtils.safe_sqrt(4.0) ≈ 2.0
@assert MathUtils.safe_division(1.0, 0.0) != Inf
println("   MathUtils OK")

# 2) IntegrationInterface gauleg test
println("2) IntegrationInterface gauleg test")
nodes, weights = IntegrationInterface.gauleg(0.0, 1.0, 5)
@assert length(nodes) == 5
@assert length(weights) == 5
println("   gauleg OK")

# 3) AutodiffInterface compute_gradient simple test
println("3) AutodiffInterface compute_gradient test")
# f(x) = x1^2 + x2^2
f = x -> x[1]^2 + x[2]^2
grad = AutodiffInterface.compute_gradient(f, [1.0, 2.0])
@assert length(grad) == 2
@assert grad[1] ≈ 2.0
@assert grad[2] ≈ 4.0
println("   compute_gradient OK")

# 4) SolverInterface quick test (solve simple system)
println("4) SolverInterface test")
# Solve x^2 - 4 = 0 (choose initial 1.0 and wrap to vector residual)
res = SolverInterface.solve_equilibrium_equations((x, c)->[x[1]^2 - 4.0], [1.0], nothing)
@assert length(res) == 1
@assert abs(res[1] - 2.0) < 1e-6
println("   solve_equilibrium_equations OK")

println("ALL core quick tests passed")
