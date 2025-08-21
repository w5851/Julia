# 轻量 smoke 测试：验证 equation_solver 最简实现

# 此测试构造简单线性方程组：
#   x1 - 2 = 0
#   x2 - 3 = 0
# 期望解：[2.0, 3.0]

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using NLsolve
using LinearAlgebra
include(joinpath(@__DIR__, "..", "src", "core", "equation_solver.jl"))

# 简单残差函数（不使用复杂 config）
function linear_residual(x, config)
    return [x[1] - 2.0, x[2] - 3.0]
end

initial_guess = [0.0, 0.0]

all_passed = true

sol = solve_equilibrium_equations(linear_residual, initial_guess, nothing)
println("求解结果: ", sol)
if norm(sol - [2.0, 3.0]) < 1e-8
    println("SMOKE TEST PASS")
else
    println("SMOKE TEST FAIL")
    all_passed = false
end

# 额外测试：残差函数使用 config 参数
println("\n运行带 config 参数的测试...")
function linear_residual_with_config(x, cfg)
    # cfg 可以包含偏置或系数
    a = cfg[:a]
    b = cfg[:b]
    return [x[1] - a, x[2] - b]
end

cfg = Dict(:a => 2.0, :b => 3.0)
initial_guess2 = [1.0, 1.0]
sol2 = solve_equilibrium_equations(linear_residual_with_config, initial_guess2, cfg)
println("带 config 的求解结果: ", sol2)
if norm(sol2 - [2.0, 3.0]) < 1e-8
    println("CONFIG SMOKE TEST PASS")
else
    println("CONFIG SMOKE TEST FAIL")
    all_passed = false
end

if all_passed
    exit(0)
else
    exit(1)
end
