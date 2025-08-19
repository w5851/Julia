#!/usr/bin/env julia
"""
测试自动微分模块的基本功能

这个测试脚本验证新的AutomaticDifferentiation模块是否正确工作。
测试包括：梯度计算、Hessian计算、平衡条件检查等功能。

运行方法：
julia --project=. test_automatic_differentiation.jl
"""

using Test
using LinearAlgebra

# 激活项目环境并加载包
using Pkg; Pkg.activate(".")
using PNJLPhysicsSimulation

println("🧪 开始测试自动微分模块...")
println("=" ^ 60)

# 测试1: 基本梯度计算
println("测试1: 基本梯度计算")
println("-" ^ 30)

# 定义一个简单的二次函数 f(x) = x₁² + 2x₁x₂ + 3x₂²
test_function = (x) -> x[1]^2 + 2*x[1]*x[2] + 3*x[2]^2

# 测试点 [1.0, 2.0]
test_point = [1.0, 2.0]

# 解析解：∇f = [2x₁ + 2x₂, 2x₁ + 6x₂] = [2*1 + 2*2, 2*1 + 6*2] = [6, 14]
expected_gradient = [6.0, 14.0]

# 使用自动微分计算梯度
computed_gradient = AutomaticDifferentiation.compute_gradient(test_function, test_point)

println("测试点: $test_point")
println("期望梯度: $expected_gradient")
println("计算梯度: $computed_gradient")
println("误差: $(norm(computed_gradient - expected_gradient))")

@test isapprox(computed_gradient, expected_gradient, atol=1e-10)
println("✅ 梯度计算测试通过")

# 测试2: Hessian矩阵计算
println("\n测试2: Hessian矩阵计算")
println("-" ^ 30)

# 对于同样的函数，Hessian矩阵是：
# H = [∂²f/∂x₁², ∂²f/∂x₁∂x₂]   [2, 2]
#     [∂²f/∂x₂∂x₁, ∂²f/∂x₂²] = [2, 6]
expected_hessian = [2.0 2.0; 2.0 6.0]

computed_hessian = AutomaticDifferentiation.compute_hessian(test_function, test_point)

println("期望Hessian:")
display(expected_hessian)
println("计算Hessian:")
display(computed_hessian)
println("误差: $(norm(computed_hessian - expected_hessian))")

@test isapprox(computed_hessian, expected_hessian, atol=1e-10)
println("✅ Hessian矩阵计算测试通过")

# 测试3: 平衡条件检查
println("\n测试3: 平衡条件检查")
println("-" ^ 30)

# 在最小值点 x = [0, 0]，梯度应该为零
equilibrium_point = [0.0, 0.0]
equilibrium_result = AutomaticDifferentiation.check_equilibrium_conditions(
    test_function, equilibrium_point, tolerance=1e-8
)

println("平衡点: $equilibrium_point")
println("是否平衡: $(equilibrium_result.is_equilibrium)")
println("最大导数: $(equilibrium_result.max_derivative)")
println("梯度: $(equilibrium_result.gradients)")

@test equilibrium_result.is_equilibrium == true
@test equilibrium_result.max_derivative < 1e-10
println("✅ 平衡条件检查测试通过")

# 测试4: 非平衡点检查  
println("\n测试4: 非平衡点检查")
println("-" ^ 30)

non_equilibrium_point = [1.0, 1.0]
non_equilibrium_result = AutomaticDifferentiation.check_equilibrium_conditions(
    test_function, non_equilibrium_point, tolerance=1e-8
)

println("非平衡点: $non_equilibrium_point")
println("是否平衡: $(non_equilibrium_result.is_equilibrium)")
println("最大导数: $(non_equilibrium_result.max_derivative)")
println("最大导数位置: $(non_equilibrium_result.max_derivative_index)")

@test non_equilibrium_result.is_equilibrium == false
@test non_equilibrium_result.max_derivative > 1e-6
println("✅ 非平衡点检查测试通过")

# 测试5: 带附加信息的梯度计算
println("\n测试5: 带附加信息的梯度计算")
println("-" ^ 30)

variable_names = ["x₁", "x₂"]
gradient_info = AutomaticDifferentiation.compute_gradient(
    test_function, test_point, 
    return_info=true, variable_names=variable_names
)

println("变量名: $(gradient_info.variable_names)")
println("梯度: $(gradient_info.gradient)")
println("梯度范数: $(gradient_info.norm)")
println("计算时间: $(gradient_info.computation_time)秒")

@test length(gradient_info.gradient) == 2
@test gradient_info.norm > 0
@test gradient_info.computation_time >= 0
println("✅ 带附加信息的梯度计算测试通过")

# 测试6: 多参数函数测试（模拟物理函数）
println("\n测试6: 多参数函数测试")
println("-" ^ 30)

# 模拟一个简化的热力学势函数
# Ω(φ, T, μ) = φ² + T*φ - μ*φ³
physics_function = (variables, T, mu) -> begin
    phi = variables[1]
    return phi^2 + T*phi - mu*phi^3
end

# 测试参数
phi_test = [0.5]
T_test = 0.15
mu_test = 0.3

# 计算对φ的偏导数  
physics_gradient = AutomaticDifferentiation.compute_gradient(
    physics_function, phi_test, T_test, mu_test
)

# 解析解：∂Ω/∂φ = 2φ + T - 3μφ²
expected_physics_gradient = [2*0.5 + 0.15 - 3*0.3*0.5^2]

println("物理函数测试点: φ=$phi_test, T=$T_test, μ=$mu_test")
println("期望梯度: $expected_physics_gradient")  
println("计算梯度: $physics_gradient")
println("误差: $(norm(physics_gradient - expected_physics_gradient))")

@test isapprox(physics_gradient, expected_physics_gradient, atol=1e-10)
println("✅ 多参数函数测试通过")

# 测试7: 温度导数计算
println("\n测试7: 温度导数计算")
println("-" ^ 30)

# 对于上面的物理函数，∂Ω/∂T = φ
phi_fixed = [0.5]
dOmega_dT = AutomaticDifferentiation.compute_temperature_derivative(
    physics_function, phi_fixed, T_test, mu_test
)

expected_dOmega_dT = phi_fixed[1]  # = 0.5

println("固定变量: φ=$phi_fixed")
println("期望温度导数: $expected_dOmega_dT")
println("计算温度导数: $dOmega_dT")
println("误差: $(abs(dOmega_dT - expected_dOmega_dT))")

@test isapprox(dOmega_dT, expected_dOmega_dT, atol=1e-10)
println("✅ 温度导数计算测试通过")

println("\n" * "=" * 60)
println("🎉 所有自动微分模块测试通过！")
println("自动微分模块已成功从统一公共接口中分离出来，")
println("遵循了单一职责原则，专注于自动微分功能。")
println("=" * 60)
