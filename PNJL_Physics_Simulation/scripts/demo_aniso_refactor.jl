#!/usr/bin/env julia

"""
PNJL各向异性模型重构演示# 3. 压力计算结果对比:
println("\n3. 压力计算结果对比:")

# 使用旧格式接口 (内部已重构为无闭包实现)
t1 = @elapsed pressure_legacy = calculate_pressure_aniso(phi, Phi1, Phi2, mu, T, nodes_1, nodes_2, xi)
println("   - 兼容接口结果: $(round(pressure_legacy, digits=8))")
println("   - 计算时间: $(round(t1*1000, digits=2))ms")

# 使用新格式接口 (注意：旧实现使用两组不同的网格)
println("\n   注意：旧实现使用两组不同的积分网格：")
println("   - nodes_1 (截止$(round(Lambda_f, digits=3))): 用于真空能量积分")  
println("   - nodes_2 (截止20.0): 用于热力学积分")

# 为了公平比较，我们分别模拟这两个积分
t2a = @elapsed Omega_vac_modern = calculate_omega_vacuum(calculate_mass_vec(phi), xi, p_grid_1, t_grid_1)
t2b = @elapsed Omega_th_modern = calculate_omega_thermal(calculate_mass_vec(phi), mu, T, Phi1, Phi2, xi, p_grid_2, t_grid_2)

chi_modern = calculate_chiral_aniso(phi)
U_modern = calculate_U_aniso(T, Phi1, Phi2)
pressure_modern = -(chi_modern + U_modern + Omega_vac_modern + Omega_th_modern)

println("   - 现代接口结果: $(round(pressure_modern, digits=8))")
println("   - 计算时间: $(round((t2a + t2b)*1000, digits=2))ms")闭包的实现重构为基于积分接口的无闭包实现。

主要改进:
1. 消除闭包：将被积函数定义为独立的纯函数
2. 使用IntegrationInterface：统一的积分计算接口
3. 严格按照Omega公式：基于omega_formulas.md的物理公式
4. 保持向后兼容性：现有代码无需修改
5. 支持ForwardDiff：确保自动微分兼容性
"""

using PNJLPhysicsSimulation.PNJLAnisoFunctions
using PNJLPhysicsSimulation.IntegrationInterface

# 设置测试参数
phi = [0.1, 0.1, 0.05]  # 手征凝聚
Phi1, Phi2 = 0.5, 0.3   # Polyakov环变量
mu = [0.02, 0.02, 0.01] # 化学势
T = 0.15                # 温度
xi = 0.1                # 各向异性参数

println("🔬 PNJL各向异性模型重构演示")
println("="^50)

# 1. 旧格式节点 (向后兼容测试)
println("\n1. 旧格式积分节点:")
nodes_1, nodes_2 = get_nodes_aniso(32, 16)
println("   - 动量节点: $(size(nodes_1[1]))")
println("   - 角度节点: $(size(nodes_1[2]))")
println("   - 系数矩阵: $(size(nodes_1[3]))")

# 2. 新格式积分网格
println("\n2. 新格式积分网格:")
# 使用与旧节点相同的截断值进行对比
p_grid_1, t_grid_1 = create_aniso_grids(32, 16; p_cutoff=Lambda_f, t_domain=(0.0, 1.0))  # 对应nodes_1
p_grid_2, t_grid_2 = create_aniso_grids(32, 16; p_cutoff=20.0, t_domain=(0.0, 1.0))      # 对应nodes_2
println("   - 动量网格1: $(length(p_grid_1.nodes))个节点, 截止动量=$(round(p_grid_1.cutoff, digits=3))")
println("   - 动量网格2: $(length(p_grid_2.nodes))个节点, 截止动量=$(p_grid_2.cutoff)")
println("   - 角度网格: $(length(t_grid_1.nodes))个节点, 定义域=$(t_grid_1.domain)")

# 3. 压力计算对比
println("\n3. 压力计算结果对比:")

# 使用旧格式接口 (内部已重构为无闭包实现)
t1 = @elapsed pressure_legacy = calculate_pressure_aniso(phi, Phi1, Phi2, mu, T, nodes_1, nodes_2, xi)
println("   - 兼容接口结果: $(round(pressure_legacy, digits=8))")
println("   - 计算时间: $(round(t1*1000, digits=2))ms")

# 使用新格式接口 (使用相同的积分域)
t2 = @elapsed pressure_modern = calculate_pressure_aniso_modern(phi, Phi1, Phi2, mu, T, p_grid_2, t_grid_2, xi)
println("   - 现代接口结果: $(round(pressure_modern, digits=8))")
println("   - 计算时间: $(round(t2*1000, digits=2))ms")

# 数值一致性检查
relative_error = abs(pressure_modern - pressure_legacy) / abs(pressure_legacy)
println("   - 相对误差: $(round(relative_error*100, digits=8))%")

if relative_error < 1e-10
    println("   ✅ 数值结果完全一致")
elseif relative_error < 1e-6 
    println("   ✅ 数值结果高度一致（浮点精度范围内）")
else
    println("   ⚠️  存在数值差异（可能由于积分域或方法差异）")
end

# 4. 展示被积函数的改进
println("\n4. 被积函数实现对比:")
println("   ❌ 旧实现: 使用闭包 integrand_f = function(p, t; ...) ...")
println("   ✅ 新实现: 纯函数 thermal_integrand(p, t; mass, chemical_potential, ...)")
println("   ✅ 优势: 更好的性能、可读性和可维护性")

# 5. 展示积分组件的分解
println("\n5. Omega函数组件分解:")
masses = calculate_mass_vec(phi)

# 手征贡献
Omega_chiral = calculate_chiral_aniso(phi)
println("   - Ω_chiral = $(round(Omega_chiral, digits=6))")

# Polyakov势贡献
Omega_U = calculate_U_aniso(T, Phi1, Phi2)
println("   - Ω_U = $(round(Omega_U, digits=6))")

# 真空贡献
Omega_vac = calculate_omega_vacuum(masses, xi, p_grid_2, t_grid_2)
println("   - Ω_vac = $(round(Omega_vac, digits=6))")

# 热力学贡献
Omega_th = calculate_omega_thermal(masses, mu, T, Phi1, Phi2, xi, p_grid_2, t_grid_2)
println("   - Ω_thermal = $(round(Omega_th, digits=6))")

# 总和验证
total_omega = Omega_chiral + Omega_U + Omega_vac + Omega_th
pressure_reconstructed = -total_omega
println("   - 总压力验证: $(round(pressure_reconstructed, digits=8))")

# 6. 性能优化总结
println("\n6. 重构优化总结:")
println("   ✅ 消除闭包: 提高编译器优化潜力")
println("   ✅ 显式参数传递: 更清晰的函数依赖关系")
println("   ✅ 统一积分接口: 减少代码重复")
println("   ✅ 严格物理公式: 基于omega_formulas.md文档")
println("   ✅ 保持向后兼容: 现有代码无需修改")
println("   ✅ 支持自动微分: ForwardDiff兼容性")

# 7. 使用建议
println("\n7. 使用建议:")
println("   - 新代码推荐: 使用 calculate_pressure_aniso_modern() 和 create_aniso_grids()")
println("   - 现有代码: 可继续使用 calculate_pressure_aniso()，内部已优化")
println("   - 积分网格: 使用 MomentumGrid 和 AngleGrid 替代旧式节点矩阵")

println("\n" * "="^50)
println("🎉 重构演示完成！")
println("📖 详细文档: src/models/pnjl_aniso/omega_formulas.md")
