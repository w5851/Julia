#!/usr/bin/env julia

"""
气液相变程序 - Function_Gas_Liquid.jl 示例用法

这个脚本展示了如何使用 Function_Gas_Liquid.jl 中的各种函数
包括：
1. 单点计算示例
2. 批量化学势扫描
3. 固定μ_B的温度扫描

作者: 基于 Function_Gas_Liquid.jl
日期: 2025-09-05
"""

# 包含主要函数文件
include("../../src/Function_Gas_Liquid.jl")

println("="^60)
println("气液相变程序 - Function_Gas_Liquid.jl 示例用法")
println("="^60)

#=============================================================================
1. 基本参数设置
=============================================================================#
println("\n1. 设置基本参数...")

# 设置基本参数
nodes = get_nodes(256)
T = 50.0/hc
gsigma = 1.25
gdelta = 0.01
fs = 10.329
fo = 5.423
fr = 3.15
fd = 2.5
b = 0.00692
c = -0.0048
couplings = [fs, fo, fr, fd, b, c]

println("积分节点数: 256")
println("温度: $(T*hc) MeV")
println("耦合常数设置完成")

#=============================================================================
2. 单点计算示例
=============================================================================#
println("\n" * "="^50)
println("2. 单点计算示例")
println("="^50)

# 单点计算示例
μ_B = 1001.0/hc
params = [T, μ_B]
mu_p = μ_B / 2.0
mu_n = μ_B / 2.0
x0 = [gsigma, gdelta, mu_p, mu_n]

println("重子化学势: $(μ_B*hc) MeV")
println("初始猜测值: [gσ=$(gsigma), gδ=$(gdelta), μ_p=$(mu_p*hc) MeV, μ_n=$(mu_n*hc) MeV]")

# 方法1: 使用通用导数计算函数
println("\n方法1: 使用通用导数计算函数")
pressure, derivatives = calculate_pressure_derivatives(μ_B, T, x0, nodes, couplings)
println("压强: ", pressure)
println("一阶导数: ", derivatives[1])
println("二阶导数: ", derivatives[2])
println("三阶导数: ", derivatives[3])
println("四阶导数: ", derivatives[4])

# 方法2: 使用高效导数计算函数
println("\n方法2: 使用高效导数计算函数")
pressure, dpre_dmu1, dpre_dmu2, dpre_dmu3, dpre_dmu4 = 
    calculate_pressure_derivatives_efficient(μ_B, T, x0, nodes, couplings)
println("压强: ", pressure)
println("一阶导数: ", dpre_dmu1)
println("二阶导数: ", dpre_dmu2)
println("三阶导数: ", dpre_dmu3)
println("四阶导数: ", dpre_dmu4)

# 方法3: 计算热力学涨落
println("\n方法3: 计算热力学涨落")
kappa1, kappa2, kappa3, kappa4, fluctuation_ratios = 
    calculate_thermodynamic_fluctuations(μ_B, T, x0, nodes, couplings)
    
println("累积量:")
println("κ₁ = ", kappa1)
println("κ₂ = ", kappa2) 
println("κ₃ = ", kappa3)
println("κ₄ = ", kappa4)
println("涨落比值: κ₂/κ₁ = ", fluctuation_ratios[1])
println("涨落比值: κ₃/κ₂ = ", fluctuation_ratios[2])
println("涨落比值: κ₄/κ₂ = ", fluctuation_ratios[3])

#=============================================================================
3. 批量化学势扫描示例
=============================================================================#
println("\n" * "="^50)
println("3. 批量化学势扫描示例")
println("="^50)

# 批量计算示例
μ_B_range = 1001.0/hc:-50.0/hc:800.0/hc  # 从1001 MeV到800 MeV，步长50 MeV（演示用）
println("化学势扫描范围: $(maximum(μ_B_range)*hc) - $(minimum(μ_B_range)*hc) MeV")
println("化学势步长: $(abs(μ_B_range[2]-μ_B_range[1])*hc) MeV")
println("计算点数: $(length(μ_B_range))")

results = calculate_derivatives_batch(μ_B_range, T, x0, nodes, couplings, 
                                    save_results=true, 
                                    output_file=joinpath(@__DIR__, "..", "..", "output", "pressure_derivatives_example.csv"))

println("批量计算完成，结果已保存为CSV格式")

#=============================================================================
4. 固定μ_B扫描温度计算涨落比值
=============================================================================#
println("\n" * "="^50)
println("4. 固定μ_B扫描温度计算涨落比值示例")
println("="^50)

# 设置温度扫描参数
μ_B_fixed = 1001.0/hc  # 固定重子化学势
T_min = 20.0/hc        # 最小温度 20 MeV
T_max = 80.0/hc        # 最大温度 80 MeV
T_step = 5.0/hc        # 温度步长 5 MeV（演示用）

println("固定重子化学势: $(μ_B_fixed*hc) MeV")
println("温度扫描范围: $(T_min*hc) - $(T_max*hc) MeV")
println("温度步长: $(T_step*hc) MeV")

# 执行温度扫描
temperature_array, kappa3_over_kappa1, kappa4_over_kappa2, results_matrix = 
    calculate_fluctuation_ratios_vs_temperature(μ_B_fixed, T_min, T_max, x0, nodes, couplings,
                                               T_step=T_step, save_results=true,
                                               output_file=joinpath(@__DIR__, "..", "..", "output", "fluctuation_ratios_vs_T_example.csv"))

# 显示部分结果
println("\n前5个温度点的结果:")
println("T(MeV)    κ₃/κ₁         κ₄/κ₂")
println("-"^40)
for i in 1:min(5, length(temperature_array))
    T_MeV = temperature_array[i] * hc
    ratio_31 = kappa3_over_kappa1[i]
    ratio_42 = kappa4_over_kappa2[i]
    if isfinite(ratio_31) && isfinite(ratio_42)
        println("$(rpad(round(T_MeV, digits=1), 8))  $(rpad(round(ratio_31, digits=6), 12))  $(round(ratio_42, digits=6))")
    else
        println("$(rpad(round(T_MeV, digits=1), 8))  $(rpad("NaN", 12))  NaN")
    end
end

println("\n结果矩阵格式 [T(MeV), κ₃/κ₁, κ₄/κ₂]:")
println("维度: $(size(results_matrix))")
println("前3行示例:")
for i in 1:min(3, size(results_matrix, 1))
    println("$(round(results_matrix[i,1], digits=1))  $(round(results_matrix[i,2], digits=6))  $(round(results_matrix[i,3], digits=6))")
end

#=============================================================================
5. 总结
=============================================================================#
println("\n" * "="^50)
println("5. 示例演示完成")
println("="^50)

println("输出文件:")
println("- 压强导数批量计算: ../../output/pressure_derivatives_example.csv")
println("- 涨落比值温度扫描: ../../output/fluctuation_ratios_vs_T_example.csv")
println("\n所有输出文件均为CSV格式，便于数据分析和绘图")

println("\n" * "="^60)
println("演示完成!")
println("="^60)
