#!/usr/bin/env julia
using Dates
"""
气液相变程序 - 温度扫描涨落比值计算演示

该脚本演示如何使用新添加的函数 calculate_fluctuation_ratios_vs_temperature
来固定重子化学势μ_B，改变温度T，计算不同温度下的κ₃/κ₁和κ₄/κ₂比值。

作者: 基于 Function_Gas_Liquid.jl
日期: $(Dates.format(now(), "yyyy-mm-dd"))
"""


include("../../src/Function_Gas_Liquid.jl")

println("="^60)
println("气液相变程序 - 温度扫描涨落比值计算演示")
println("="^60)

# 设置基本参数
println("1. 设置基本参数...")
nodes = get_nodes(256)
gsigma = 1.25
gdelta = 0.01
fs = 10.329
fo = 5.423
fr = 3.15
fd = 2.5
b = 0.00692
c = -0.0048
couplings = [fs, fo, fr, fd, b, c]

# 设置温度扫描参数
μ_B_fixed = 1001.0/hc  # 固定重子化学势 1001 MeV
T_min = 20.0/hc        # 最小温度 20 MeV
T_max = 50.0/hc        # 最大温度 50 MeV
T_step = 2.0/hc        # 温度步长 2 MeV (为了演示，使用较大步长)

println("固定重子化学势: μ_B = $(μ_B_fixed*hc) MeV")
println("温度扫描范围: $(T_min*hc) - $(T_max*hc) MeV")
println("温度步长: $(T_step*hc) MeV")

# 初始猜测值
mu_p = μ_B_fixed / 2.0
mu_n = μ_B_fixed / 2.0
x0 = [gsigma, gdelta, mu_p, mu_n]

println("\n2. 开始温度扫描计算...")

# 执行温度扫描
temperature_array, kappa3_over_kappa1, kappa4_over_kappa2, results_matrix = 
    calculate_fluctuation_ratios_vs_temperature(μ_B_fixed, T_min, T_max, x0, nodes, couplings,
                                               T_step=T_step, save_results=true,
                                               output_file=joinpath(@__DIR__, "..", "..", "output", "demo_fluctuation_ratios_vs_T.csv"))

println("\n3. 结果分析:")
println("-"^40)

# 转换温度到MeV单位显示
T_MeV_array = temperature_array .* hc

println("温度范围: $(minimum(T_MeV_array)) - $(maximum(T_MeV_array)) MeV")
println("计算点数: $(length(T_MeV_array))")

# 显示所有结果
println("\n完整结果表:")
println("T(MeV)    κ₃/κ₁        κ₄/κ₂")
println("-"^40)
for i in eachindex(temperature_array)
    T_MeV = T_MeV_array[i]
    ratio_31 = kappa3_over_kappa1[i]
    ratio_42 = kappa4_over_kappa2[i]
    if isfinite(ratio_31) && isfinite(ratio_42)
        println("$(rpad(round(T_MeV, digits=1), 8)) $(rpad(round(ratio_31, digits=6), 12)) $(round(ratio_42, digits=6))")
    else
        println("$(rpad(round(T_MeV, digits=1), 8)) $(rpad("NaN", 12)) NaN")
    end
end

# 统计分析
finite_31 = filter(isfinite, kappa3_over_kappa1)
finite_42 = filter(isfinite, kappa4_over_kappa2)

println("\n4. 统计信息:")
println("-"^40)
println("κ₃/κ₁ 有效值数量: $(length(finite_31))/$(length(kappa3_over_kappa1))")
if length(finite_31) > 0
    println("κ₃/κ₁ 范围: $(round(minimum(finite_31), digits=6)) 到 $(round(maximum(finite_31), digits=6))")
end

println("κ₄/κ₂ 有效值数量: $(length(finite_42))/$(length(kappa4_over_kappa2))")
if length(finite_42) > 0
    println("κ₄/κ₂ 范围: $(round(minimum(finite_42), digits=6)) 到 $(round(maximum(finite_42), digits=6))")
end

println("\n5. 结果矩阵格式:")
println("-"^40)
println("维度: $(size(results_matrix))")
println("列含义: [T(MeV), κ₃/κ₁, κ₄/κ₂]")

println("\n6. 输出文件:")
println("-"^40)
output_file = joinpath(@__DIR__, "..", "..", "output", "demo_fluctuation_ratios_vs_T.csv")
println("结果已保存到: $output_file")

println("\n" * "="^60)
println("演示完成!")
println("="^60)
