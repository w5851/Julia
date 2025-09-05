#!/usr/bin/env julia

"""
迭代初解策略性能对比测试

比较固定初解策略和迭代初解策略在温度扫描中的性能差异
"""

using BenchmarkTools
include("../../src/Function_Gas_Liquid.jl")

println("="^60)
println("迭代初解策略性能对比测试")
println("="^60)

# 设置基本参数
nodes = get_nodes(128)  # 减少节点数以加快测试
gsigma = 1.25
gdelta = 0.01
fs = 10.329
fo = 5.423
fr = 3.15
fd = 2.5
b = 0.00692
c = -0.0048
couplings = [fs, fo, fr, fd, b, c]

# 测试参数
μ_B_fixed = 1001.0/hc
T_min = 30.0/hc
T_max = 50.0/hc
T_step = 5.0/hc  # 较大步长用于快速测试
mu_p = μ_B_fixed / 2.0
mu_n = μ_B_fixed / 2.0
x0 = [gsigma, gdelta, mu_p, mu_n]

println("测试参数:")
println("- 化学势: $(μ_B_fixed*hc) MeV")
println("- 温度范围: $(T_min*hc) - $(T_max*hc) MeV")
println("- 温度步长: $(T_step*hc) MeV")
println("- 积分节点数: 128")

println("\n1. 测试固定初解策略:")
println("-"^40)
@time begin
    temp_array_fixed, kappa31_fixed, kappa42_fixed, results_fixed = 
        calculate_fluctuation_ratios_vs_temperature_advanced(μ_B_fixed, T_min, T_max, x0, nodes, couplings,
                                                            T_step=T_step, save_results=false,
                                                            use_iterative_guess=false)
end

println("\n2. 测试迭代初解策略:")
println("-"^40)
@time begin
    temp_array_iter, kappa31_iter, kappa42_iter, results_iter, solutions = 
        calculate_fluctuation_ratios_vs_temperature_advanced(μ_B_fixed, T_min, T_max, x0, nodes, couplings,
                                                            T_step=T_step, save_results=false,
                                                            use_iterative_guess=true, extrapolation_weight=0.3,
                                                            return_solution_history=true)
end

println("\n3. 结果对比:")
println("-"^40)
println("温度点数: $(length(temp_array_fixed))")

println("\n收敛解演化（迭代策略）:")
for i in eachindex(solutions)
    sol = solutions[i]
    T_val = temp_array_iter[i] * hc
    println("T = $(round(T_val, digits=1)) MeV: [gσ=$(round(sol[1], digits=4)), gδ=$(round(sol[2], digits=4)), μ_p=$(round(sol[3]*hc, digits=2)), μ_n=$(round(sol[4]*hc, digits=2))]")
end

println("\nκ₃/κ₁ 比较:")
println("T(MeV)  固定初解      迭代初解      差异")
for i in eachindex(temp_array_fixed)
    T_val = temp_array_fixed[i] * hc
    diff = abs(kappa31_fixed[i] - kappa31_iter[i])
    println("$(rpad(round(T_val, digits=1), 7)) $(rpad(round(kappa31_fixed[i], digits=6), 12)) $(rpad(round(kappa31_iter[i], digits=6), 12)) $(round(diff, digits=8))")
end

println("\nκ₄/κ₂ 比较:")
println("T(MeV)  固定初解      迭代初解      差异")
for i in eachindex(temp_array_fixed)
    T_val = temp_array_fixed[i] * hc
    diff = abs(kappa42_fixed[i] - kappa42_iter[i])
    println("$(rpad(round(T_val, digits=1), 7)) $(rpad(round(kappa42_fixed[i], digits=6), 12)) $(rpad(round(kappa42_iter[i], digits=6), 12)) $(round(diff, digits=8))")
end

# 计算总体差异
total_diff_31 = sum(abs.(kappa31_fixed .- kappa31_iter))
total_diff_42 = sum(abs.(kappa42_fixed .- kappa42_iter))

println("\n4. 总体评估:")
println("-"^40)
println("κ₃/κ₁ 总差异: $(round(total_diff_31, digits=8))")
println("κ₄/κ₂ 总差异: $(round(total_diff_42, digits=8))")

if total_diff_31 < 1e-10 && total_diff_42 < 1e-10
    println("✅ 两种策略得到几乎相同的结果")
else
    println("⚠️  两种策略存在可观察的差异")
end

println("\n" * "="^60)
println("性能测试完成!")
println("="^60)
