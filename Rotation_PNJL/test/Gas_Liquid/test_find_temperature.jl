# test_find_temperature.jl
# 测试温度反向查找功能
# 测试Advanced_FindTforDiff.jl模块中的温度查找功能

using Printf
using Dates

# 引入温度查找模块
include("../../src/Gas_Liquid/Advanced_FindTforDiff.jl")

println("="^80)
println("温度反向查找功能测试")
println("测试时间: $(Dates.now())")
println("="^80)

# 测试参数设置
μ_B_test = 697.0 / hc  # 697 MeV 重子化学势
T_min = 20.0 / hc      # 20 MeV 
T_max = 200.0 / hc     # 200 MeV
T_step = 2.0 / hc      # 2 MeV 步长（加快测试速度）

# 目标κ比值（示例值）
target_kappa3_kappa1 = 0.5  # 目标 κ₃/κ₁ 值
target_kappa4_kappa2 = 1.2  # 目标 κ₄/κ₂ 值

println("\n测试配置:")
println("μ_B = $(μ_B_test * hc) MeV")
println("T范围: $(T_min * hc) - $(T_max * hc) MeV")
println("T步长: $(T_step * hc) MeV")
println("目标 κ₃/κ₁ = $target_kappa3_kappa1")
println("目标 κ₄/κ₂ = $target_kappa4_kappa2")

# 执行温度查找测试
println("\n开始温度反向查找测试...")
try
    T_kappa3_kappa1, T_kappa4_kappa2 = find_temperature_for_kappa_ratios(
        target_kappa3_kappa1, target_kappa4_kappa2, μ_B_test, 
        T_min, T_max, T_step;
        gsigma=1.25, gdelta=0.01,
        fs=17.28476, fo=11.66174, fr=0.89363, fd=0.0,
        b=0.00210, c=-0.00297, n_nodes=128,  # 减少积分点以加快测试
        verbose=true
    )
    
    println("\n" * "="^80)
    println("温度反向查找测试完成!")
    println("="^80)
    
    # 输出结果
    println("\n最终结果:")
    if !isnan(T_kappa3_kappa1)
        @printf("κ₃/κ₁ = %.3f 对应温度: T = %.2f MeV\n", target_kappa3_kappa1, T_kappa3_kappa1 * hc)
    else
        println("κ₃/κ₁ = $target_kappa3_kappa1: 未找到对应温度")
    end
    
    if !isnan(T_kappa4_kappa2)
        @printf("κ₄/κ₂ = %.3f 对应温度: T = %.2f MeV\n", target_kappa4_kappa2, T_kappa4_kappa2 * hc)
    else
        println("κ₄/κ₂ = $target_kappa4_kappa2: 未找到对应温度")
    end
    
    # 计算温度差异
    if !isnan(T_kappa3_kappa1) && !isnan(T_kappa4_kappa2)
        T_diff = abs(T_kappa3_kappa1 - T_kappa4_kappa2) * hc
        @printf("\n两个κ比值对应温度差异: %.2f MeV\n", T_diff)
        
        if T_diff < 5.0  # 5 MeV容差
            println("✅ 两个κ比值在相近温度下同时满足")
        else
            println("⚠️  两个κ比值对应的温度差异较大")
        end
    end
    
catch e
    println("\n✗ 温度反向查找测试失败:")
    println("错误信息: $e")
    println("错误位置: $(stacktrace()[1])")
end

println("\n测试结束时间: $(Dates.now())")
println("="^80)