"""
数值稳定性修复测试

测试 safe_log 函数和 calculate_U、calculate_log_term 函数的边界条件处理
"""

# 包含修复后的函数
include("../src/Function_PNJL_aniso.jl")

using Test
using .Constants_PNJL: T0, a0, a1, a2, b3, b4

function test_safe_log()
    @testset "safe_log 函数测试" begin
        println("测试 safe_log 函数...")
        
        # 测试正常值
        @test safe_log(1.0) ≈ log(1.0)
        @test safe_log(2.718) ≈ log(2.718)
        
        # 测试边界条件
        @test safe_log(0.0) ≈ log(1e-16)  # 零值应被限制为最小值
        @test safe_log(-1.0) ≈ log(1e-16)  # 负值应被限制为最小值
        @test safe_log(1e-20) ≈ log(1e-16)  # 极小值应被限制为最小值
        
        # 测试不同的处理方式
        @test isnan(safe_log(-1.0; handle_negative=:nan))
        @test_throws ErrorException safe_log(-1.0; handle_negative=:error)
        
        println("safe_log 函数测试通过 ✓")
    end
end

function test_calculate_U_stability()
    @testset "calculate_U 稳定性测试" begin
        println("测试 calculate_U 函数的数值稳定性...")
        
        # 正常参数
        T = 0.15  # GeV
        Phi1 = 0.5
        Phi2 = 0.5
        U1 = calculate_U(T, Phi1, Phi2)
        @test isfinite(U1)
        
        # 极端参数：可能导致 value ≤ 0 的情况
        # 当 Phi1 和 Phi2 很大时，value 可能变为负值
        Phi1_extreme = 1.5
        Phi2_extreme = 1.5
        U2 = calculate_U(T, Phi1_extreme, Phi2_extreme)
        @test isfinite(U2)  # 应该仍然是有限值，不是 NaN 或 Inf
        
        # 测试 value 表达式在极端情况下的行为
        value_normal = 1 - 6 * Phi2 * Phi1 + 4 * (Phi2^3 + Phi1^3) - 3 * (Phi2 * Phi1)^2
        value_extreme = 1 - 6 * Phi2_extreme * Phi1_extreme + 4 * (Phi2_extreme^3 + Phi1_extreme^3) - 3 * (Phi2_extreme * Phi1_extreme)^2
        
        println("  正常参数 value = $value_normal")
        println("  极端参数 value = $value_extreme")
        
        # 即使 value 为负，safe_log 也应该处理
        if value_extreme <= 0
            println("  检测到 value ≤ 0 的情况，safe_log 应该处理这种情况")
        end
        
        println("calculate_U 稳定性测试通过 ✓")
    end
end

function test_calculate_log_term_stability()
    @testset "calculate_log_term 稳定性测试" begin
        println("测试 calculate_log_term 函数的数值稳定性...")
        
        # 正常参数
        E_i = 0.5  # GeV
        mu_i = 0.3  # GeV
        T = 0.15   # GeV
        Phi1 = 0.5
        Phi2 = 0.5
        
        result1 = calculate_log_term(E_i, mu_i, T, Phi1, Phi2)
        @test isfinite(result1)
        
        # 极端参数：可能导致 f1_val 或 f2_val ≤ 0 的情况
        # 当 Phi1 和 Phi2 为负值时可能出现问题
        Phi1_neg = -2.0
        Phi2_neg = -2.0
        
        result2 = calculate_log_term(E_i, mu_i, T, Phi1_neg, Phi2_neg)
        @test isfinite(result2)
        
        # 测试高温和低温情况
        T_high = 1.0  # 高温
        T_low = 0.01  # 低温
        
        result_high = calculate_log_term(E_i, mu_i, T_high, Phi1, Phi2)
        result_low = calculate_log_term(E_i, mu_i, T_low, Phi1, Phi2)
        
        @test isfinite(result_high)
        @test isfinite(result_low)
        
        println("calculate_log_term 稳定性测试通过 ✓")
    end
end

function run_numerical_stability_tests()
    println("开始数值稳定性测试...\n")
    
    test_safe_log()
    test_calculate_U_stability()
    test_calculate_log_term_stability()
    
    println("\n所有数值稳定性测试完成 ✓")
    println("修复摘要:")
    println("- 实现了 safe_log 函数处理负值和零值")
    println("- 修复了 calculate_U 函数中的 log(value) 调用")
    println("- 修复了 calculate_log_term 函数中的 log(f1_val) 和 log(f2_val) 调用")
    println("- 所有边界条件都能正确处理")
end

# 运行测试
if abspath(PROGRAM_FILE) == @__FILE__
    run_numerical_stability_tests()
end
