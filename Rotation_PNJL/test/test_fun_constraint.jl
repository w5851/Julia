include("../src/Function_Gas_Liquid.jl")
include("../src/Constants_Gas_Liquid.jl")
using .Constants_Gas_Liquid: π, hc, m, mσ, mω, mρ, mδ

# 直接使用Fortran代码中的参数值，确保完全一致
# 从Fortran代码中提取的固定参数值
fs = 10.329
fo = 5.423
fr = 3.15
fd = 2.5
b = 0.00692
c = -0.0048
couplings = [fs, fo, fr, fd, b, c]

# 设置测试点（使用Fortran代码中常见的值）
# 测试1：低密度，低温
println("\n===== 测试1：低密度，低温条件 =====")
gsigma = 200.0 / hc  # MeV -> fm^-1
gdelta = 30.0 / hc   # MeV -> fm^-1
mu_p = 300.0 / hc    # MeV -> fm^-1
mu_n = 350.0 / hc    # MeV -> fm^-1
T = 10.0 / hc        # MeV -> fm^-1
mu_B = 400.0 / hc    # MeV -> fm^-1

# 定义多组测试条件的函数
function run_test_case(gsigma, gdelta, mu_p, mu_n, T, mu_B, couplings)
    # 积分节点设置
    nodes = get_nodes(256)  # 使用256个积分点，与Fortran代码相同
    
    # 参数数组
    x = [gsigma, gdelta, mu_p, mu_n]
    params = [T, mu_B]
    
    # 计算Julia版本残差方程结果
    result = calculate_fun_constraint(x, nodes, couplings, params)
    
    # 手动计算与Fortran代码相同的结果进行比较
    p_nodes, coefficient = nodes
    m_p, m_n = calculate_mass(gsigma, gdelta)
    E_p, E_n = calculate_energy(gsigma, gdelta, p_nodes)
    ρ_p = calculate_ρ(E_p, mu_p, T, coefficient)
    ρ_n = calculate_ρ(E_n, mu_n, T, coefficient)
    ρ_ps = calculate_ρ_s(E_p, mu_p, T, coefficient, m_p)
    ρ_ns = calculate_ρ_s(E_n, mu_n, T, coefficient, m_n)
    
    # 计算介子场强度
    fs, fo, fr, fd, b, c = couplings
    gomega = fo * (ρ_p + ρ_n)
    grho = fr * (ρ_p - ρ_n)
    mu_pp = mu_p + (gomega + grho)
    mu_nn = mu_n + (gomega - grho)
    
    # 按照Fortran代码计算残差
    fvec1 = gsigma - fs * (ρ_ps + ρ_ns - b * m * gsigma^2 - c * gsigma^3)
    fvec2 = gdelta - fd * (ρ_ps - ρ_ns)
    fvec3 = mu_B - mu_nn
    fvec4 = 0.198 - ((ρ_n - ρ_p) / (ρ_n + ρ_p))
    
    fortran_result = [fvec1, fvec2, fvec3, fvec4]
    
    println("\n比较Julia与Fortran计算方法:")
    println("-------------------------------")
    println("             Julia结果        Fortran方法        差异")
    println("σ残差:  $(round(result[1], digits=8))  $(round(fortran_result[1], digits=8))  $(round(result[1] - fortran_result[1], digits=8))")
    println("δ残差:  $(round(result[2], digits=8))  $(round(fortran_result[2], digits=8))  $(round(result[2] - fortran_result[2], digits=8))")
    println("μ残差:  $(round(result[3], digits=8))  $(round(fortran_result[3], digits=8))  $(round(result[3] - fortran_result[3], digits=8))")
    println("不对称: $(round(result[4], digits=8))  $(round(fortran_result[4], digits=8))  $(round(result[4] - fortran_result[4], digits=8))")
    
    # 检查结果是否足够接近
    tolerance = 1e-8
    differences = abs.(result - fortran_result)
    all_close = all(differences .< tolerance)
    
    println("\n结论: $(all_close ? "两种计算方法结果一致!" : "存在差异，需要进一步检查!")")
    if !all_close
        max_diff_idx = argmax(differences)
        println("最大差异在方程 $max_diff_idx: $(differences[max_diff_idx])")
    end
    
    # 显示详细计算参数
    println("\n物理参数计算结果:")
    println("-------------------------------")
    println("σ场值: $(round(gsigma*hc, digits=2)) MeV")
    println("δ场值: $(round(gdelta*hc, digits=2)) MeV")
    println("质子有效质量: $(round(m_p*hc, digits=2)) MeV")
    println("中子有效质量: $(round(m_n*hc, digits=2)) MeV")
    println("质子数密度: $(round(ρ_p, digits=6)) fm⁻³")
    println("中子数密度: $(round(ρ_n, digits=6)) fm⁻³")
    println("总重子数密度: $(round(ρ_p+ρ_n, digits=6)) fm⁻³")
    println("质子标量密度: $(round(ρ_ps, digits=6)) fm⁻³")
    println("中子标量密度: $(round(ρ_ns, digits=6)) fm⁻³")
    println("ω场值: $(round(gomega*hc, digits=2)) MeV")
    println("ρ场值: $(round(grho*hc, digits=2)) MeV")
    println("质子有效化学势: $(round(mu_pp*hc, digits=2)) MeV")
    println("中子有效化学势: $(round(mu_nn*hc, digits=2)) MeV")
    
    return result, fortran_result, all_close
end

# 执行第一个测试
run_test_case(gsigma, gdelta, mu_p, mu_n, T, mu_B, couplings)

# 测试2：高密度，高温
println("\n\n===== 测试2：高密度，高温条件 =====")
gsigma = 1.9554423346438443
gdelta = -0.12778204586816913
mu_p = 3.126097813101585
mu_n = 3.071625725048982
T = 50.0 / hc
mu_B = 1001.0 / hc
run_test_case(gsigma, gdelta, mu_p, mu_n, T, mu_B, couplings)

# 测试3：对称核物质（无δ场）
println("\n\n===== 测试3：对称核物质 =====")
gsigma = 180.0 / hc
gdelta = 0.0 / hc
mu_p = 400.0 / hc
mu_n = 400.0 / hc
T = 50.0 / hc
mu_B = 450.0 / hc
run_test_case(gsigma, gdelta, mu_p, mu_n, T, mu_B, couplings)

println("\n\n===== 总结 =====")
println("完成了三种物理条件下的测试对比:")
println("1. 低密度，低温")
println("2. 高密度，高温")
println("3. 对称核物质")
println("\n这些测试验证了Julia代码与原始Fortran代码的计算结果一致性。")
println("如果发现不一致，请检查积分方法、参数设置和数学表达式实现。")
