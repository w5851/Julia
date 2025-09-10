#!/usr/bin/env julia
"""
更精确模拟Fortran逐阶求导算法的实现

基于提供的Fortran代码：
- derivation subroutine: 逐阶求导主程序
- HCHIR subroutine: 自适应多项式拟合算法

重要参数：NN=5, M=4, M1=5
"""

include("../../src/Function_Gas_Liquid.jl")

function hchir_julia(x, y, n, m, m1)
    """
    Julia实现的HCHIR子程序
    自适应多项式拟合算法，选择最佳插值点
    
    参数:
    - x: 自变量数组
    - y: 因变量数组  
    - n: 数据点数
    - m: 多项式最高次数
    - m1: m+1
    
    返回:
    - a: 多项式系数数组
    """
    
    # 初始化
    a = zeros(m1)
    if m >= n
        m = n - 1
    end
    if m >= 20
        m = 19
    end
    m1 = m + 1
    
    ha = 0.0
    ix = zeros(Int, 20)
    h = zeros(20)
    
    # 初始点选择
    ix[1] = 1
    ix[m1] = n
    l = div(n-1, m)
    j = l
    for i in 2:m
        ix[i] = j + 1
        j = j + l
    end
    
    max_iterations = 50  # 防止无限循环
    iteration = 0
    
    while iteration < max_iterations
        iteration += 1
        
        # 开始拟合过程
        hh = 1.0
        for i in 1:m1
            a[i] = y[ix[i]]
            h[i] = -hh
            hh = -hh
        end
        
        # 计算差商表
        for j in 1:m
            ii = m1
            y2 = a[ii]
            h2 = h[ii]
            for i in j:m
                d = x[ix[ii]] - x[ix[m1-i]]
                y1 = a[m-i+j]
                h1 = h[m-i+j]
                a[ii] = (y2 - y1) / d
                h[ii] = (h2 - h1) / d
                ii = m - i + j
                y2 = y1
                h2 = h1
            end
        end
        
        # 调整系数
        hh = -a[m1] / h[m1]
        for i in 1:m1
            a[i] = a[i] + h[i] * hh
        end
        
        # 反向计算最终系数
        for j in 1:(m-1)
            ii = m - j
            d = x[ix[ii]]
            y2 = a[ii]
            for k in (m1-j):m
                y1 = a[k]
                a[ii] = y2 - d * y1
                y2 = y1
                ii = k
            end
        end
        
        # 检查拟合质量
        hm = abs(hh)
        if hm <= ha
            a[m1] = -hm
            break
        end
        
        a[m1] = hm
        ha = hm
        im = ix[1]
        h1 = hh
        j = 1
        
        # 寻找最大误差点
        for i in 1:n
            if i == ix[j] && j < m1
                j = j + 1
            else
                # 计算多项式值
                h2 = a[m]
                for k in (m-1):-1:1
                    h2 = h2 * x[i] + a[k]
                end
                h2 = h2 - y[i]
                
                if abs(h2) > hm
                    hm = abs(h2)
                    h1 = h2
                    im = i
                end
            end
        end
        
        if im == ix[1]
            break
        end
        
        # 更新插值点
        i = 1
        while i <= m1 && im >= ix[i]
            i = i + 1
        end
        
        if i > m1
            i = m1
        end
        
        h2 = (i % 2 == 0) ? hh : -hh
        
        if h1 * h2 >= 0.0
            ix[i] = im
            continue
        end
        
        if im < ix[1]
            for j in m:-1:1
                ix[j+1] = ix[j]
            end
            ix[1] = im
        elseif im > ix[m1]
            for j in 2:m1
                ix[j-1] = ix[j]
            end
            ix[m1] = im
        else
            ix[i-1] = im
        end
    end
    
    return a[1:m1]  # 返回多项式系数
end

function evaluate_polynomial_hchir(coeffs, x)
    """
    使用HCHIR得到的系数计算多项式值
    """
    result = coeffs[end]
    for i in (length(coeffs)-1):-1:1
        result = result * x + coeffs[i]
    end
    return result
end

function evaluate_polynomial_derivative_hchir(coeffs, x, order)
    """
    使用HCHIR得到的系数计算多项式导数
    """
    if order == 0
        return evaluate_polynomial_hchir(coeffs, x)
    end
    
    # 计算导数系数
    deriv_coeffs = copy(coeffs)
    for ord in 1:order
        new_coeffs = zeros(length(deriv_coeffs)-1)
        for i in 1:(length(deriv_coeffs)-1)
            new_coeffs[i] = deriv_coeffs[i+1] * i
        end
        deriv_coeffs = new_coeffs
        if length(deriv_coeffs) == 0
            return 0.0
        end
    end
    
    return evaluate_polynomial_hchir(deriv_coeffs, x)
end

function fortran_exact_derivative(x_points, y_values, target_order)
    """
    精确实现Fortran的derivation子程序
    
    参数:
    - x_points: 自变量数组（对应Fortran中的b）
    - y_values: 因变量数组（对应Fortran中的a）
    - target_order: 目标导数阶数
    
    返回:
    - 各阶导数数组
    """
    
    n = length(x_points)
    current_values = copy(y_values)
    
    # Fortran参数设置
    NN = 5
    M = 4
    M1 = 5
    
    println("  使用Fortran精确参数: NN=$NN, M=$M, M1=$M1")
    
    for order in 1:target_order
        println("    计算第$(order)阶导数...")
        
        new_derivatives = zeros(n)
        
        # 内部点：使用中心差分（对应Fortran的 Do i=2,n-1）
        for i in 2:(n-1)
            new_derivatives[i] = (current_values[i+1] - current_values[i-1]) / (x_points[i+1] - x_points[i-1])
        end
        
        # 左边界点：使用HCHIR拟合（对应 CALL HCHIR(B(2:2+m),c(2:2+m),NN,AA,M,M1)）
        if n >= 6  # 确保有足够的点
            start_idx = 2
            end_idx = min(2 + M, n)
            x_fit = x_points[start_idx:end_idx]
            y_fit = current_values[start_idx:end_idx]
            
            try
                coeffs = hchir_julia(x_fit, y_fit, length(x_fit), M, M1)
                new_derivatives[1] = evaluate_polynomial_derivative_hchir(coeffs, x_points[1], 1)
            catch e
                println("      警告: 左边界HCHIR拟合失败，使用前向差分")
                if n >= 2
                    new_derivatives[1] = (current_values[2] - current_values[1]) / (x_points[2] - x_points[1])
                end
            end
        end
        
        # 右边界点：使用HCHIR拟合（对应 CALL HCHIR(B(n-m1:n-1),c(n-m1:n-1),NN,AA,M,M1)）
        if n >= 6
            start_idx = max(n - M1, 1)
            end_idx = n - 1
            x_fit = x_points[start_idx:end_idx]
            y_fit = current_values[start_idx:end_idx]
            
            try
                coeffs = hchir_julia(x_fit, y_fit, length(x_fit), M, M1)
                new_derivatives[n] = evaluate_polynomial_derivative_hchir(coeffs, x_points[n], 1)
            catch e
                println("      警告: 右边界HCHIR拟合失败，使用后向差分")
                if n >= 2
                    new_derivatives[n] = (current_values[n] - current_values[n-1]) / (x_points[n] - x_points[n-1])
                end
            end
        end
        
        current_values = new_derivatives
        
        # 显示当前阶导数的统计信息
        finite_vals = filter(isfinite, current_values)
        if length(finite_vals) > 0
            println("      第$(order)阶导数范围: $(round(minimum(finite_vals), digits=8)) 到 $(round(maximum(finite_vals), digits=8))")
        end
    end
    
    return current_values
end

function fortran_exact_pressure_derivatives(μ_B_center, T, x0, nodes, couplings; 
                                          h_base=1e-3, n_points=11)
    """
    使用Fortran精确算法计算压强导数
    """
    
    println("使用Fortran精确逐阶求导方法")
    println("基础步长: h = $(h_base), 网格点数: $(n_points)")
    
    # 构建网格点
    half_range = (n_points - 1) ÷ 2
    μ_points = [μ_B_center + (i - half_range - 1) * h_base for i in 1:n_points]
    
    # 计算无量纲压强函数
    normalized_pressure_func = μ -> begin
        pressure = calculate_pressure_solved(μ, T, x0, nodes, couplings)
        return pressure / T^4
    end
    
    println("计算$(n_points)个网格点的无量纲压强...")
    pressure_values = zeros(n_points)
    for i in 1:n_points
        try
            pressure_values[i] = normalized_pressure_func(μ_points[i])
            if i == 1 || i == n_points || i == (n_points+1)÷2
                println("  μ/T = $(round(μ_points[i]/T, digits=4)), P/T⁴ = $(round(pressure_values[i], digits=8))")
            end
        catch e
            println("  网格点 $(i) 计算失败: $e")
            pressure_values[i] = NaN
        end
    end
    
    # 检查数据有效性
    if !all(isfinite.(pressure_values))
        error("某些网格点计算失败，无法进行导数计算")
    end
    
    # 转换为无量纲变量 μ/T
    normalized_mu = μ_points ./ T
    
    println("\n开始Fortran精确逐阶求导计算...")
    
    # 1阶导数
    dpre_dmu1_array = fortran_exact_derivative(normalized_mu, pressure_values, 1)
    
    # 2阶导数  
    dpre_dmu2_array = fortran_exact_derivative(normalized_mu, pressure_values, 2)
    
    # 3阶导数
    dpre_dmu3_array = fortran_exact_derivative(normalized_mu, pressure_values, 3)
    
    # 4阶导数
    dpre_dmu4_array = fortran_exact_derivative(normalized_mu, pressure_values, 4)
    
    # 取中心点的导数值
    center_idx = (n_points + 1) ÷ 2
    center_pressure = pressure_values[center_idx]
    center_dpre_dmu1 = dpre_dmu1_array[center_idx]
    center_dpre_dmu2 = dpre_dmu2_array[center_idx]
    center_dpre_dmu3 = dpre_dmu3_array[center_idx]
    center_dpre_dmu4 = dpre_dmu4_array[center_idx]
    
    println("\n中心点 μ/T = $(round(μ_B_center/T, digits=4)) 处的Fortran精确结果:")
    println("P/T⁴ = $(round(center_pressure, digits=8))")
    println("∂(P/T⁴)/∂(μ/T) = $(round(center_dpre_dmu1, digits=8))")
    println("∂²(P/T⁴)/∂(μ/T)² = $(round(center_dpre_dmu2, digits=8))")
    println("∂³(P/T⁴)/∂(μ/T)³ = $(round(center_dpre_dmu3, digits=8))")
    println("∂⁴(P/T⁴)/∂(μ/T)⁴ = $(round(center_dpre_dmu4, digits=8))")
    
    return center_pressure, center_dpre_dmu1, center_dpre_dmu2, center_dpre_dmu3, center_dpre_dmu4
end

function fortran_exact_thermodynamic_fluctuations(μ_B, T, x0, nodes, couplings; 
                                                 h_base=1e-3, n_points=11)
    """
    使用Fortran精确方法计算热力学涨落
    """
    
    # 计算无量纲压强及其导数
    pressure_norm, dpre_dmu1, dpre_dmu2, dpre_dmu3, dpre_dmu4 = 
        fortran_exact_pressure_derivatives(μ_B, T, x0, nodes, couplings, 
                                         h_base=h_base, n_points=n_points)
    
    # 累积量就是导数本身（无量纲化处理已包含）
    kappa1 = dpre_dmu1
    kappa2 = dpre_dmu2
    kappa3 = dpre_dmu3
    kappa4 = dpre_dmu4
    
    println("\n热力学累积量（Fortran精确方法）:")
    println("κ₁ = $(round(kappa1, digits=8))")
    println("κ₂ = $(round(kappa2, digits=8))")
    println("κ₃ = $(round(kappa3, digits=8))")
    println("κ₄ = $(round(kappa4, digits=8))")
    
    return kappa1, kappa2, kappa3, kappa4
end

# 主程序测试
println("="^60)
println("Fortran精确逐阶求导方法测试")
println("="^60)

# 设置参数（与demo_temperature_scan.jl一致）
nodes = get_nodes(256)
gsigma = 1.25
gdelta = 0.01
fs = 10.329
fo = 5.423
fr = 3.15
fd = 0.0
b = 0.00692
c = -0.0048
couplings = [fs, fo, fr, fd, b, c]

μ_B_test = 697.0/hc  # 697 MeV
T_test = 150.0/hc     # 50 MeV

mu_p = μ_B_test / 2.0
mu_n = μ_B_test / 2.0
x0 = [gsigma, gdelta, mu_p, mu_n]

println("测试参数:")
println("μ_B = $(μ_B_test*hc) MeV")
println("T = $(T_test*hc) MeV")
println("μ_B/T = $(round(μ_B_test/T_test, digits=2))")

try
    # 首先求解自洽方程
    params = [T_test, μ_B_test]
    converged_solution = solve_fun_constraints(x0, nodes, couplings, params)
    
    println("\n收敛解:")
    println("gσ = $(round(converged_solution[1], digits=6))")
    println("gδ = $(round(converged_solution[2], digits=6))")
    println("μ_p = $(round(converged_solution[3]*hc, digits=2)) MeV")
    println("μ_n = $(round(converged_solution[4]*hc, digits=2)) MeV")
    
    # 计算热力学涨落
    println("\n" * "="^40)
    kappa1, kappa2, kappa3, kappa4 = fortran_exact_thermodynamic_fluctuations(
        μ_B_test, T_test, converged_solution, nodes, couplings, 
        h_base=1e-3, n_points=11)
    
    # 计算涨落比值
    println("\n涨落比值（Fortran精确方法）:")
    if abs(kappa1) > 1e-15
        ratio_31 = kappa3 / kappa1
        println("κ₃/κ₁ = $(round(ratio_31, digits=6))")
    else
        println("κ₃/κ₁ = NaN (κ₁接近零)")
    end
    
    if abs(kappa2) > 1e-15
        ratio_42 = kappa4 / kappa2
        println("κ₄/κ₂ = $(round(ratio_42, digits=6))")
    else
        println("κ₄/κ₂ = NaN (κ₂接近零)")
    end
    
catch e
    println("计算失败: $e")
    println("错误详情: ", sprint(showerror, e))
end

println("\n" * "="^60)
println("Fortran精确方法测试完成!")
println("="^60)
