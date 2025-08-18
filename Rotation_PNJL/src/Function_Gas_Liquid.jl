include("Constants_Gas_Liquid.jl")
include("init3.jl")
using .Constants_Gas_Liquid: π, hc, m, mσ, mω, mρ, mδ, calculate_couplings
using .init3: gauleg
using NLsolve
using FiniteDifferences

using BenchmarkTools

function get_nodes(n_p)
    p_nodes,p_weights = gauleg(0,20.0,n_p)
    coefficient = @. p_nodes^2 * p_weights / π^2
    return [p_nodes,coefficient]
end

@inline function fermion(E,μ,T)
    """费米子分布函数"""
    return 1 / (exp((E - μ) / T) + 1)
end

@inline function fermion_anti(E,μ,T)
    """反费米子分布函数"""
    return 1 / (exp((E + μ) / T) + 1)
end

@inline function calculate_log(E,μ,T)
    x = E - μ
    x_anti = E + μ
    term1 = 1 + exp(-x / T)
    term2 = 1 + exp(-x_anti / T)
    return log(term1) + log(term2)
end

@inline function calculate_mass(gσ,gδ)
    m_p = m-gσ-gδ
    m_n = m-gσ+gδ
    return m_p,m_n
end

@inline function calculate_energy(gσ,gδ,p_nodes)
    m_p,m_n = calculate_mass(gσ,gδ)
    p2 = @. p_nodes^2
    E_p = @. sqrt(p2 + m_p^2)
    E_n = @. sqrt(p2 + m_n^2)
    return E_p,E_n
end

@inline function calculate_ρ(E,μ,T,coef)
    # 避免创建中间数组，使用 map-reduce 模式
    return mapreduce(i -> (fermion(E[i],μ,T) - fermion_anti(E[i],μ,T))*coef[i], +, eachindex(E))
end

@inline function calculate_ρ_s(E,μ,T,coef,m)
    # 避免创建中间数组，使用 map-reduce 模式
    return mapreduce(i -> (fermion(E[i],μ,T) + fermion_anti(E[i],μ,T))*coef[i]*m/E[i], +, eachindex(E))
end

@inline function calculate_σ_term(gσ,ρ_ps,ρ_ns,couplings)
    fσ, _, _, _, b, c = couplings
    return - gσ + fσ * (ρ_ps + ρ_ns - b * m * gσ^2 - c * gσ^3)
end

@inline function calculate_δ_term(gδ,ρ_ps,ρ_ns,couplings)
    fδ = couplings[4]
    return -gδ + fδ * (ρ_ps - ρ_ns)
end

@inline function calculate_ρ_term(gρ,ρ_p,ρ_n,couplings)
    fρ = couplings[3]
    return -gρ + fρ * (ρ_p - ρ_n)
end

@inline function calculate_ω_term(gω,ρ_p,ρ_n,couplings)
    fω = couplings[2]
    return -gω + fω * (ρ_p + ρ_n)
end

@inline function calculate_chemical_constraint(μ_B, μ_n,ρ_p, ρ_n,couplings)
    """计算重子化学势约束条件(同位旋不对称:质子中子数密度不等)"""
    gω = calculate_ω_term(0.0, ρ_p, ρ_n, couplings)
    gρ = calculate_ρ_term(0.0, ρ_p, ρ_n, couplings)
    return μ_B - μ_n - gω + gρ
end

@inline function calculate_asymmetry_constraint(ρ_n, ρ_p, target_asymmetry=0.198)
    """计算同位旋不对称度约束条件"""
    return target_asymmetry - (ρ_n - ρ_p)/(ρ_n + ρ_p)
end

function calculate_fun_constraint(x,nodes,couplings,params)
    """计算化学势约束条件下的残差方程"""
    gσ,gδ,μ_p,μ_n = x
    p_nodes, coefficient = nodes
    T = params[1]
    μ_B = params[2]
    
    m_p, m_n = calculate_mass(gσ, gδ)
    E_p, E_n = calculate_energy(gσ, gδ, p_nodes)
    ρ_p = calculate_ρ(E_p, μ_p, T, coefficient)
    ρ_n = calculate_ρ(E_n, μ_n, T, coefficient)
    ρ_ps = calculate_ρ_s(E_p, μ_p, T, coefficient, m_p)
    ρ_ns = calculate_ρ_s(E_n, μ_n, T, coefficient, m_n)

    σ_term = calculate_σ_term(gσ, ρ_ps, ρ_ns, couplings)
    δ_term = calculate_δ_term(gδ, ρ_ps, ρ_ns, couplings)
    # 计算重子化学势约束
    chem_constraint = calculate_chemical_constraint(μ_B, μ_n, ρ_p, ρ_n, couplings)
    
    # 计算同位旋不对称度约束
    asymmetry_constraint = calculate_asymmetry_constraint(ρ_n, ρ_p)

    return [σ_term, δ_term, chem_constraint, asymmetry_constraint]
end

function solve_fun_constraints(x0, nodes,couplings, params)
    """求解化学势约束条件"""
       
    # 使用NLsolve求解
    result = nlsolve(x -> calculate_fun_constraint(x, nodes, couplings, params), x0)
    
    return result.zero
end

@inline function calculate_init_term(E,μ,T,nodes)
    """计算积分项"""
    p, coef = nodes
    return mapreduce(i -> (fermion(E[i],μ,T) + fermion_anti(E[i],μ,T))*coef[i]*p[i]^2/E[i], +, eachindex(E))/3.0
end

@inline function calculate_pressure(gσ, gδ, gω, gρ, μ_p, μ_n, T, nodes, couplings)
    """计算压强"""
    fσ, fω, fρ, fδ, b, c = couplings
    p_nodes, _ = nodes
    
    # 计算质子和中子的有效质量和能量
    E_p, E_n = calculate_energy(gσ, gδ, p_nodes)

    # 计算p_p和p_n（动量积分项）
    p_p = calculate_init_term(E_p, μ_p, T, nodes)
    p_n = calculate_init_term(E_n, μ_n, T, nodes)
    
    # 计算压强
    pressure = -(1.0/3.0) * b * m * gσ^3 - 
               (1.0/4.0) * c * gσ^4 - 
               (1.0/(2.0*fσ)) * gσ^2 + 
               (1.0/(2.0*fω)) * gω^2 + 
               p_p + p_n + 
               (1.0/(2.0*fρ)) * gρ^2 - 
               (1.0/(2.0*fδ)) * gδ^2
               
    return pressure
end

@inline function calculate_pressure_wrapper(x, nodes, couplings, params)
    """计算压强（自动计算场量gω和gρ）"""
    p_nodes, coef = nodes
    T, _ = params
    gσ, gδ, μ_p, μ_n = x

    # 计算密度用于确定场量
    E_p, E_n = calculate_energy(gσ, gδ, p_nodes)
    ρ_p = calculate_ρ(E_p, μ_p, T, coef)
    ρ_n = calculate_ρ(E_n, μ_n, T, coef)
    
    # 根据自洽条件计算场量
    gω = calculate_ω_term(0.0, ρ_p, ρ_n, couplings)
    gρ = calculate_ρ_term(0.0, ρ_p, ρ_n, couplings)

    return calculate_pressure(gσ, gδ, gω, gρ, μ_p, μ_n, T, nodes, couplings)
end

function calculate_pressure_solved(μ_B,T, x0, nodes, couplings)
    """计算压强，使用求解后的场量"""
    params = [T, μ_B]
    x = solve_fun_constraints(x0, nodes, couplings, params)
    return calculate_pressure_wrapper(x, nodes, couplings, params)
end

function calculate_pressure_derivatives(μ_B, T, x0, nodes, couplings; order=4, method=central_fdm(5, 1))
    """
    计算压强对重子化学势μ_B的一到四阶导数
    
    参数:
    - μ_B: 重子化学势
    - T: 温度
    - x0: 初始猜测值 [gσ, gδ, μ_p, μ_n]
    - nodes: 积分节点和权重
    - couplings: 耦合常数 [fσ, fω, fρ, fδ, b, c]
    - order: 计算导数的最高阶数 (默认4)
    - method: 有限差分方法 (默认5点中心差分)
    
    返回:
    - pressure: 压强值
    - derivatives: 包含一到四阶导数的数组
    """
    
    # 定义压强函数，只依赖于μ_B
    pressure_func = μ -> calculate_pressure_solved(μ, T, x0, nodes, couplings)
    
    # 计算压强值
    pressure = pressure_func(μ_B)
    
    # 计算各阶导数
    derivatives = zeros(order)
    
    for i in 1:order
        # 创建相应阶数的有限差分方法
        fdm_method = central_fdm(5, i)
        derivatives[i] = fdm_method(pressure_func, μ_B)
    end
    
    return pressure, derivatives
end

function calculate_pressure_derivatives_efficient(μ_B, T, x0, nodes, couplings; h=1e-5)
    """
    高效计算压强对重子化学势μ_B的一到四阶导数
    使用预定义的有限差分方法以提高性能
    
    参数:
    - μ_B: 重子化学势
    - T: 温度  
    - x0: 初始猜测值 [gσ, gδ, μ_p, μ_n]
    - nodes: 积分节点和权重
    - couplings: 耦合常数 [fσ, fω, fρ, fδ, b, c]
    - h: 步长 (默认1e-5)
    
    返回:
    - pressure: 压强值
    - dpre_dmu1: 一阶导数 ∂P/∂μ_B
    - dpre_dmu2: 二阶导数 ∂²P/∂μ_B²
    - dpre_dmu3: 三阶导数 ∂³P/∂μ_B³
    - dpre_dmu4: 四阶导数 ∂⁴P/∂μ_B⁴
    """
    
    # 定义压强函数
    pressure_func = μ -> calculate_pressure_solved(μ, T, x0, nodes, couplings)
    
    # 预定义有限差分方法
    fdm1 = central_fdm(5, 1)  # 一阶导数，5点中心差分
    fdm2 = central_fdm(5, 2)  # 二阶导数，5点中心差分
    fdm3 = central_fdm(7, 3)  # 三阶导数，7点中心差分
    fdm4 = central_fdm(7, 4)  # 四阶导数，7点中心差分
    
    # 计算压强和各阶导数
    pressure = pressure_func(μ_B)
    dpre_dmu1 = fdm1(pressure_func, μ_B)
    dpre_dmu2 = fdm2(pressure_func, μ_B)
    dpre_dmu3 = fdm3(pressure_func, μ_B)
    dpre_dmu4 = fdm4(pressure_func, μ_B)
    
    return pressure, dpre_dmu1, dpre_dmu2, dpre_dmu3, dpre_dmu4
end

function calculate_thermodynamic_fluctuations(μ_B, T, x0, nodes, couplings)
    """
    计算热力学涨落相关量
    基于压强对化学势的导数计算累积量和涨落
    
    参数:
    - μ_B: 重子化学势
    - T: 温度
    - x0: 初始猜测值 [gσ, gδ, μ_p, μ_n] 
    - nodes: 积分节点和权重
    - couplings: 耦合常数 [fσ, fω, fρ, fδ, b, c]
    
    返回:
    - kappa1: 第一累积量 (数密度)
    - kappa2: 第二累积量 (方差)
    - kappa3: 第三累积量 (偏度)
    - kappa4: 第四累积量 (峰度)
    - fluctuation_ratios: 涨落比值 [κ₂/κ₁, κ₃/κ₂, κ₄/κ₂]
    """
    
    # 计算压强和导数
    pressure, dpre_dmu1, dpre_dmu2, dpre_dmu3, dpre_dmu4 = 
        calculate_pressure_derivatives_efficient(μ_B, T, x0, nodes, couplings)
    
    # 根据Fortran代码中的公式计算累积量
    # 注意：这里假设压强已经按T⁴归一化
    kappa1 = dpre_dmu1  # 第一累积量：数密度
    kappa2 = dpre_dmu2  # 第二累积量：方差
    kappa3 = dpre_dmu3  # 第三累积量：偏度
    kappa4 = dpre_dmu4  # 第四累积量：峰度
    
    # 计算涨落比值
    fluctuation_ratios = [
        kappa2 / kappa1,   # κ₂/κ₁ (归一化方差)
        kappa3 / kappa2,   # κ₃/κ₂ (偏度相关)
        kappa4 / kappa2    # κ₄/κ₂ (峰度相关)
    ]
    
    return kappa1, kappa2, kappa3, kappa4, fluctuation_ratios
end

function calculate_derivatives_batch(μ_B_array, T, x0, nodes, couplings; save_results=false, output_file=joinpath(@__DIR__, "..", "output", "derivatives_output.dat"))
    """
    批量计算多个化学势点的压强导数
    模拟Fortran代码中的循环计算过程
    
    参数:
    - μ_B_array: 化学势数组
    - T: 温度
    - x0: 初始猜测值 [gσ, gδ, μ_p, μ_n]
    - nodes: 积分节点和权重
    - couplings: 耦合常数
    - save_results: 是否保存结果到文件
    - output_file: 输出文件名
    
    返回:
    - results: 包含所有计算结果的NamedTuple
    """
    
    n_points = length(μ_B_array)
    
    # 预分配结果数组
    pressure_array = zeros(n_points)
    dpre_dmu1_array = zeros(n_points)
    dpre_dmu2_array = zeros(n_points)
    dpre_dmu3_array = zeros(n_points)
    dpre_dmu4_array = zeros(n_points)
    kappa1_array = zeros(n_points)
    kappa2_array = zeros(n_points)
    kappa3_array = zeros(n_points)
    kappa4_array = zeros(n_points)
    fluctuation_ratios_array = zeros(n_points, 3)
    
    println("开始批量计算 $(n_points) 个化学势点的压强导数...")
    
    # 循环计算每个化学势点
    for (i, μ_B) in enumerate(μ_B_array)
        try
            # 计算压强和导数
            pressure, dpre_dmu1, dpre_dmu2, dpre_dmu3, dpre_dmu4 = 
                calculate_pressure_derivatives_efficient(μ_B, T, x0, nodes, couplings)
            
            # 计算热力学涨落
            kappa1, kappa2, kappa3, kappa4, fluctuation_ratios = 
                calculate_thermodynamic_fluctuations(μ_B, T, x0, nodes, couplings)
            
            # 存储结果
            pressure_array[i] = pressure
            dpre_dmu1_array[i] = dpre_dmu1
            dpre_dmu2_array[i] = dpre_dmu2
            dpre_dmu3_array[i] = dpre_dmu3
            dpre_dmu4_array[i] = dpre_dmu4
            kappa1_array[i] = kappa1
            kappa2_array[i] = kappa2
            kappa3_array[i] = kappa3
            kappa4_array[i] = kappa4
            fluctuation_ratios_array[i, :] = fluctuation_ratios
            
            # 进度报告
            if i % 10 == 0 || i == n_points
                println("已完成: $(i)/$(n_points) ($(round(i/n_points*100, digits=1))%)")
            end
            
        catch e
            println("警告: 化学势 μ_B = $(μ_B*hc) MeV 处计算失败: $e")
            # 填入NaN值
            pressure_array[i] = NaN
            dpre_dmu1_array[i] = NaN
            dpre_dmu2_array[i] = NaN
            dpre_dmu3_array[i] = NaN
            dpre_dmu4_array[i] = NaN
            kappa1_array[i] = NaN
            kappa2_array[i] = NaN
            kappa3_array[i] = NaN
            kappa4_array[i] = NaN
            fluctuation_ratios_array[i, :] .= NaN
        end
    end
    
    # 创建结果NamedTuple
    results = (
        μ_B = μ_B_array,
        T = T,
        pressure = pressure_array,
        dpre_dmu1 = dpre_dmu1_array,
        dpre_dmu2 = dpre_dmu2_array,
        dpre_dmu3 = dpre_dmu3_array,
        dpre_dmu4 = dpre_dmu4_array,
        kappa1 = kappa1_array,
        kappa2 = kappa2_array,
        kappa3 = kappa3_array,
        kappa4 = kappa4_array,
        fluctuation_ratios = fluctuation_ratios_array
    )
    
    # 保存结果到文件
    if save_results
        save_derivatives_results(results, output_file)
        println("结果已保存到文件: $output_file")
    end
    
    println("批量计算完成!")
    return results
end

function save_derivatives_results(results, filename)
    """保存导数计算结果到文件"""
    println("正在保存结果到文件: $filename")
    
    # 确保输出目录存在
    output_dir = dirname(filename)
    if !isdir(output_dir)
        println("创建输出目录: $output_dir")
        mkpath(output_dir)
    end
    
    try
        open(filename, "w") do io
            # 写入文件头
            write(io, "# 压强导数计算结果\n")
            write(io, "# 列: μ_B(MeV) T(MeV) P ∂P/∂μ ∂²P/∂μ² ∂³P/∂μ³ ∂⁴P/∂μ⁴ κ₁ κ₂ κ₃ κ₄ κ₂/κ₁ κ₃/κ₂ κ₄/κ₂\n")
            
            # 写入数据
            for i in eachindex(results.μ_B)
                write(io, "$(results.μ_B[i]*hc) $(results.T*hc) $(results.pressure[i]) ")
                write(io, "$(results.dpre_dmu1[i]) $(results.dpre_dmu2[i]) $(results.dpre_dmu3[i]) $(results.dpre_dmu4[i]) ")
                write(io, "$(results.kappa1[i]) $(results.kappa2[i]) $(results.kappa3[i]) $(results.kappa4[i]) ")
                write(io, "$(results.fluctuation_ratios[i,1]) $(results.fluctuation_ratios[i,2]) $(results.fluctuation_ratios[i,3])\n")
            end
        end
        println("文件保存成功!")
    catch e
        println("保存文件时出错: $e")
        rethrow(e)
    end
end


#示例用法:

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

# 单点计算示例
μ_B = 1001.0/hc
params = [T, μ_B]
mu_p = μ_B / 2.0
mu_n = μ_B / 2.0
x0 = [gsigma, gdelta, mu_p, mu_n]

# 方法1: 使用通用导数计算函数
pressure, derivatives = calculate_pressure_derivatives(μ_B, T, x0, nodes, couplings)
println("压强: ", pressure)
println("一阶导数: ", derivatives[1])
println("二阶导数: ", derivatives[2])
println("三阶导数: ", derivatives[3])
println("四阶导数: ", derivatives[4])

# 方法2: 使用高效导数计算函数
pressure, dpre_dmu1, dpre_dmu2, dpre_dmu3, dpre_dmu4 = 
    calculate_pressure_derivatives_efficient(μ_B, T, x0, nodes, couplings)

# 方法3: 计算热力学涨落
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

# 批量计算示例
μ_B_range = 1001.0/hc:-10.0/hc:600.0/hc  # 从1001 MeV开始
results = calculate_derivatives_batch(μ_B_range, T, x0, nodes, couplings, 
                                    save_results=true, output_file=joinpath(@__DIR__, "..", "output", "pressure_derivatives.dat"))
