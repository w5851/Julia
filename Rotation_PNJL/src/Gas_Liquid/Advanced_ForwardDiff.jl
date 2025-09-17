# Advanced_ForwardDiff.jl
# ForwardDiff自动微分高阶函数模块：用于PNJL模型的高精度导数计算
# 依赖于 Function_Gas_Liquid.jl 中的基础函数

include("Function_Gas_Liquid.jl")
#include("Constants_Gas_Liquid.jl")
using .Constants_Gas_Liquid: hc
using NLsolve
using ForwardDiff
using CSV
using DataFrames
using Dates

function NewQuark_mu_pnjl_fixed(x_array, T, mu_B, model_params)
    """
    修正的PNJL模型类型保护约束方程求解器
    
    参数:
    - x_array: 待求解的场量和化学势 [gσ, gδ, μ_p, μ_n]
    - T: 温度
    - mu_B: 重子化学势
    - model_params: 模型参数 (nodes, couplings)
    
    返回:
    - fvec: 约束方程残差向量
    """
    nodes, couplings = model_params
    
    # 动态创建与输入类型兼容的数组
    T_out = promote_type(eltype(x_array), typeof(T), typeof(mu_B))
    fvec = zeros(T_out, 4)

    # 调用PNJL核心约束方程计算
    params = [T, mu_B]
    result = calculate_fun_constraint(x_array, nodes, couplings, params)
    
    for i in 1:4
        fvec[i] = result[i]
    end
    
    return fvec
end

function SolveOmega_pnjl_fixed(X0, T, mu_B, model_params)
    """
    修正的PNJL模型Omega求解器
    使用ForwardDiff兼容的求解方式
    
    参数:
    - X0: 初始猜测值 [gσ, gδ, μ_p, μ_n]
    - T: 温度
    - mu_B: 重子化学势
    - model_params: 模型参数 (nodes, couplings)
    
    返回:
    - pressure: 压强值
    """
    # 正确更新初值中的化学势部分
    X0_updated = copy(X0)
    
    # X0 = [gsigma, gdelta, mu_p, mu_n]
    # 在PNJL模型中，应该有 mu_p = mu_n = mu_B/2
    mu_u = mu_B / 2.0  # 对应质子化学势
    mu_d = mu_B / 2.0  # 对应中子化学势
    
    X0_updated[3] = mu_u
    X0_updated[4] = mu_d
    
    # 类型转换
    X0_typed = convert.(promote_type(eltype(X0_updated), typeof(T), typeof(mu_B)), X0_updated)
    
    fWrapper(Xs) = NewQuark_mu_pnjl_fixed(Xs, T, mu_B, model_params)

    res = nlsolve(fWrapper, X0_typed, autodiff=:forward)
    NewX = res.zero
    
    # 计算压强
    nodes, couplings = model_params
    params_typed = [T, mu_B]
    pressure = calculate_pressure_wrapper(NewX, nodes, couplings, params_typed)
    return pressure
end

function create_pressure_function(gsigma, gdelta, T, model_params)
    """
    创建压强函数的闭包，避免类型推断问题
    
    参数:
    - gsigma: σ场初值
    - gdelta: δ场初值  
    - T: 温度
    - model_params: 模型参数 (nodes, couplings)
    
    返回:
    - 压强函数 μ -> P(μ)/T⁴
    """
    return μ -> begin
        # 每次调用时在lambda内部创建新数组，确保类型一致性
        X0_new = [gsigma, gdelta, μ/2.0, μ/2.0]
        SolveOmega_pnjl_fixed(X0_new, T, μ, model_params) / T^4
    end
end

function D1_Pressure_mu(gsigma, gdelta, T, mu_B, model_params)
    """
    一阶导数：∂(P/T⁴)/∂μ_B
    
    参数:
    - gsigma: σ场初值
    - gdelta: δ场初值
    - T: 温度
    - mu_B: 重子化学势
    - model_params: 模型参数
    
    返回:
    - 一阶导数值
    """
    pressure_func = create_pressure_function(gsigma, gdelta, T, model_params)
    return ForwardDiff.derivative(pressure_func, mu_B)
end

function D2_Pressure_mu(gsigma, gdelta, T, mu_B, model_params)
    """
    二阶导数：∂²(P/T⁴)/∂μ_B²
    
    参数:
    - gsigma: σ场初值
    - gdelta: δ场初值
    - T: 温度
    - mu_B: 重子化学势
    - model_params: 模型参数
    
    返回:
    - 二阶导数值
    """
    d1_func = μ -> D1_Pressure_mu(gsigma, gdelta, T, μ, model_params)
    return ForwardDiff.derivative(d1_func, mu_B)
end

function D3_Pressure_mu(gsigma, gdelta, T, mu_B, model_params)
    """
    三阶导数：∂³(P/T⁴)/∂μ_B³
    
    参数:
    - gsigma: σ场初值
    - gdelta: δ场初值
    - T: 温度
    - mu_B: 重子化学势
    - model_params: 模型参数
    
    返回:
    - 三阶导数值
    """
    d2_func = μ -> D2_Pressure_mu(gsigma, gdelta, T, μ, model_params)
    return ForwardDiff.derivative(d2_func, mu_B)
end

function D4_Pressure_mu_enhanced(gsigma, gdelta, T, mu_B, model_params)
    """
    增强的四阶导数计算：使用三阶导数的数值微分
    采用5点中心差分公式提高精度
    
    参数:
    - gsigma: σ场初值
    - gdelta: δ场初值
    - T: 温度
    - mu_B: 重子化学势
    - model_params: 模型参数
    
    返回:
    - 四阶导数值
    """
    h = 1e-3  # 基础步长
    
    # 使用5点中心差分公式计算四阶导数
    # f''''(x) ≈ [f'''(x-2h) - 8f'''(x-h) + 8f'''(x+h) - f'''(x+2h)] / (12h)
    
    try
        d3_m2h = D3_Pressure_mu(gsigma, gdelta, T, mu_B - 2*h, model_params)
        d3_m1h = D3_Pressure_mu(gsigma, gdelta, T, mu_B - h, model_params)
        d3_p1h = D3_Pressure_mu(gsigma, gdelta, T, mu_B + h, model_params)
        d3_p2h = D3_Pressure_mu(gsigma, gdelta, T, mu_B + 2*h, model_params)
        
        # 5点中心差分公式
        d4_enhanced = (d3_m2h - 8*d3_m1h + 8*d3_p1h - d3_p2h) / (12*h)
        
        return d4_enhanced
    catch e
        # 如果计算失败，返回NaN
        return NaN
    end
end

function calculate_forwarddiff_derivatives(gsigma, gdelta, T, mu_B, model_params)
    """
    使用ForwardDiff计算所有阶导数
    
    参数:
    - gsigma: σ场初值
    - gdelta: δ场初值
    - T: 温度
    - mu_B: 重子化学势
    - model_params: 模型参数
    
    返回:
    - (d1, d2, d3, d4): 一到四阶导数元组
    """
    d1 = D1_Pressure_mu(gsigma, gdelta, T, mu_B, model_params)
    d2 = D2_Pressure_mu(gsigma, gdelta, T, mu_B, model_params)
    d3 = D3_Pressure_mu(gsigma, gdelta, T, mu_B, model_params)
    d4 = D4_Pressure_mu_enhanced(gsigma, gdelta, T, mu_B, model_params)
    
    return d1, d2, d3, d4
end

function calculate_forwarddiff_thermodynamic_fluctuations(gsigma, gdelta, T, mu_B, model_params)
    """
    使用ForwardDiff计算热力学涨落量
    
    参数:
    - gsigma: σ场初值
    - gdelta: δ场初值
    - T: 温度
    - mu_B: 重子化学势
    - model_params: 模型参数
    
    返回:
    - (κ1, κ2, κ3, κ4, κ3/κ1, κ4/κ2): 累积量和涨落比值
    """
    # 计算基础压强
    X0 = [gsigma, gdelta, mu_B/2.0, mu_B/2.0]
    pressure = SolveOmega_pnjl_fixed(X0, T, mu_B, model_params)
    pressure_normalized = pressure / T^4
    
    # 计算一到四阶导数
    d1_raw, d2_raw, d3_raw, d4_raw = calculate_forwarddiff_derivatives(gsigma, gdelta, T, mu_B, model_params)
    
    # 转换为对μ/T的导数
    d1_norm = d1_raw * T     # ∂(P/T⁴)/∂(μ/T)
    d2_norm = d2_raw * T^2   # ∂²(P/T⁴)/∂(μ/T)²
    d3_norm = d3_raw * T^3   # ∂³(P/T⁴)/∂(μ/T)³
    d4_norm = d4_raw * T^4   # ∂⁴(P/T⁴)/∂(μ/T)⁴
    
    # 计算累积量
    κ1 = d1_norm
    κ2 = d2_norm  
    κ3 = d3_norm
    κ4 = d4_norm
    
    # 计算涨落比值
    κ3_κ1 = if abs(κ1) > 1e-15
        κ3 / κ1
    else
        NaN
    end
    
    κ4_κ2 = if abs(κ2) > 1e-15 && isfinite(κ4)
        κ4 / κ2
    else
        NaN
    end
    
    return κ1, κ2, κ3, κ4, κ3_κ1, κ4_κ2
end

function forwarddiff_temperature_scan(μ_B, T_min, T_max, T_step, output_file; 
                                     gsigma=1.25, gdelta=0.01, 
                                     fs=17.28476, fo=11.66174, fr=0.89363, fd=0.0, 
                                     b=0.00210, c=-0.00297, n_nodes=256)
    """
    使用ForwardDiff方法进行温度扫描
    
    参数:
    - μ_B: 固定重子化学势
    - T_min: 最小温度
    - T_max: 最大温度
    - T_step: 温度步长
    - output_file: 输出文件路径
    - gsigma, gdelta: 场初值 (可选)
    - fs, fo, fr, fd, b, c: 耦合常数 (可选)
    - n_nodes: 积分节点数 (可选)
    
    返回:
    - DataFrame: 包含所有计算结果的数据框
    """
    
    # 设置模型参数
    nodes = get_nodes(n_nodes)
    couplings = [fs, fo, fr, fd, b, c]
    model_params = (nodes, couplings)
    
    # 生成温度数组
    T_array = T_min:T_step:T_max
    n_temps = length(T_array)
    
    println("="^60)
    println("ForwardDiff温度扫描")
    println("="^60)
    println("固定μ_B = $(μ_B*hc) MeV，$(n_temps) 个温度点")
    println("温度范围: $(T_min*hc) - $(T_max*hc) MeV，步长: $(T_step*hc) MeV")
    println("使用固定初解，不进行迭代更新")
    println("\n模型参数:")
    println("  gsigma = $gsigma")
    println("  gdelta = $gdelta")
    println("  fs = $fs")
    println("  fo = $fo")
    println("  fr = $fr")
    println("  fd = $fd")
    println("  b = $b")
    println("  c = $c")
    println("  nodes = $n_nodes")
    
    # 预分配结果数组
    results = []
    
    # 循环计算每个温度点
    for (i, T) in enumerate(T_array)
        println("\n$(i)/$(n_temps): 计算 T = $(round(T*hc, digits=1)) MeV")
        
        try
            # 计算热力学涨落量
            κ1, κ2, κ3, κ4, κ3_κ1, κ4_κ2 = calculate_forwarddiff_thermodynamic_fluctuations(
                gsigma, gdelta, T, μ_B, model_params)
            
            # 计算基础压强
            X0 = [gsigma, gdelta, μ_B/2.0, μ_B/2.0]
            pressure = SolveOmega_pnjl_fixed(X0, T, μ_B, model_params)
            pressure_normalized = pressure / T^4
            
            # 存储结果
            result_row = (
                T_MeV = T * hc,
                P_T4 = pressure_normalized,
                kappa1 = κ1,
                kappa2 = κ2,
                kappa3 = κ3,
                kappa4 = κ4,
                kappa3_over_kappa1 = κ3_κ1,
                kappa4_over_kappa2 = κ4_κ2,
                mu_over_T = μ_B / T
            )
            
            push!(results, result_row)
            
            println("  P/T⁴ = $(round(pressure_normalized, digits=6))")
            println("  κ₃/κ₁ = $(isfinite(κ3_κ1) ? round(κ3_κ1, digits=6) : "NaN")")
            println("  κ₄/κ₂ = $(isfinite(κ4_κ2) ? round(κ4_κ2, digits=6) : "NaN")")
            
        catch e
            println("  计算失败: $e")
            
            # 添加失败的记录
            result_row = (
                T_MeV = T * hc,
                P_T4 = NaN,
                kappa1 = NaN,
                kappa2 = NaN,
                kappa3 = NaN,
                kappa4 = NaN,
                kappa3_over_kappa1 = NaN,
                kappa4_over_kappa2 = NaN,
                mu_over_T = μ_B / T
            )
            
            push!(results, result_row)
        end
    end
    
    # 保存结果到CSV文件
    save_forwarddiff_results(results, μ_B, T_min, T_max, T_step, output_file, 
                            gsigma, gdelta, fs, fo, fr, fd, b, c, n_nodes)
    
    return DataFrame(results)
end

function save_forwarddiff_results(results, μ_B, T_min, T_max, T_step, output_file,
                                  gsigma, gdelta, fs, fo, fr, fd, b, c, n_nodes)
    """
    保存ForwardDiff计算结果到CSV文件
    
    参数:
    - results: 结果数组
    - μ_B: 重子化学势
    - T_min, T_max, T_step: 温度扫描参数
    - output_file: 输出文件路径
    - gsigma, gdelta, fs, fo, fr, fd, b, c, n_nodes: 模型参数
    """
    println("\n" * "="^60)
    println("保存结果到CSV文件")
    println("="^60)
    
    # 创建DataFrame
    df = DataFrame(results)
    
    # 确保输出目录存在
    output_dir = dirname(output_file)
    if !isdir(output_dir)
        mkpath(output_dir)
    end
    
    # 保存CSV文件 - 包含元数据头部
    try
        # 手动写入带有元数据的CSV文件
        open(output_file, "w") do io
            # 写入元数据头部（以#开头的注释行）
            println(io, "# ForwardDiff Temperature Scan Results")
            println(io, "# Generated on: $(Dates.now())")
            println(io, "# Model Parameters:")
            println(io, "# gsigma = $gsigma")
            println(io, "# gdelta = $gdelta")
            println(io, "# fs = $fs")
            println(io, "# fo = $fo")
            println(io, "# fr = $fr")
            println(io, "# fd = $fd")
            println(io, "# b = $b")
            println(io, "# c = $c")
            println(io, "# mu_B = $(μ_B*hc) MeV")
            println(io, "# T_range = $(T_min*hc) - $(T_max*hc) MeV")
            println(io, "# T_step = $(T_step*hc) MeV")
            println(io, "# nodes = $n_nodes")
            println(io, "#")
            
            # 写入CSV数据部分
            # 首先写入列名
            col_names = names(df)
            println(io, join(col_names, ","))
            
            # 然后写入数据行
            for row in eachrow(df)
                values = [string(row[col]) for col in col_names]
                println(io, join(values, ","))
            end
        end
        
        println("✓ 结果已保存到: $output_file")
        
        # 显示统计信息
        successful_points = sum(.!isnan.(df.P_T4))
        finite_κ3κ1 = sum(isfinite.(df.kappa3_over_kappa1))
        finite_κ4κ2 = sum(isfinite.(df.kappa4_over_kappa2))
        
        println("\n统计信息:")
        println("成功计算点数: $(successful_points)/$(nrow(df))")
        println("κ₃/κ₁ 有效值: $(finite_κ3κ1)/$(nrow(df))")
        println("κ₄/κ₂ 有效值: $(finite_κ4κ2)/$(nrow(df))")
        
        if finite_κ3κ1 > 0
            valid_κ3κ1 = filter(isfinite, df.kappa3_over_kappa1)
            println("κ₃/κ₁ 范围: $(round(minimum(valid_κ3κ1), digits=6)) 到 $(round(maximum(valid_κ3κ1), digits=6))")
        end
        
        if finite_κ4κ2 > 0
            valid_κ4κ2 = filter(isfinite, df.kappa4_over_kappa2)
            println("κ₄/κ₂ 范围: $(round(minimum(valid_κ4κ2), digits=6)) 到 $(round(maximum(valid_κ4κ2), digits=6))")
        end
        
    catch e
        println("✗ 保存CSV文件失败: $e")
    end
end