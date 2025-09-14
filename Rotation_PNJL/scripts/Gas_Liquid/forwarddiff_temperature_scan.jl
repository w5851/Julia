#!/usr/bin/env julia
"""
使用ForwardDiff方法进行温度扫描
基于test_forwarddiff_high_order.jl的逻辑，扫描温度从20-200 MeV
"""

using NLsolve
using ForwardDiff
using CSV
using DataFrames
using Dates

# 加载PNJL模型的核心函数
include("test_fortran_exact_derivative.jl")

function NewQuark_mu_pnjl_fixed(x_array, T, mu_B, model_params)
    """
    修正的PNJL模型类型保护约束方程求解器
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

# 修正的高阶导数函数，使用闭包方式避免类型问题
function create_pressure_function(gsigma, gdelta, T, model_params)
    """创建压强函数的闭包，避免类型推断问题"""
    return μ -> begin
        # 每次调用时在lambda内部创建新数组，确保类型一致性
        X0_new = [gsigma, gdelta, μ/2.0, μ/2.0]
        SolveOmega_pnjl_fixed(X0_new, T, μ, model_params) / T^4
    end
end

function D1_Pressure_mu(gsigma, gdelta, T, mu_B, model_params)
    """一阶导数：∂(P/T⁴)/∂μ_B"""
    pressure_func = create_pressure_function(gsigma, gdelta, T, model_params)
    return ForwardDiff.derivative(pressure_func, mu_B)
end

function D2_Pressure_mu(gsigma, gdelta, T, mu_B, model_params)
    """二阶导数：∂²(P/T⁴)/∂μ_B²"""
    d1_func = μ -> D1_Pressure_mu(gsigma, gdelta, T, μ, model_params)
    return ForwardDiff.derivative(d1_func, mu_B)
end

function D3_Pressure_mu(gsigma, gdelta, T, mu_B, model_params)
    """三阶导数：∂³(P/T⁴)/∂μ_B³"""
    d2_func = μ -> D2_Pressure_mu(gsigma, gdelta, T, μ, model_params)
    return ForwardDiff.derivative(d2_func, mu_B)
end

function D4_Pressure_mu_enhanced(gsigma, gdelta, T, mu_B, model_params)
    """增强的四阶导数计算：使用三阶导数的数值微分"""
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

function forwarddiff_temperature_scan(μ_B, T_min, T_max, T_step, output_file)
    """
    使用ForwardDiff方法进行温度扫描
    """
    
    # 设置模型参数
    nodes = get_nodes(256)
    gsigma = 1.25
    gdelta = 0.01
    fs = 17.28476
    fo = 11.66174
    fr = 0.89363
    fd = 0.0
    b = 0.00210
    c = -0.00297
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
    
    # 预分配结果数组
    results = []
    
    # 固定初始猜测值
    mu_p = μ_B / 2.0
    mu_n = μ_B / 2.0
    X0 = [gsigma, gdelta, mu_p, mu_n]
    
    # 循环计算每个温度点
    for (i, T) in enumerate(T_array)
        println("\n$(i)/$(n_temps): 计算 T = $(round(T*hc, digits=1)) MeV")
        
        try
            # 1. 验证基础计算
            pressure = SolveOmega_pnjl_fixed(X0, T, μ_B, model_params)
            pressure_normalized = pressure / T^4
            
            # 2. 计算一到四阶导数
            d1_raw = D1_Pressure_mu(gsigma, gdelta, T, μ_B, model_params)
            d1_norm = d1_raw * T  # 转换为对μ/T的导数
            
            d2_raw = D2_Pressure_mu(gsigma, gdelta, T, μ_B, model_params)
            d2_norm = d2_raw * T^2  # 转换为对(μ/T)²的导数
            
            d3_raw = D3_Pressure_mu(gsigma, gdelta, T, μ_B, model_params)
            d3_norm = d3_raw * T^3  # 转换为对(μ/T)³的导数
            
            d4_raw = D4_Pressure_mu_enhanced(gsigma, gdelta, T, μ_B, model_params)
            d4_norm = d4_raw * T^4  # 转换为对(μ/T)⁴的导数
            
            # 3. 计算热力学涨落量比值
            κ1 = d1_norm
            κ2 = d2_norm  
            κ3 = d3_norm
            κ4 = d4_norm
            
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
            
            # 4. 存储结果
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
    
    # 5. 保存结果到CSV文件
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
            println(io, "# nodes = $(length(nodes))")
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
    
    return df
end

# 主程序
println("开始ForwardDiff温度扫描计算")

# 设置扫描参数
μ_B_fixed = 697.0/hc  # 固定重子化学势 697 MeV
T_min = 20.0/hc       # 最小温度 20 MeV
T_max = 200.0/hc      # 最大温度 200 MeV
T_step = 1.0/hc       # 温度步长 1 MeV

# 设置输出文件路径
output_file = joinpath(@__DIR__, "..", "..", "output", "Gas_Liquid", "forwarddiff_temperature_scan.csv")

println("扫描参数:")
println("μ_B = $(μ_B_fixed*hc) MeV (固定)")
println("温度范围: $(T_min*hc) - $(T_max*hc) MeV")
println("步长: $(T_step*hc) MeV")
println("输出文件: $output_file")

# 执行温度扫描
try
    results_df = forwarddiff_temperature_scan(μ_B_fixed, T_min, T_max, T_step, output_file)
    
    println("\n" * "="^60)
    println("🎉 ForwardDiff温度扫描完成!")
    println("="^60)
    
catch e
    println("\n" * "="^60)
    println("✗ ForwardDiff温度扫描失败: $e")
    println("错误详情: ", sprint(showerror, e))
    println("="^60)
end
