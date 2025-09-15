# Advanced_FindTforDiff_Bayesian.jl
# 基于BayesianOptimization.jl的智能温度查找模块
# 使用贝叶斯优化代替等间距扫描，显著提高搜索效率
# 依赖于 Advanced_ForwardDiff.jl 中的核心计算函数

include("Advanced_ForwardDiff.jl")
using BayesianOptimization
using GaussianProcesses

function find_temperature_bayesian(target_kappa3_kappa1, target_kappa4_kappa2, μ_B, 
                                  T_min, T_max;
                                  gsigma=1.25, gdelta=0.01,
                                  fs=17.28476, fo=11.66174, fr=0.89363, fd=0.0,
                                  b=0.00210, c=-0.00297, n_nodes=256,
                                  max_iterations=30, verbose=true)
    """
    使用贝叶斯优化寻找给定κ₃/κ₁和κ₄/κ₂值对应的温度
    
    参数:
    - target_kappa3_kappa1: 目标κ₃/κ₁值
    - target_kappa4_kappa2: 目标κ₄/κ₂值
    - μ_B: 重子化学势
    - T_min: 温度搜索下限
    - T_max: 温度搜索上限
    - 其他参数: 模型参数
    - max_iterations: 贝叶斯优化最大迭代次数
    - verbose: 是否打印详细信息
    
    返回:
    - (T_kappa3_kappa1, T_kappa4_kappa2): 对应的温度值
    """
    
    if verbose
        println("="^60)
        println("贝叶斯优化温度查找")
        println("="^60)
        println("目标值:")
        println("  κ₃/κ₁ = $target_kappa3_kappa1")
        println("  κ₄/κ₂ = $target_kappa4_kappa2")
        println("搜索范围:")
        println("  μ_B = $(μ_B*hc) MeV")
        println("  T ∈ [$(T_min*hc), $(T_max*hc)] MeV")
        println("  最大迭代次数: $max_iterations")
    end
    
    # 设置模型参数
    nodes = get_nodes(n_nodes)
    couplings = [fs, fo, fr, fd, b, c]
    model_params = (nodes, couplings)
    
    # 创建温度计算函数（返回κ值）
    function calculate_kappa_ratios(T)
        try
            κ1, κ2, κ3, κ4, κ3_κ1, κ4_κ2 = calculate_forwarddiff_thermodynamic_fluctuations(
                gsigma, gdelta, T, μ_B, model_params)
            
            if isfinite(κ3_κ1) && isfinite(κ4_κ2)
                return κ3_κ1, κ4_κ2
            else
                return NaN, NaN
            end
        catch e
            return NaN, NaN
        end
    end
    
    # 定义目标函数（最小化距离）
    function objective_kappa3_kappa1(T)
        κ3_κ1, _ = calculate_kappa_ratios(T)
        if isfinite(κ3_κ1)
            return -(target_kappa3_kappa1 - κ3_κ1)^2  # 负号因为BayesianOptimization最大化
        else
            return -1e6  # 惩罚无效值
        end
    end
    
    function objective_kappa4_kappa2(T)
        _, κ4_κ2 = calculate_kappa_ratios(T)
        if isfinite(κ4_κ2)
            return -(target_kappa4_kappa2 - κ4_κ2)^2  # 负号因为BayesianOptimization最大化
        else
            return -1e6  # 惩罚无效值
        end
    end
    
    # 寻找κ₃/κ₁对应的温度
    if verbose
        println("\n第一步：贝叶斯优化寻找κ₃/κ₁ = $target_kappa3_kappa1 对应的温度...")
    end
    
    T_kappa3_kappa1 = bayesian_optimize_temperature(
        objective_kappa3_kappa1, T_min, T_max, max_iterations, verbose, "κ₃/κ₁")
    
    # 寻找κ₄/κ₂对应的温度
    if verbose
        println("\n第二步：贝叶斯优化寻找κ₄/κ₂ = $target_kappa4_kappa2 对应的温度...")
    end
    
    T_kappa4_kappa2 = bayesian_optimize_temperature(
        objective_kappa4_kappa2, T_min, T_max, max_iterations, verbose, "κ₄/κ₂")
    
    if verbose
        println("\n" * "="^60)
        println("贝叶斯优化结果:")
        println("  κ₃/κ₁ = $target_kappa3_kappa1 → T = $(round(T_kappa3_kappa1*hc, digits=2)) MeV")
        println("  κ₄/κ₂ = $target_kappa4_kappa2 → T = $(round(T_kappa4_kappa2*hc, digits=2)) MeV")
        println("="^60)
    end
    
    return T_kappa3_kappa1, T_kappa4_kappa2
end

function bayesian_optimize_temperature(objective_func, T_min, T_max, max_iterations, verbose, target_name)
    """
    使用贝叶斯优化寻找单个目标函数的最优温度
    """
    
    # 设置高斯过程
    kern = SEArd([1.0], 5.0)  # Squared Exponential kernel
    gp = GP(Float64[], Float64[], MeanZero(), kern, -2.0)
    
    # 设置贝叶斯优化器
    model = BOpt(objective_func,
                 ElType = Float64,
                 GP = gp,
                 sense = Max,  # 最大化目标函数
                 verbosity = verbose ? Progress : Silent,
                 maxiterations = max_iterations,
                 initializer_iterations = 5,
                 acquisitionoptions = (method = :UpperConfidenceBound,
                                     beta = 2.0))
    
    # 设置搜索域
    boptimize!(model, UpperConfidenceBound(),
              [(T_min, T_max)], 
              initializer = CLHSampling())
    
    # 获取最优温度
    best_T = model.opt[1]
    best_value = model.opt[2]
    
    if verbose
        println("  $target_name 优化完成:")
        println("    最优温度: $(round(best_T*hc, digits=2)) MeV")
        println("    目标函数值: $(round(-best_value, digits=6))")  # 负号转回距离
        println("    总评估次数: $(length(model.GP.x))")
    end
    
    return best_T
end

function compare_bayesian_vs_grid_search(target_kappa3_kappa1, target_kappa4_kappa2, μ_B, 
                                        T_min, T_max;
                                        grid_points=50, max_bayesian_iter=30,
                                        gsigma=1.25, gdelta=0.01,
                                        fs=17.28476, fo=11.66174, fr=0.89363, fd=0.0,
                                        b=0.00210, c=-0.00297, n_nodes=256)
    """
    比较贝叶斯优化与网格搜索的性能
    
    返回:
    - (bayesian_results, grid_results, comparison): 贝叶斯结果、网格结果、性能对比
    """
    
    println("="^80)
    println("贝叶斯优化 vs 网格搜索 性能对比")
    println("="^80)
    
    # 贝叶斯优化
    println("\n🚀 贝叶斯优化方法:")
    bayesian_start = time()
    T_bay_k3, T_bay_k4 = find_temperature_bayesian(
        target_kappa3_kappa1, target_kappa4_kappa2, μ_B, T_min, T_max;
        max_iterations=max_bayesian_iter, verbose=true,
        gsigma=gsigma, gdelta=gdelta, fs=fs, fo=fo, fr=fr, fd=fd, b=b, c=c, n_nodes=n_nodes)
    bayesian_time = time() - bayesian_start
    
    # 网格搜索（简化版本，使用等间距扫描）
    println("\n📊 网格搜索方法:")
    grid_start = time()
    T_grid_k3, T_grid_k4 = find_temperature_grid_search(
        target_kappa3_kappa1, target_kappa4_kappa2, μ_B, T_min, T_max, grid_points;
        gsigma=gsigma, gdelta=gdelta, fs=fs, fo=fo, fr=fr, fd=fd, b=b, c=c, n_nodes=n_nodes)
    grid_time = time() - grid_start
    
    # 性能对比
    println("\n" * "="^80)
    println("性能对比结果:")
    println("="^80)
    
    println("计算时间:")
    println("  贝叶斯优化: $(round(bayesian_time, digits=2)) 秒")
    println("  网格搜索:   $(round(grid_time, digits=2)) 秒")
    println("  加速比:     $(round(grid_time/bayesian_time, digits=2))x")
    
    println("\n计算次数:")
    println("  贝叶斯优化: ~$(max_bayesian_iter*2) 次函数评估")
    println("  网格搜索:   $grid_points 次函数评估")
    
    println("\n温度结果对比:")
    println("  κ₃/κ₁目标值: $target_kappa3_kappa1")
    println("    贝叶斯: $(round(T_bay_k3*hc, digits=2)) MeV")
    println("    网格:   $(round(T_grid_k3*hc, digits=2)) MeV")
    println("    差异:   $(round(abs(T_bay_k3-T_grid_k3)*hc, digits=2)) MeV")
    
    println("  κ₄/κ₂目标值: $target_kappa4_kappa2")
    println("    贝叶斯: $(round(T_bay_k4*hc, digits=2)) MeV")
    println("    网格:   $(round(T_grid_k4*hc, digits=2)) MeV")
    println("    差异:   $(round(abs(T_bay_k4-T_grid_k4)*hc, digits=2)) MeV")
    
    bayesian_results = (T_kappa3_kappa1=T_bay_k3, T_kappa4_kappa2=T_bay_k4, time=bayesian_time)
    grid_results = (T_kappa3_kappa1=T_grid_k3, T_kappa4_kappa2=T_grid_k4, time=grid_time)
    comparison = (speedup=grid_time/bayesian_time, 
                 accuracy_k3=abs(T_bay_k3-T_grid_k3)*hc,
                 accuracy_k4=abs(T_bay_k4-T_grid_k4)*hc)
    
    return bayesian_results, grid_results, comparison
end

function find_temperature_grid_search(target_kappa3_kappa1, target_kappa4_kappa2, μ_B, 
                                     T_min, T_max, n_points;
                                     gsigma=1.25, gdelta=0.01,
                                     fs=17.28476, fo=11.66174, fr=0.89363, fd=0.0,
                                     b=0.00210, c=-0.00297, n_nodes=256)
    """
    网格搜索方法（用于对比）
    """
    
    # 设置模型参数
    nodes = get_nodes(n_nodes)
    couplings = [fs, fo, fr, fd, b, c]
    model_params = (nodes, couplings)
    
    T_array = range(T_min, T_max, length=n_points)
    kappa3_kappa1_array = Float64[]
    kappa4_kappa2_array = Float64[]
    T_valid_array = Float64[]
    
    println("  网格搜索进行中...")
    
    for T in T_array
        try
            κ1, κ2, κ3, κ4, κ3_κ1, κ4_κ2 = calculate_forwarddiff_thermodynamic_fluctuations(
                gsigma, gdelta, T, μ_B, model_params)
            
            if isfinite(κ3_κ1) && isfinite(κ4_κ2)
                push!(kappa3_kappa1_array, κ3_κ1)
                push!(kappa4_kappa2_array, κ4_κ2)
                push!(T_valid_array, T)
            end
        catch e
            # 跳过计算失败的点
        end
    end
    
    # 简单的最近邻搜索
    function find_closest_temperature(target_value, kappa_array, T_array)
        min_dist = Inf
        best_T = T_array[1]
        
        for (i, kappa_val) in enumerate(kappa_array)
            dist = abs(target_value - kappa_val)
            if dist < min_dist
                min_dist = dist
                best_T = T_array[i]
            end
        end
        
        return best_T
    end
    
    T_kappa3_kappa1 = find_closest_temperature(target_kappa3_kappa1, kappa3_kappa1_array, T_valid_array)
    T_kappa4_kappa2 = find_closest_temperature(target_kappa4_kappa2, kappa4_kappa2_array, T_valid_array)
    
    println("  网格搜索完成，共计算 $(length(T_valid_array)) 个有效点")
    
    return T_kappa3_kappa1, T_kappa4_kappa2
end

# 导出主要函数
export find_temperature_bayesian, compare_bayesian_vs_grid_search