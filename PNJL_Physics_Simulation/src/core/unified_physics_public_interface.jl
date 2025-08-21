"""
统一物理计算公共接口

这个模块提供了适用于所有物理模型（PNJL、PNJL各向异性、旋转、气液相变等）的统一计算接口。
实现了自动微分、方程求解、相图计算等核心功能。

主要功能:
1. 自动微分接口 - 计算热力学势对各变量的偏导数
2. 方程组求解接口 - 使用nlsolve求解平衡态方程
3. 物理量计算接口 - 在给定解下计算各种物理量
4. 相图扫描接口 - 扫描T-μ或T-ρ相图并输出结果

作者: AI助手
日期: 2025年8月19日
"""
module UnifiedPhysicsPublicInterface

using ..ModelConfiguration
using ..UnifiedPhysicsInterface
using ..MathUtils
using ..AutodiffInterface
using ..EquationSolver: solve_equilibrium_equations

export 
    # Interface 1: 方程组求解接口  
    solve_equilibrium_equations,
    # Interface 2: 物理量计算接口
    calculate_physical_properties, PhysicalProperties,
    # Interface 3: 相图扫描接口
    scan_phase_diagram, PhasePoint, save_phase_diagram,
    # 重新导出自动微分接口
    compute_gradient, compute_hessian, check_equilibrium_conditions

# ============================================================================
# Interface 1: 通用方程组求解接口
# ============================================================================

"""
平衡态求解结果结构
"""
# EquilibriumSolution 已被移除，替换为独立文件 `equation_solver.jl` 中的纯解向量接口
include("solver_interface.jl")  # 现在将 solver_interface 作为求解器接口

"""
    solve_equilibrium_equations(equation_system, initial_guess, config::ModelConfig; options...)

使用nlsolve求解平衡态方程组 ∂Ω/∂xᵢ = 0

# 参数
- `equation_system`: 方程组函数，计算各变量的偏导数（残差）
- `initial_guess`: 初始猜测解向量
- `config`: 模型配置对象
- `options`: nlsolve的可选参数

# 返回
- `EquilibriumSolution` 对象，包含解、收敛状态等信息

# 示例
```julia
# 定义方程组（残差函数）
function equilibrium_equations(x, config)
    omega_func = variables -> calculate_omega_total(variables, config)
    return calculate_derivatives(omega_func, x, config)
end

# 求解
initial_guess = [0.1, 0.1, 0.2, 0.5, 0.5]
solution = solve_equilibrium_equations(equilibrium_equations, initial_guess, config)

if solution.converged
    println("找到平衡解: ", solution.solution)
else
    println("求解未收敛，残差: ", solution.residual_norm)
end
```
"""
# 低层实现已移至 `equation_solver.jl`，该文件提供了返回解向量的 `solve_equilibrium_equations`。
# 此处保留文档接口以便从热力学势直接调用（see below）。

"""
    solve_equilibrium_equations(thermodynamic_potential, initial_guess, config::ModelConfig; options...)

直接从热力学势函数求解平衡态（便捷接口）

# 参数
- `thermodynamic_potential`: 热力学势函数
- 其他参数同上述版本
"""
function solve_equilibrium_equations(thermodynamic_potential, initial_guess::Vector{T}, 
                                   config::ModelConfig; kwargs...) where T<:Real
    # 创建方程组函数，使用新的自动微分模块
    equation_system = (x, conf) -> compute_gradient(thermodynamic_potential, x, conf)

    # 调用低层求解器（返回解向量或在未收敛时抛出错误）
    return solve_equilibrium_equations(equation_system, initial_guess, config; kwargs...)
end

# ============================================================================
# Interface 3: 通用物理量计算接口
# ============================================================================

"""
物理性质结构，包含各种热力学量和观测量
"""
struct PhysicalProperties{T<:Real}
    # 基础热力学量
    pressure::T
    energy_density::T
    entropy_density::T
    baryon_density::T
    
    # 序参量
    chiral_condensates::Vector{T}    # 手征凝聚 [φᵤ, φₐ, φₛ]
    polyakov_loops::Tuple{T, T}      # Polyakov环 (Φ₁, Φ₂)
    
    # 敏感性系数
    compressibility::T               # 等温压缩系数
    heat_capacity::T                 # 比热容
    susceptibilities::Dict{String, T}  # 各种敏感性系数
    
    # 计算元数据
    temperature::T
    chemical_potentials::Vector{T}
    model_type::String
    calculation_time::Float64
end

"""
    calculate_physical_properties(solution_variables, config::ModelConfig, thermodynamic_potential; 
                                compute_susceptibilities=true)

在给定平衡解下计算各种物理量

# 参数
- `solution_variables`: 平衡态解向量 [φᵤ, φₐ, φₛ, Φ₁, Φ₂, ...]
- `config`: 模型配置对象
- `thermodynamic_potential`: 热力学势函数
- `compute_susceptibilities`: 是否计算敏感性系数

# 返回
- `PhysicalProperties` 对象，包含所有计算的物理量

# 示例
```julia
# 在平衡解下计算物理量
equilibrium = solve_equilibrium_equations(omega_func, initial_guess, config)
if equilibrium.converged
    properties = calculate_physical_properties(equilibrium.solution, config, omega_func)
    println("压强: $(properties.pressure)")
    println("能量密度: $(properties.energy_density)")
end
```
"""
function calculate_physical_properties(solution_variables::Vector{T}, config::ModelConfig,
                                     thermodynamic_potential; 
                                     compute_susceptibilities::Bool=true) where T<:Real
    start_time = time()
    
    try
        # 基础热力学量计算
        omega_value = thermodynamic_potential(solution_variables)
        pressure = -omega_value
        
        # 从解中提取序参量（假设前5个是标准变量）
        n_vars = length(solution_variables)
        if n_vars >= 5
            chiral_condensates = solution_variables[1:3]  # φᵤ, φₐ, φₛ
            polyakov_loops = (solution_variables[4], solution_variables[5])  # Φ₁, Φ₂
        else
            chiral_condensates = zeros(T, 3)
            polyakov_loops = (zero(T), zero(T))
        end
        
        # 通过自动微分计算其他热力学量
        # 能量密度: ε = -T²∂(Ω/T)/∂T
        energy_density = zero(T)  # 需要温度导数，暂时设为0
        
        # 熵密度: s = -∂Ω/∂T  
        entropy_density = zero(T)  # 需要温度导数，暂时设为0
        
        # 重子密度: ρ_B = -∂Ω/∂μ_B
        if hasfield(typeof(config), :chemical_potentials)
            # 通过化学势导数计算密度，使用新的自动微分模块
            try
                # 使用新的化学势导数计算接口
                dOmega_dmu = compute_chemical_potential_derivatives(
                    thermodynamic_potential, solution_variables, config.chemical_potentials, config
                )
                baryon_density = -dOmega_dmu[1]  # 第一个化学势对应重子密度
            catch
                baryon_density = zero(T)
            end
        else
            baryon_density = zero(T)
        end
        
        # 敏感性系数计算
        compressibility = zero(T)
        heat_capacity = zero(T)
        susceptibilities = Dict{String, T}()
        
        if compute_susceptibilities
            # 等温压缩系数: κ_T = -1/V (∂V/∂P)_T = 1/ρ (∂ρ/∂P)_T
            # 比热容: C_V = T(∂S/∂T)_V
            # 各种敏感性系数...
            susceptibilities["baryon"] = zero(T)
            susceptibilities["chiral"] = zero(T)
            susceptibilities["polyakov"] = zero(T)
        end
        
        calculation_time = time() - start_time
        
        model_type = string(typeof(config))
        
        return PhysicalProperties(
            pressure, energy_density, entropy_density, baryon_density,
            chiral_condensates, polyakov_loops,
            compressibility, heat_capacity, susceptibilities,
            config.temperature, config.chemical_potentials,
            model_type, calculation_time
        )
        
    catch e
        error("物理量计算失败: $e")
    end
end

# 辅助函数：更新化学势（需要针对不同模型实现）
function _update_chemical_potential(config::ModelConfig, new_mu::Real)
    # 这是一个简化的实现，实际需要针对每种模型配置进行特化
    return config  # 暂时返回原配置
end

# ============================================================================
# Interface 4: 通用相图扫描接口
# ============================================================================

"""
相图上的一个点
"""
struct PhasePoint{T<:Real}
    temperature::T
    chemical_potential::T        # 或密度，取决于扫描类型
    solution_variables::Vector{T}  # 平衡解
    properties::PhysicalProperties{T}  # 物理量
    converged::Bool             # 该点是否收敛
    scan_index::Tuple{Int, Int} # 在扫描网格中的索引
end

"""
    scan_phase_diagram(thermodynamic_potential, base_config::ModelConfig;
                      temperature_range=(0.1, 0.3), temperature_points=20,
                      chemical_potential_range=(0.0, 0.5), chemical_potential_points=20,
                      initial_guess_func=nothing, parallel=true, 
                      compute_properties=true, show_progress=true)

扫描T-μ相图并计算各点的平衡态和物理量

# 参数
- `thermodynamic_potential`: 热力学势函数
- `base_config`: 基础模型配置（温度和化学势会被扫描值覆盖）
- `temperature_range`: 温度扫描范围 (T_min, T_max)
- `temperature_points`: 温度点数
- `chemical_potential_range`: 化学势扫描范围 (μ_min, μ_max) 
- `chemical_potential_points`: 化学势点数
- `initial_guess_func`: 初始猜测函数 (T, μ) -> Vector，如为nothing则使用默认
- `parallel`: 是否使用并行计算
- `compute_properties`: 是否计算详细的物理量
- `show_progress`: 是否显示进度

# 返回
- `Vector{PhasePoint}`: 相图上所有点的结果

# 示例
```julia
# 扫描PNJL各向异性模型的T-μ相图
config = PNJLAnisoConfig()
omega_func = x -> calculate_omega_total(x, config)

phase_points = scan_phase_diagram(
    omega_func, config;
    temperature_range=(0.05, 0.25),
    temperature_points=15,
    chemical_potential_range=(0.0, 0.4), 
    chemical_potential_points=20
)

# 保存结果到文件
save_phase_diagram(phase_points, "pnjl_aniso_phase_diagram.dat")
```
"""
function scan_phase_diagram(thermodynamic_potential, base_config::ModelConfig;
                          temperature_range::Tuple{Float64, Float64}=(0.1, 0.3),
                          temperature_points::Int=20,
                          chemical_potential_range::Tuple{Float64, Float64}=(0.0, 0.5),
                          chemical_potential_points::Int=20,
                          initial_guess_func=nothing,
                          parallel::Bool=false,  # 暂时设为false，避免并行复杂性
                          compute_properties::Bool=true,
                          show_progress::Bool=true) where T<:Real
    
    println("🔍 开始相图扫描...")
    println("   温度范围: $(temperature_range) ($(temperature_points)点)")
    println("   化学势范围: $(chemical_potential_range) ($(chemical_potential_points)点)")
    println("   总计算点数: $(temperature_points * chemical_potential_points)")
    
    # 创建网格
    T_grid = range(temperature_range[1], temperature_range[2], length=temperature_points)
    mu_grid = range(chemical_potential_range[1], chemical_potential_range[2], length=chemical_potential_points)
    
    # 默认初始猜测函数
    if initial_guess_func === nothing
        initial_guess_func = (T, mu) -> [0.1, 0.1, 0.2, 0.5, 0.5]  # 标准5变量猜测
    end
    
    results = Vector{PhasePoint}()
    total_points = temperature_points * chemical_potential_points
    completed_points = 0
    
    start_time = time()
    
    for (i, T) in enumerate(T_grid)
        for (j, mu) in enumerate(mu_grid)
            # 更新配置
            current_config = _update_config_for_scan(base_config, T, mu)
            
            # 获取初始猜测
            initial_guess = initial_guess_func(T, mu)
            
            try
                # 求解平衡态（返回解向量或抛出错误）
                sol_vec = solve_equilibrium_equations(
                    thermodynamic_potential, initial_guess, current_config;
                    ftol=1e-8, show_trace=false
                )

                converged_flag = true

                # 计算物理量（仅在收敛时计算详细物理量）
                properties = if compute_properties && converged_flag
                    calculate_physical_properties(sol_vec, current_config, thermodynamic_potential)
                else
                    PhysicalProperties(0.0, 0.0, 0.0, 0.0, [0.0, 0.0, 0.0], (0.0, 0.0),
                                     0.0, 0.0, Dict{String, Float64}(), T, [mu], 
                                     string(typeof(current_config)), 0.0)
                end

                # 创建相图点
                point = PhasePoint(T, mu, sol_vec, properties, converged_flag, (i, j))
                push!(results, point)

            catch e
                # 如果计算失败或未收敛，创建未收敛的点
                empty_solution = zeros(length(initial_guess))
                empty_properties = PhysicalProperties(0.0, 0.0, 0.0, 0.0, [0.0, 0.0, 0.0], (0.0, 0.0),
                                                    0.0, 0.0, Dict{String, Float64}(), T, [mu],
                                                    string(typeof(current_config)), 0.0)

                point = PhasePoint(T, mu, empty_solution, empty_properties, false, (i, j))
                push!(results, point)

                if show_progress
                    println("   ⚠️  点 (T=$(round(T,digits=3)), μ=$(round(mu,digits=3))) 计算失败或未收敛: $e")
                end
            end
            
            completed_points += 1
            
            # 显示进度
            if show_progress && completed_points % max(1, div(total_points, 20)) == 0
                progress = completed_points / total_points * 100
                elapsed = time() - start_time
                estimated_total = elapsed / (completed_points / total_points)
                remaining = estimated_total - elapsed
                println("   📊 进度: $(round(progress, digits=1))% ($(completed_points)/$(total_points)), 预计剩余: $(round(remaining, digits=1))s")
            end
        end
    end
    
    total_time = time() - start_time
    converged_points = count(p -> p.converged, results)
    
    println("✅ 相图扫描完成!")
    println("   总用时: $(round(total_time, digits=2))s")
    println("   收敛点数: $(converged_points)/$(total_points) ($(round(converged_points/total_points*100, digits=1))%)")
    
    return results
end

"""
    save_phase_diagram(phase_points::Vector{PhasePoint}, filename::String; 
                      format=:dat, include_header=true)

将相图结果保存到文件

# 参数
- `phase_points`: 相图点数组
- `filename`: 输出文件名
- `format`: 输出格式 (:dat, :csv, :json)
- `include_header`: 是否包含表头

# 输出格式 (dat/csv)
列: T, μ, pressure, energy_density, chiral_condensate_u, chiral_condensate_d, chiral_condensate_s, 
    polyakov_1, polyakov_2, converged, ...
"""
function save_phase_diagram(phase_points::Vector{PhasePoint}, filename::String; 
                          format::Symbol=:dat, include_header::Bool=true)
    
    if format == :dat || format == :csv
        delimiter = format == :csv ? ',' : '\t'
        
        # 准备数据矩阵
        n_points = length(phase_points)
        if n_points == 0
            error("没有数据点可保存")
        end
        
        # 确定列数（基础列 + 变量数）
        n_variables = length(phase_points[1].solution_variables)
        n_cols = 10 + n_variables  # T, μ, P, ε, s, ρ, φ1,φ2,φ3, Φ1, Φ2, converged + variables
        
        data = Matrix{Float64}(undef, n_points, n_cols)
        
        for (i, point) in enumerate(phase_points)
            col = 1
            data[i, col] = point.temperature; col += 1
            data[i, col] = point.chemical_potential; col += 1  
            data[i, col] = point.properties.pressure; col += 1
            data[i, col] = point.properties.energy_density; col += 1
            data[i, col] = point.properties.entropy_density; col += 1
            data[i, col] = point.properties.baryon_density; col += 1
            
            # 手征凝聚
            for j in 1:3
                if length(point.properties.chiral_condensates) >= j
                    data[i, col] = point.properties.chiral_condensates[j]
                else
                    data[i, col] = 0.0
                end
                col += 1
            end
            
            # Polyakov环
            data[i, col] = point.properties.polyakov_loops[1]; col += 1
            data[i, col] = point.properties.polyakov_loops[2]; col += 1
            
            # 收敛状态 (1.0 = 收敛, 0.0 = 未收敛)
            data[i, col] = point.converged ? 1.0 : 0.0; col += 1
            
            # 解变量
            for j in 1:n_variables
                if col <= n_cols
                    data[i, col] = point.solution_variables[j]
                    col += 1
                end
            end
        end
        
        # 写入文件
        open(filename, "w") do io
            if include_header
                header_parts = ["T", "mu", "pressure", "energy_density", "entropy_density", "baryon_density",
                               "chiral_u", "chiral_d", "chiral_s", "polyakov_1", "polyakov_2", "converged"]
                for i in 1:n_variables
                    push!(header_parts, "var_$i")
                end
                
                if format == :dat
                    println(io, "# ", join(header_parts, "\t"))
                    # 省略生成时间以避免对 Dates 的依赖
                    println(io, "# Total points: ", n_points)
                    println(io, "# Converged points: ", count(p -> p.converged, phase_points))
                else
                    println(io, join(header_parts, ","))
                end
            end
            
            # 写入数据（手动实现，避免对 DelimitedFiles 的依赖问题）
            for i in 1:size(data, 1)
                row = join([string(data[i, j]) for j in 1:size(data, 2)], string(delimiter))
                write(io, row * "\n")
            end
        end
        
        println("📁 相图数据已保存至: $filename")
        println("   格式: $format")
        println("   数据点: $(n_points)")
        println("   列数: $(n_cols)")
        
    else
        error("不支持的输出格式: $format (支持: :dat, :csv)")
    end
end

# 辅助函数：为扫描更新配置
function _update_config_for_scan(config::ModelConfig, T::Float64, mu::Float64)
    # 这需要针对每种配置类型进行特化实现
    # 暂时返回原配置作为占位符
    return config
end

# 为不同模型配置特化更新函数
function _update_config_for_scan(config::PNJLAnisoConfig, T::Float64, mu::Float64)
    return PNJLAnisoConfig(
        cutoff=config.momentum_cutoff,
        n_p=config.n_momentum_points,
        n_theta=config.n_angle_points,
        T=T,
        mu=[mu, mu, mu],  # 假设所有味道的化学势相同
        Phi=config.polyakov_fields,
        xi=config.anisotropy_parameter
    )
end

function _update_config_for_scan(config::PNJLConfig, T::Float64, mu::Float64)
    return PNJLConfig(
        cutoff=config.momentum_cutoff,
        n_points=config.n_momentum_points,
        T=T,
        mu=[mu, mu, mu],
        Phi=config.polyakov_fields
    )
end

end  # module UnifiedPhysicsPublicInterface
