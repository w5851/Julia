"""
统一的物理计算接口

这个模块提供了所有物理模型的统一计算接口，简化了复杂计算的调用。
"""
module UnifiedPhysicsInterface

using ..IntegrationInterface
using ..ModelConfiguration
using ..PNJLFunctions
using ..MathUtils: safe_log, safe_exp

export PhysicsCalculation, calculate_thermodynamics, calculate_phase_diagram
export calculate_susceptibilities, calculate_fluctuations

"""
统一物理计算结果结构
"""
struct ThermodynamicResult{T<:Real}
    pressure::T
    energy_density::T
    entropy_density::T
    baryon_density::T
    chiral_condensate::Vector{T}
    polyakov_loop::Tuple{T, T}
    temperature::T
    chemical_potential::Vector{T}
end

"""
为PNJL模型计算热力学量
"""
function calculate_thermodynamics(config::PNJLConfig, phi::Vector{T}) where T<:Real
    grid = get_grid_config(config)
    masses = PNJLFunctions.calculate_mass_vec(phi)
    
    # 使用新的统一积分接口
    thermal_contribution = omega_thermal_integral(
        masses, 
        config.chemical_potentials, 
        config.temperature,
        config.polyakov_fields[1],
        config.polyakov_fields[2],
        grid
    )
    
    vacuum_contribution = vacuum_energy_integral(masses, grid)
    
    # 计算其他热力学量
    chiral_condensate = phi  # 简化，实际应该计算condensate
    polyakov_loop = config.polyakov_fields
    
    # 计算压强
    pressure = -(thermal_contribution + vacuum_contribution)
    
    # 其他量的计算...
    energy_density = 0.0  # 待实现
    entropy_density = 0.0  # 待实现  
    baryon_density = 0.0   # 待实现
    
    return ThermodynamicResult(
        pressure,
        energy_density,
        entropy_density, 
        baryon_density,
        chiral_condensate,
        polyakov_loop,
        config.temperature,
        config.chemical_potentials
    )
end

"""
为PNJL Aniso模型计算热力学量
"""
function calculate_thermodynamics(config::PNJLAnisoConfig, phi::Vector{T}) where T<:Real
    grids = get_grid_config(config)
    masses = PNJLFunctions.calculate_mass_vec(phi)
    
    # 使用2D积分接口
    method = GaussLegendreIntegration()
    
    # 计算各向异性的热力学积分
    thermal_integrand = function(p, theta)
        # 实现各向异性的被积函数
        E = sqrt(p^2 + masses[1]^2)  # 简化实现
        return safe_log(1 + safe_exp(-(E - config.chemical_potentials[1])/config.temperature))
    end
    
    thermal_contribution = integrate_2d(method, grids.momentum, grids.angle, thermal_integrand)
    thermal_contribution *= -config.temperature
    
    # 真空贡献
    vacuum_contribution = vacuum_energy_integral(collect(masses), grids.momentum)
    
    pressure = -(thermal_contribution + vacuum_contribution)
    
    return ThermodynamicResult(
        pressure, 0.0, 0.0, 0.0,
        phi, config.polyakov_fields,
        config.temperature, config.chemical_potentials
    )
end

"""
为Rotation模型计算热力学量
"""
function calculate_thermodynamics(config::RotationConfig, phi::Vector{T}) where T<:Real
    grids = get_grid_config(config)
    masses = PNJLFunctions.calculate_mass_vec(phi)
    
    method = GaussLegendreIntegration()
    
    # 包含旋转效应的积分
    thermal_integrand = function(p, l)
        E = sqrt(p^2 + masses[1]^2)  # 简化实现
        # 添加角动量量子化效应
        E_rot = E + l * config.angular_velocity
        return safe_log(1 + safe_exp(-(E_rot - config.chemical_potentials[1])/config.temperature))
    end
    
    thermal_contribution = integrate_2d(method, grids.momentum, grids.angular, thermal_integrand)
    thermal_contribution *= -config.temperature
    
    vacuum_contribution = vacuum_energy_integral(collect(masses), grids.momentum)
    
    pressure = -(thermal_contribution + vacuum_contribution)
    
    return ThermodynamicResult(
        pressure, 0.0, 0.0, 0.0,
        phi, (0.0, 0.0),  # Rotation model doesn't have Polyakov loops
        config.temperature, config.chemical_potentials
    )
end

"""
计算相图数据点
"""
function calculate_phase_diagram(model_type::Symbol, T_range::Vector{T}, 
                                mu_range::Vector{T}) where T<:Real
    results = Matrix{ThermodynamicResult{T}}(undef, length(T_range), length(mu_range))
    
    for (i, T) in enumerate(T_range)
        for (j, mu) in enumerate(mu_range)
            config = create_default_config(model_type)
            # 更新温度和化学势
            if hasfield(typeof(config), :temperature)
                config = @set config.temperature = T
            end
            if hasfield(typeof(config), :chemical_potentials)
                config = @set config.chemical_potentials = [mu, mu, mu]
            end
            
    # 简化的场参数，实际应该通过自洽求解得到
    phi = [-0.1, -0.1, -1.7]
            
            results[i, j] = calculate_thermodynamics(config, phi)
        end
    end
    
    return results
end

"""
计算热力学涨落和磁化率
"""
function calculate_susceptibilities(config::ModelConfig, phi::Vector{T}) where T<:Real
    # 通过有限差分计算各种磁化率
    epsilon = 1e-6
    
    base_result = calculate_thermodynamics(config, phi)
    
    # 重子数磁化率
    config_mu_plus = @set config.chemical_potentials[1] += epsilon
    result_mu_plus = calculate_thermodynamics(config_mu_plus, phi)
    
    chi_B = (result_mu_plus.baryon_density - base_result.baryon_density) / epsilon
    
    # 其他磁化率...
    
    return (chi_B=chi_B,)
end

end # module UnifiedPhysicsInterface
