"""
统一的物理模型配置管理

这个模块提供了所有物理模型的统一配置接口，简化模型使用和参数管理。
"""
module ModelConfiguration

using ..IntegrationInterface

export ModelConfig, PNJLConfig, PNJLAnisoConfig, RotationConfig, GasLiquidConfig
export create_default_config, get_grid_config, optimize_config

"""
抽象模型配置基类
"""
abstract type ModelConfig end

"""
PNJL模型配置
"""
struct PNJLConfig <: ModelConfig
    momentum_cutoff::Float64
    n_momentum_points::Int
    temperature::Float64
    chemical_potentials::Vector{Float64}
    polyakov_fields::Tuple{Float64, Float64}
    
    function PNJLConfig(;cutoff=10.0, n_points=64, T=0.15, 
                       mu=[0.32, 0.32, 0.32], Phi=(0.5, 0.5))
        new(cutoff, n_points, T, mu, Phi)
    end
end

"""
PNJL Anisotropic模型配置
"""
struct PNJLAnisoConfig <: ModelConfig
    momentum_cutoff::Float64
    n_momentum_points::Int
    n_angle_points::Int
    temperature::Float64
    chemical_potentials::Vector{Float64}
    polyakov_fields::Tuple{Float64, Float64}
    anisotropy_parameter::Float64
    
    function PNJLAnisoConfig(;cutoff=10.0, n_p=64, n_theta=32, T=0.15,
                           mu=[0.32, 0.32, 0.32], Phi=(0.5, 0.5), xi=0.1)
        new(cutoff, n_p, n_theta, T, mu, Phi, xi)
    end
end

"""
Rotation模型配置
"""
struct RotationConfig <: ModelConfig
    momentum_cutoff::Float64
    n_momentum_points::Int
    n_angular_points::Int
    temperature::Float64
    chemical_potentials::Vector{Float64}
    angular_velocity::Float64
    
    function RotationConfig(;cutoff=10.0, n_p=64, n_l=32, T=0.15,
                          mu=[0.32, 0.32, 0.32], Omega=0.05)
        new(cutoff, n_p, n_l, T, mu, Omega)
    end
end

"""
Gas-Liquid模型配置
"""
struct GasLiquidConfig <: ModelConfig
    momentum_cutoff::Float64
    n_momentum_points::Int
    temperature::Float64
    baryon_density::Float64
    asymmetry_parameter::Float64
    coupling_constants::NamedTuple
    
    function GasLiquidConfig(;cutoff=20.0, n_points=64, T=0.01, 
                           rho_B=0.15, Y=0.5, couplings=(fσ=10.0, fω=13.0, fρ=7.0, fδ=3.0, b=0.002, c=0.001))
        new(cutoff, n_points, T, rho_B, Y, couplings)
    end
end

"""
为指定模型创建默认配置
"""
function create_default_config(model_type::Symbol)
    if model_type == :PNJL
        return PNJLConfig()
    elseif model_type == :PNJL_aniso  
        return PNJLAnisoConfig()
    elseif model_type == :Rotation
        return RotationConfig()
    elseif model_type == :GasLiquid
        return GasLiquidConfig()
    else
        error("Unknown model type: $model_type")
    end
end

"""
获取积分网格配置
"""
function get_grid_config(config::PNJLConfig)
    return create_momentum_grid(config.n_momentum_points, config.momentum_cutoff)
end

function get_grid_config(config::PNJLAnisoConfig)
    p_grid = create_momentum_grid(config.n_momentum_points, config.momentum_cutoff)
    theta_grid = create_angle_grid(config.n_angle_points)
    return (momentum=p_grid, angle=theta_grid)
end

function get_grid_config(config::RotationConfig)
    p_grid = create_momentum_grid(config.n_momentum_points, config.momentum_cutoff)
    l_grid = create_angular_momentum_grid(config.n_angular_points)
    return (momentum=p_grid, angular=l_grid)
end

function get_grid_config(config::GasLiquidConfig)
    return create_momentum_grid(config.n_momentum_points, config.momentum_cutoff)
end

"""
优化配置参数以平衡精度和性能
"""
function optimize_config(config::ModelConfig, target_precision::Float64)
    # 基于目标精度调整积分点数
    # 这里可以实现自适应网格优化算法
    return config  # 简化实现，返回原配置
end

end # module ModelConfiguration
