"""
PNJL各向异性模型专用积分工具

本模块包含PNJL各向异性模型特定的积分函数和工具，
使用核心积分接口提供高层次的物理计算接口。
"""
module PNJLAnisoIntegrationUtils

using ..IntegrationInterface: GaussLegendreIntegration, integrate_2d, MomentumGrid, AngleGrid
using ..MathUtils: safe_log

export omega_thermal_integral_aniso, vacuum_energy_integral_aniso

# ============================================================================
# PNJL各向异性模型专用积分函数
# ============================================================================

"""
    omega_thermal_integral_aniso(masses::Vector, mu::Vector, T::Float64,
                                Phi1::Float64, Phi2::Float64,
                                p_grid::MomentumGrid, t_grid::AngleGrid,
                                xi::Float64=0.0, method=GaussLegendreIntegration())

计算各向异性PNJL模型中的热力学Omega贡献积分。

# Arguments
- `masses`: 有效夸克质量向量
- `mu`: 化学势向量
- `T`: 温度
- `Phi1, Phi2`: Polyakov loop参数
- `p_grid`: 动量积分网格
- `t_grid`: 角度积分网格
- `xi`: 各向异性参数
- `method`: 积分方法

# Returns
- 热力学Omega积分贡献

# Physics
包含各向异性效应的能量：E(p,t) = √(p² + m² + ξ(pt)²)
"""
function omega_thermal_integral_aniso(masses::Vector, mu::Vector, T::Float64,
                                    Phi1::Float64, Phi2::Float64,
                                    p_grid::MomentumGrid, t_grid::AngleGrid, 
                                    xi::Float64=0.0, method=GaussLegendreIntegration())
    
    total_contribution = 0.0
    
    # 计算每个夸克味的贡献
    @inbounds for (i, mass_i) in enumerate(masses)
        mu_i = mu[i]
        
        # 定义各向异性被积函数
        integrand = function(p, t; mass=mass_i, chemical_potential=mu_i, temperature=T,
                           polyakov1=Phi1, polyakov2=Phi2, anisotropy=xi)
            # 各向异性能量
            E_aniso = sqrt(p^2 + mass^2 + anisotropy * (p*t)^2)
            
            # Polyakov loop对数项
            log_term = _calculate_polyakov_log_term_aniso(E_aniso, chemical_potential, temperature, 
                                                        polyakov1, polyakov2)
            
            return p^2 * log_term  # 包含球坐标系数
        end
        
        # 执行2D积分
        contribution = integrate_2d(method, p_grid, t_grid, integrand;
                                  mass=mass_i, chemical_potential=mu_i, temperature=T,
                                  polyakov1=Phi1, polyakov2=Phi2, anisotropy=xi)
        total_contribution += contribution
    end
    
    # 包含温度和几何因子
    return total_contribution * (-T) / (3.0 * π^2)
end

"""
    vacuum_energy_integral_aniso(masses::Vector, p_grid::MomentumGrid, t_grid::AngleGrid,
                                xi::Float64=0.0, method=GaussLegendreIntegration())

计算各向异性模型的真空能量积分。

# Arguments
- `masses`: 有效夸克质量向量
- `p_grid`: 动量积分网格
- `t_grid`: 角度积分网格
- `xi`: 各向异性参数
- `method`: 积分方法

# Returns
- 真空能量积分结果
"""
function vacuum_energy_integral_aniso(masses::Vector, p_grid::MomentumGrid, t_grid::AngleGrid,
                                    xi::Float64=0.0, method=GaussLegendreIntegration())
    
    total_energy = 0.0
    
    @inbounds for mass_i in masses
        # 定义各向异性真空能被积函数
        integrand = function(p, t; mass=mass_i, anisotropy=xi)
            E_aniso = sqrt(p^2 + mass^2 + anisotropy * (p*t)^2)
            return p^2 * E_aniso
        end
        
        # 执行2D积分
        contribution = integrate_2d(method, p_grid, t_grid, integrand;
                                  mass=mass_i, anisotropy=xi)
        total_energy += contribution
    end
    
    # 包含几何因子和颜色因子
    return total_energy * (-3) / (3.0 * π^2)  # -N_c factor
end

# ============================================================================
# 各向异性模型专用辅助函数
# ============================================================================

"""
    _calculate_polyakov_log_term_aniso(E, mu, temperature, Phi1, Phi2)

计算各向异性模型中包含Polyakov loop效应的对数项。
"""
function _calculate_polyakov_log_term_aniso(E::Float64, mu::Float64, temperature::Float64, 
                                          Phi1::Float64, Phi2::Float64)
    # 计算费米子和反费米子的能量参数
    x = (E - mu) / temperature
    x_anti = (E + mu) / temperature
    
    # 费米子项（完整的3色群论展开）
    f_term = safe_log(1.0 + 3.0*Phi1*exp(-x) + 3.0*Phi2*exp(-2*x) + exp(-3*x))
    
    # 反费米子项
    f_anti_term = safe_log(1.0 + 3.0*Phi2*exp(-x_anti) + 3.0*Phi1*exp(-2*x_anti) + exp(-3*x_anti))
    
    return f_term + f_anti_term
end

end # module PNJLAnisoIntegrationUtils
