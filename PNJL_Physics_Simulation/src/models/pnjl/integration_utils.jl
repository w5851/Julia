"""
PNJL模型专用积分工具

本模块包含PNJL物理模型特定的积分函数和工具，
使用核心积分接口提供高层次的物理计算接口。
"""
module PNJLIntegrationUtils

using ..IntegrationInterface: GaussLegendreIntegration, integrate, MomentumGrid
using ..MathUtils: safe_log

export omega_thermal_integral, vacuum_energy_integral

# ============================================================================
# PNJL模型专用积分函数
# ============================================================================

"""
    omega_thermal_integral(masses::Vector{T}, mu::Vector{T}, temperature::T, 
                          Phi1::T, Phi2::T, grid::MomentumGrid, 
                          method::IntegrationMethod=GaussLegendreIntegration()) where T<:Real -> T

计算PNJL模型中的热力学Omega贡献积分。

这个函数替代了原来分散在各模型中的 `calculate_log_sum` 函数，
提供统一的热力学积分计算接口。

# Arguments
- `masses::Vector{T}`: 有效夸克质量 [m_u, m_d, m_s]
- `mu::Vector{T}`: 化学势 [μ_u, μ_d, μ_s]
- `temperature::T`: 温度
- `Phi1, Phi2::T`: Polyakov loop参数
- `grid::MomentumGrid`: 动量积分网格
- `method::IntegrationMethod`: 积分方法

# Returns
- `T`: 热力学Omega积分贡献

# Physics
计算公式：Ω_thermal = -T ∑ᵢ ∫ dp p² log[1 + exp(-E/T)] * Polyakov_factor
其中 E = √(p² + mᵢ²)，Polyakov_factor包含Polyakov loop的影响。

# Example
```julia
using .PNJLIntegrationUtils
using .IntegrationInterface

masses = [0.3, 0.3, 0.5]  # GeV
mu = [0.0, 0.0, 0.0]      # GeV  
T = 0.15                  # GeV
Phi1, Phi2 = 0.5, 0.5
grid = create_momentum_grid(64, 2.0)

result = omega_thermal_integral(masses, mu, T, Phi1, Phi2, grid)
```
"""
function omega_thermal_integral(masses::Vector{T}, mu::Vector{T}, temperature::T,
                               Phi1::T, Phi2::T, grid::MomentumGrid,
                               method=GaussLegendreIntegration()) where T<:Real
    
    total_contribution = zero(T)
    
    # 对每种夸克味道计算积分贡献
    @inbounds for (i, mass_i) in enumerate(masses)
        mu_i = mu[i]
        
        # 定义被积函数：包含能量、分布函数和Polyakov loop效应
        integrand = function(p; mass=mass_i, chemical_potential=mu_i, temp=temperature, 
                           polyakov1=Phi1, polyakov2=Phi2)
            E = sqrt(p^2 + mass^2)  # 夸克能量
            
            # 计算包含Polyakov loop效应的对数项
            log_term = _calculate_polyakov_log_term(E, chemical_potential, temp, polyakov1, polyakov2)
            
            # 返回 p²log_term 作为被积函数
            return p^2 * log_term
        end
        
        # 执行积分（使用关键字参数传递）
        contribution = integrate(method, grid, integrand; 
                               mass=mass_i, chemical_potential=mu_i, temp=temperature,
                               polyakov1=Phi1, polyakov2=Phi2)
        total_contribution += contribution
    end
    
    # 返回热力学贡献，包含温度因子和几何因子
    return total_contribution * (-temperature) / (3.0 * π^2)
end

"""
    vacuum_energy_integral(masses::Vector{T}, grid::MomentumGrid, 
                          method=GaussLegendreIntegration()) where T<:Real -> T

计算真空能量积分贡献。

替代原来的 `calculate_energy_sum` 函数，提供统一的真空能量计算。

# Arguments
- `masses::Vector{T}`: 有效夸克质量
- `grid::MomentumGrid`: 动量积分网格
- `method::IntegrationMethod`: 积分方法

# Returns
- `T`: 真空能量积分结果

# Physics
计算公式：E_vacuum = ∑ᵢ ∫ dp p² √(p² + mᵢ²)

# Example
```julia
masses = [0.3, 0.3, 0.5]  # GeV
grid = create_momentum_grid(64, 2.0)
vacuum_energy = vacuum_energy_integral(masses, grid)
```
"""
function vacuum_energy_integral(masses::Vector{T}, grid::MomentumGrid,
                               method=GaussLegendreIntegration()) where T<:Real
    
    total_energy = zero(T)
    
    @inbounds for mass_i in masses
        # 定义被积函数：p² * E(p)
        integrand = function(p; mass=mass_i)
            E = sqrt(p^2 + mass^2)
            return p^2 * E
        end
        
        # 执行积分
        energy_contribution = integrate(method, grid, integrand; mass=mass_i)
        total_energy += energy_contribution
    end
    
    # 返回真空能量，包含几何因子
    return total_energy / (3.0 * π^2)
end

# ============================================================================
# PNJL模型专用辅助函数
# ============================================================================

"""
    _calculate_polyakov_log_term(E, mu, temperature, Phi1, Phi2)

计算包含Polyakov loop效应的对数项。

这是内部辅助函数，实现PNJL模型中费米子分布函数与Polyakov loop的耦合。
"""
function _calculate_polyakov_log_term(E::T, mu::T, temperature::T, Phi1::T, Phi2::T) where T<:Real
    # 计算费米子和反费米子的能量参数
    x = (E - mu) / temperature
    x_anti = (E + mu) / temperature
    
    # 使用安全对数函数避免数值问题
    term1 = safe_log(1 + Phi1 * exp(-x) + Phi2 * exp(-2*x))
    term2 = safe_log(1 + Phi2 * exp(-x_anti) + Phi1 * exp(-2*x_anti))
    
    return term1 + term2
end

end # module PNJLIntegrationUtils
