"""
旋转PNJL模型中求和接口的使用示例

展示如何使用新的求和接口函数来正确处理角动量量子化。
"""

using StaticArrays
using ..IntegrationInterface: GaussLegendreIntegration, angular_momentum_sum, 
                             discrete_sum, mixed_integral_sum,
                             create_momentum_grid
using ..MathUtils: safe_log

# 示例：旋转模型的热力学积分
function rotation_thermal_integral_example(masses::Vector, mu::Vector, T::Float64,
                                         Phi1::Float64, Phi2::Float64, omega::Float64,
                                         p_nodes, p_weights, n_max::Int)
    """
    使用新求和接口计算旋转模型的热力学积分
    
    物理公式：Ω_thermal = -T ∫ dp Σₙ (2n+1) log[ℱ(E(p,n))]
    其中 E(p,n) = sqrt(p² + m² + (ωn/c)²) 是旋转系统中的能量
    """
    
    # 创建动量积分网格
    p_grid = create_momentum_grid(length(p_nodes), maximum(p_nodes))
    method = GaussLegendreIntegration()
    
    total_contribution = 0.0
    
    # 对每个夸克味道计算
    @inbounds for (i, mass_i) in enumerate(masses)
        mu_i = mu[i]
        
        # 定义被积函数（关于动量p和角动量量子数n）
        integrand = function(p, n)
            # 旋转系统中的能量
            E_rot = sqrt(p^2 + mass_i^2 + (omega * n)^2)
            
            # 包含Polyakov loop的对数项
            log_term = _calculate_rotation_log_term(E_rot, mu_i, T, Phi1, Phi2)
            
            return p^2 * log_term  # 包含球坐标中的p²因子
        end
        
        # 使用角动量求和函数（自动包含(2n+1)权重）
        contribution = angular_momentum_sum(method, p_grid, n_max, integrand)
        total_contribution += contribution
    end
    
    # 包含温度因子和几何因子
    return total_contribution * (-T) / (3.0 * π^2)
end

# 辅助函数：计算旋转系统的Polyakov loop对数项
function _calculate_rotation_log_term(E_rot::Float64, mu::Float64, T::Float64, 
                                    Phi1::Float64, Phi2::Float64)
    """计算旋转系统中包含Polyakov loop的对数项"""
    x = (E_rot - mu) / T
    x_anti = (E_rot + mu) / T
    
    # 费米子项
    f_term = safe_log(1 + 3*Phi1*exp(-x) + 3*Phi2*exp(-2*x) + exp(-3*x))
    # 反费米子项  
    f_anti_term = safe_log(1 + 3*Phi2*exp(-x_anti) + 3*Phi1*exp(-2*x_anti) + exp(-3*x_anti))
    
    return f_term + f_anti_term
end

# 示例：纯离散求和（不涉及积分）
function discrete_sum_example()
    """
    展示纯离散求和的用法
    
    计算：Σₙ₌₀¹⁰ (2n+1) exp(-n²/T)
    """
    T = 1.0
    indices = 0:10
    
    result = discrete_sum(indices) do n
        statistical_weight = 2*n + 1
        boltzmann_factor = exp(-(n^2)/T)
        return statistical_weight * boltzmann_factor
    end
    
    return result
end

# 示例：混合积分-求和
function mixed_integral_sum_example()
    """
    展示混合积分-求和的用法
    
    计算：∫₀²⁰ dp Σₙ₌₀¹⁰ p² exp(-√(p² + n²))
    """
    p_grid = create_momentum_grid(64, 20.0)
    n_indices = 0:10
    method = GaussLegendreIntegration()
    
    integrand_sum = function(p, n)
        energy = sqrt(p^2 + n^2)
        return p^2 * exp(-energy)
    end
    
    result = mixed_integral_sum(method, p_grid, n_indices, integrand_sum)
    
    return result
end

"""
总结：求和接口的使用场景

1. discrete_sum(): 
   - 纯离散求和
   - 用于角动量态统计、能级求和等

2. mixed_integral_sum():
   - 混合积分-求和计算
   - 用于连续积分变量 + 离散求和变量的情况

3. angular_momentum_sum():
   - 旋转系统的专用函数
   - 自动包含角动量统计权重 (2n+1)
   - 是 mixed_integral_sum() 的特化版本

这些函数完善了积分接口，使其能够处理：
- 连续积分（已有）
- 离散求和（新增）
- 混合积分-求和（新增）

从而覆盖了PNJL物理模型中的所有数值计算需求。
"""
