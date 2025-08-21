module AutodiffInterface

using ForwardDiff

export compute_gradient,
       compute_temperature_derivative,
       compute_chemical_potential_derivatives

"""
通用自动微分接口（简化版本，重命名为 AutodiffInterface）

保留函数：
- `compute_gradient`
- `compute_temperature_derivative`
- `compute_chemical_potential_derivatives`
"""

function compute_gradient(thermodynamic_potential, variables::Vector{T}, params...) where T<:Real
    try
        gradient = ForwardDiff.gradient(x -> thermodynamic_potential(x, params...), variables)
        return gradient
    catch e
        error("自动微分梯度计算失败: $e\n请检查热力学势函数是否支持ForwardDiff，以及变量维度是否正确")
    end
end

function compute_temperature_derivative(thermodynamic_potential, variables::Vector{T},
                                        temperature::Float64, other_params...) where T<:Real
    try
        temp_func = T_val -> thermodynamic_potential(variables, T_val, other_params...)
        dOmega_dT = ForwardDiff.derivative(temp_func, temperature)
        return dOmega_dT
    catch e
        error("温度导数计算失败: $e\n请检查热力学势函数是否正确支持温度参数")
    end
end

function compute_chemical_potential_derivatives(thermodynamic_potential, variables::Vector{T},
                                                chemical_potentials::Vector{Float64}, other_params...) where T<:Real
    try
        mu_func = mu_vec -> thermodynamic_potential(variables, mu_vec, other_params...)
        dOmega_dmu = ForwardDiff.gradient(mu_func, chemical_potentials)
        return dOmega_dmu
    catch e
        error("化学势导数计算失败: $e\n请检查热力学势函数的化学势参数接口")
    end
end

end # module AutodiffInterface

# 兼容绑定
const compute_gradient = AutodiffInterface.compute_gradient
const compute_temperature_derivative = AutodiffInterface.compute_temperature_derivative
const compute_chemical_potential_derivatives = AutodiffInterface.compute_chemical_potential_derivatives
