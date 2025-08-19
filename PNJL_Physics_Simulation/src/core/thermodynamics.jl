"""
通用热力学和统计分布函数。

此模块为不同物理模型提供共享的热力学函数。
"""
module Thermodynamics

using SpecialFunctions: exp, log

export fermion_distribution, fermion_anti_distribution, 
       bosonic_distribution, calculate_log_term

"""
    fermion_distribution(E, μ, T)

费米-狄拉克分布函数。

# 参数
- `E`: 能量
- `μ`: 化学势
- `T`: 温度

# 返回值
- 费米-狄拉克分布值
"""
@inline function fermion_distribution(E, μ, T)
    return 1 / (exp((E - μ) / T) + 1)
end

"""
    fermion_anti_distribution(E, μ, T)

反费米子分布函数。
"""
@inline function fermion_anti_distribution(E, μ, T)
    return 1 / (exp((E + μ) / T) + 1)
end

"""
    bosonic_distribution(E, μ, T)

玻色-爱因斯坦分布函数。
"""
@inline function bosonic_distribution(E, μ, T)
    return 1 / (exp((E - μ) / T) - 1)
end

"""
    calculate_log_term(E, μ, T)

计算热力学势的对数项。
"""
@inline function calculate_log_term(E, μ, T)
    x = E - μ
    x_anti = E + μ
    term1 = 1 + exp(-x / T)
    term2 = 1 + exp(-x_anti / T)
    return log(term1) + log(term2)
end

end  # module Thermodynamics
