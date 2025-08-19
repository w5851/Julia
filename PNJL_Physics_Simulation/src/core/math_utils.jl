"""
安全数学函数模块。

此模块提供常用数学函数的数值稳定版本，
用于避免溢出、下溢和定义域错误。
"""
module MathUtils

using SpecialFunctions: log, exp, sqrt

export safe_log, safe_exp, safe_sqrt, safe_division

"""
    safe_log(x; min_val=1e-16, handle_negative=:clamp)

安全对数函数，处理负值和接近零的值。

# 参数
- `x`: 输入值
- `min_val`: 要限制到的最小值 (默认: 1e-16)
- `handle_negative`: 如何处理负值
  - `:clamp`: 限制到 min_val
  - `:error`: 抛出错误
  - `:nan`: 返回 NaN

# 返回值
- 安全的对数值

# 示例
```julia
result = safe_log(-1.0)  # 返回 log(1e-16)
result = safe_log(1e-20) # 返回 log(1e-16)
result = safe_log(2.0)   # 返回 log(2.0)
```

# 注意事项
- 通过防止负值或零值的对数来确保数值稳定性
- 默认行为将值限制到最小正值
"""
@inline function safe_log(x; min_val=1e-16, handle_negative=:clamp)
    if x <= 0
        if handle_negative == :error
            error("safe_log: 输入必须为正数, 得到 $x")
        elseif handle_negative == :nan
            return NaN
        else  # :clamp
            return log(min_val)
        end
    elseif x < min_val
        return log(min_val)
    else
        return log(x)
    end
end

"""
    safe_exp(x; max_val=700.0)

安全指数函数，防止溢出。

# 参数
- `x`: 输入值
- `max_val`: 防止溢出的最大值 (默认: 700.0)

# 返回值
- 安全的指数值

# 示例
```julia
result = safe_exp(800.0)  # 返回 exp(700.0)
result = safe_exp(1.0)    # 返回 exp(1.0)
```
"""
@inline function safe_exp(x; max_val=700.0)
    if x > max_val
        return exp(max_val)
    else
        return exp(x)
    end
end

"""
    safe_sqrt(x; min_val=0.0)

安全平方根函数，处理负值。

# 参数
- `x`: 输入值
- `min_val`: 要限制到的最小值 (默认: 0.0)

# 返回值
- 安全的平方根值
"""
@inline function safe_sqrt(x; min_val=0.0)
    if x < min_val
        return sqrt(min_val)
    else
        return sqrt(x)
    end
end

"""
    safe_division(numerator, denominator; epsilon=1e-16)

安全除法函数，防止除零错误。

# 参数
- `numerator`: 分子值
- `denominator`: 分母值
- `epsilon`: 防止除零的小值 (默认: 1e-16)

# 返回值
- 安全的除法结果
"""
@inline function safe_division(numerator, denominator; epsilon=1e-16)
    if abs(denominator) < epsilon
        return sign(numerator) * sign(denominator) * abs(numerator) / epsilon
    else
        return numerator / denominator
    end
end

end # module MathUtils
