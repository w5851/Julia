"""
Safe mathematical functions module.

This module provides numerically stable versions of common mathematical
functions to avoid overflow, underflow, and domain errors.
"""
module MathUtils

using SpecialFunctions: log, exp, sqrt

export safe_log, safe_exp, safe_sqrt, safe_division

"""
    safe_log(x; min_val=1e-16, handle_negative=:clamp)

Safe logarithm function that handles negative and near-zero values.

# Arguments
- `x`: Input value
- `min_val`: Minimum value to clamp to (default: 1e-16)
- `handle_negative`: How to handle negative values
  - `:clamp`: Clamp to min_val
  - `:error`: Throw error
  - `:nan`: Return NaN

# Returns
- Safe logarithm value

# Examples
```julia
result = safe_log(-1.0)  # Returns log(1e-16)
result = safe_log(1e-20) # Returns log(1e-16)
result = safe_log(2.0)   # Returns log(2.0)
```

# Notes
- Ensures numerical stability by preventing log of negative or zero values
- Default behavior clamps values to minimum positive value
"""
@inline function safe_log(x; min_val=1e-16, handle_negative=:clamp)
    if x <= 0
        if handle_negative == :error
            error("safe_log: Input must be positive, got $x")
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

Safe exponential function that prevents overflow.

# Arguments
- `x`: Input value
- `max_val`: Maximum value to prevent overflow (default: 700.0)

# Returns
- Safe exponential value

# Examples
```julia
result = safe_exp(800.0)  # Returns exp(700.0)
result = safe_exp(1.0)    # Returns exp(1.0)
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

Safe square root function that handles negative values.

# Arguments
- `x`: Input value
- `min_val`: Minimum value to clamp to (default: 0.0)

# Returns
- Safe square root value
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

Safe division function that prevents division by zero.

# Arguments
- `numerator`: Numerator value
- `denominator`: Denominator value
- `epsilon`: Small value to prevent division by zero (default: 1e-16)

# Returns
- Safe division result
"""
@inline function safe_division(numerator, denominator; epsilon=1e-16)
    if abs(denominator) < epsilon
        return sign(numerator) * sign(denominator) * abs(numerator) / epsilon
    else
        return numerator / denominator
    end
end

end # module MathUtils
