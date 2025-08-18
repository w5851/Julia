"""
Common thermodynamic and statistical distribution functions.

This module provides shared thermodynamic functions used across
different physics models.
"""
module Thermodynamics

using SpecialFunctions: exp, log

export fermion_distribution, fermion_anti_distribution, 
       bosonic_distribution, calculate_log_term

"""
    fermion_distribution(E, μ, T)

Fermi-Dirac distribution function.

# Arguments
- `E`: Energy
- `μ`: Chemical potential
- `T`: Temperature

# Returns
- Fermi-Dirac distribution value
"""
@inline function fermion_distribution(E, μ, T)
    return 1 / (exp((E - μ) / T) + 1)
end

"""
    fermion_anti_distribution(E, μ, T)

Anti-fermion distribution function.
"""
@inline function fermion_anti_distribution(E, μ, T)
    return 1 / (exp((E + μ) / T) + 1)
end

"""
    bosonic_distribution(E, μ, T)

Bose-Einstein distribution function.
"""
@inline function bosonic_distribution(E, μ, T)
    return 1 / (exp((E - μ) / T) - 1)
end

"""
    calculate_log_term(E, μ, T)

Calculate logarithmic terms for thermodynamic potentials.
"""
@inline function calculate_log_term(E, μ, T)
    x = E - μ
    x_anti = E + μ
    term1 = 1 + exp(-x / T)
    term2 = 1 + exp(-x_anti / T)
    return log(term1) + log(term2)
end

end  # module Thermodynamics
