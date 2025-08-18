"""
Unified Constants Management

This module provides centralized constant definitions to avoid redefinition warnings
and ensure consistency across all physics models.
"""

module UnifiedConstants

# Export all physical constants
export π, hc, rho0, T0, Nc
export Lambda, G_Lam2, K_Lam5
export m0_q, m0_s, m0_q_f, m0_s_f
export Lambda_f, G_f, K_f
export get_physical_constants, get_model_constants

# Universal physical constants (avoid redefinition of Base.π)
import Base: π as base_pi
const physics_pi = Float64(base_pi)  # Use Julia's built-in π

const hc = 197.33  # ħc (MeV⋅fm)

# Nuclear constants
const rho0 = 0.16  # Nuclear saturation density (fm⁻³)
const T0 = 270 / hc  # Reference temperature (fm⁻¹)

# Color dynamics
const Nc = 3.0  # Number of colors

# Model-specific constants (PNJL)
const Lambda = 650.0    # Cutoff (MeV)
const G_Lam2 = 4.93 / 1e6  # Self-coupling (MeV⁻¹)
const K_Lam5 = 12.36      # Three-flavor coupling (dimensionless)

# Quark masses
const m0_q = 5.0       # Light quark current mass (MeV)
const m0_s = 5.0       # Strange quark current mass (MeV)

# Converted units (natural units, fm⁻¹)
const m0_q_f = m0_q / hc
const m0_s_f = m0_s / hc
const Lambda_f = Lambda / hc
const G_f = G_Lam2 * hc^2
const K_f = K_Lam5 / Lambda_f^5

"""
    get_physical_constants() -> Dict{Symbol, Float64}

Return a dictionary of all physical constants.
"""
function get_physical_constants()
    return Dict(
        :pi => physics_pi,
        :hc => hc,
        :rho0 => rho0,
        :T0 => T0,
        :Nc => Nc
    )
end

"""
    get_model_constants(model_type::Symbol) -> Dict{Symbol, Float64}

Return model-specific constants.

# Arguments
- `model_type`: :pnjl, :pnjl_aniso, :rotation, or :gas_liquid
"""
function get_model_constants(model_type::Symbol)
    base_constants = get_physical_constants()
    
    if model_type in [:pnjl, :pnjl_aniso, :rotation]
        merge!(base_constants, Dict(
            :Lambda => Lambda,
            :G_Lam2 => G_Lam2,
            :K_Lam5 => K_Lam5,
            :m0_q => m0_q,
            :m0_s => m0_s,
            :m0_q_f => m0_q_f,
            :m0_s_f => m0_s_f,
            :Lambda_f => Lambda_f,
            :G_f => G_f,
            :K_f => K_f
        ))
    elseif model_type == :gas_liquid
        # Gas-liquid model may have different parameters
        # Add specific constants if needed
    end
    
    return base_constants
end

"""
    @with_constants(model_type, block)

Macro to make constants available in a local scope without importing conflicts.
"""
macro with_constants(model_type, block)
    quote
        local constants = get_model_constants($(esc(model_type)))
        local π = constants[:pi]
        local hc = constants[:hc]
        local rho0 = constants[:rho0]
        local T0 = constants[:T0]
        local Nc = constants[:Nc]
        
        $(esc(block))
    end
end

end # module UnifiedConstants
