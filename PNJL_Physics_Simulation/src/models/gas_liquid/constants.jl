"""
Constants and parameters for the Gas-Liquid phase transition model.

This module contains only Gas-Liquid model-specific constants and parameters.
Universal physical constants are imported from PhysicalConstants.
"""
module GasLiquidConstants

using ...PhysicalConstants: hc, MeV_to_fm_inv

export m, mσ, mω, mρ, mδ, rho0_default, calculate_couplings

# Particle masses (MeV converted to fm⁻¹)
const m = MeV_to_fm_inv(939.0)    # Nucleon mass
const mσ = MeV_to_fm_inv(550.0)   # Sigma meson mass
const mω = MeV_to_fm_inv(783.0)   # Omega meson mass  
const mρ = MeV_to_fm_inv(775.0)   # Rho meson mass
const mδ = MeV_to_fm_inv(980.0)   # Delta meson mass

# Nuclear matter properties specific to Gas-Liquid model
const rho0_default = 0.16  # Nuclear saturation density (fm⁻³)

"""
    calculate_couplings(ρ0, B_A, K, m_ratio, E_sym)

Calculate coupling constants for the Gas-Liquid model based on nuclear matter properties.

# Arguments
- `ρ0`: Nuclear density (fm⁻³)
- `B_A`: Binding energy per nucleon (MeV, already converted by /hc)
- `K`: Incompressibility (MeV, already converted by /hc)  
- `m_ratio`: Effective mass ratio (meff/m)
- `E_sym`: Symmetry energy (MeV, already converted by /hc)

# Returns
- `(fσ, fω, fρ, fδ, b, c)`: Tuple of coupling constants and parameters
"""
function calculate_couplings(ρ0, B_A, K, m_ratio, E_sym)
    # Calculate effective mass
    meff = m_ratio * m
    fδ = 0.0  # Delta coupling fixed at 0.0
        
    # Fermi momentum and energy
    kF = (1.5π^2 * ρ0)^(1/3)
    EF = sqrt(kF^2 + meff^2)
    
    # Sigma field strength
    gσ = m - meff
    
    # Omega coupling constant
    fω = (m + B_A - EF) / ρ0
    
    # Intermediate variables
    x = kF / meff
    t = sqrt(1 + x^2)

    term1 = 0.5x * t + x / t - 1.5log(x + t)
    I1 = (2 / π^2) * meff^2 * term1
    
    term2 = 0.25(x * t^3 - 0.5x * t - 0.5log(x + t))
    I2 = (2 / π^2) * meff^4 * term2
    
    term3 = 0.5(x * t - log(x + t))
    I3 = (2 / π^2) * meff^3 * term3
    
    # Calculate α, β, γ, δ coefficients
    α1 = K - fω * (6kF^3 / π^2) - 3kF^2 / EF
    β1 = 2m * gσ * α1
    γ1 = 3gσ^2 * α1
    δ1 = -(6kF^3 / π^2) * (meff / EF)^2 - α1 * I1
    
    α2 = 0.5gσ^2
    β2 = (1/3)m * gσ^3
    γ2 = 0.25gσ^4
    δ2 = ρ0 * (m + B_A) - 0.5fω * ρ0^2 - I2
    
    α3 = gσ
    β3 = m * gσ^2
    γ3 = gσ^3
    δ3 = I3
    
    # Calculate parameters c and b
    denom1 = (γ1 * α2 - γ2 * α1) * (β2 * α3 - β3 * α2) -
             (γ2 * α3 - γ3 * α2) * (β1 * α2 - β2 * α1)
             
    num1 = (δ1 * α2 - δ2 * α1) * (β2 * α3 - β3 * α2) -
           (δ2 * α3 - δ3 * α2) * (β1 * α2 - β2 * α1)
           
    c = num1 / denom1
    
    denom2 = β1 * α2 - β2 * α1
    num2 = (δ1 * α2 - δ2 * α1) - (γ1 * α2 - γ2 * α1) * c
    b = num2 / denom2
    
    # Sigma coupling constant
    fσ = α1 / (δ1 - β1 * b - γ1 * c)
    
    # Rho coupling constant
    term_fρ1 = (meff / EF)^2
    term_fρ2 = 1 + 0.25fδ * I1
    fρ = (2 / ρ0) * (E_sym - (1/6)*(kF^2 / EF) + (ρ0 / 8) * fδ * term_fρ1 / term_fρ2)
    
    return (fσ, fω, fρ, fδ, b, c)
end

end  # module GasLiquidConstants
