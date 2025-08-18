"""
Constants for the Rotation Model

This module contains rotation model-specific constants and parameters.
Universal physical constants are imported from PhysicalConstants.
Based on the original Constants_Rotation.jl file.
"""
module RotationConstants

using ...PhysicalConstants: hc, Nc, MeV_to_fm_inv

export rho0, T0, Lambda_f, G_f, K_f, m0_q_f, m0_s_f, r0, C, coefficients
export a0, a1, a2, a3, b3, b4

# Nuclear matter properties for Rotation model
const rho0 = 0.16               # Nuclear saturation density (fm⁻³)
const T0 = MeV_to_fm_inv(270)   # Critical temperature (fm⁻¹)

# PNJL model constants (Rotation-specific parameter set)
const Lambda = 650.0             # Cutoff (MeV)
const G_Lam2 = 4.93 / 1e6       # Self-coupling (MeV⁻¹)
const K_Lam5 = 12.36            # Three-flavor coupling

const m0_q = 5.0    # Up/down quark mass (MeV)
const m0_s = 5.0    # Strange quark mass (MeV)

# Convert to natural units (fm⁻¹)
const m0_q_f = MeV_to_fm_inv(m0_q)
const m0_s_f = MeV_to_fm_inv(m0_s)

const Lambda_f = MeV_to_fm_inv(Lambda)
const G_f = G_Lam2 * hc^2   # (fm²)
const K_f = K_Lam5 / Lambda_f^5  # (fm⁵)

# Polyakov loop parameters (Rotation-specific values)
const a0 = 6.75
const a1 = -1.95
const a2 = 2.625
const a3 = -7.44
const b3 = 0.75
const b4 = 7.5

# Rotation-specific parameters
const r0 = 0.1 / 1000 * hc
const C = 4.0

# Rotation-dependent coefficients for Polyakov potential
const coefficients = Dict{String, Vector{Float64}}(
    "a" => [0.0454431, -1.27942e-5, -5.43339e-9],
    "b" => [46.8263, -0.0210165, -2.15394e-5],
    "c" => [1.00298, 1.55157e-4, -5.99032e-8],
    "d" => [0.0600157, -5.74388e-6, -8.24192]
)

end  # module RotationConstants
