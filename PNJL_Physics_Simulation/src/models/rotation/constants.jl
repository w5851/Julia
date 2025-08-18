"""
Constants for the Rotation Model

Physical constants and model parameters for the PNJL model with rotation effects.
Based on the original Constants_Rotation.jl file.
"""
module RotationConstants

# Physical constants
const π = 3.141592653589793
const hc = 197.33

# Nuclear matter
const rho0 = 0.16
const T0 = 270 / hc

# PNJL model constants
const Nc = 3.0  # Number of colors

const Lambda = 650.0    # (MeV)
const G_Lam2 = 4.93 / 1e6  # (MeV), self-coupling
const K_Lam5 = 12.36      # (MeV), three-flavor coupling

const m0_q = 5.0       # (MeV)
const m0_s = 5.0       # (MeV)
const m0_q_f = m0_q / hc   # (fm⁻¹)
const m0_s_f = m0_s / hc   # (fm⁻¹)

const Lambda_f = Lambda / hc      # (fm⁻¹)
const G_f = G_Lam2 * hc^2   # (fm²)
const K_f = K_Lam5 / Lambda_f^5  # (fm⁵)

# Polyakov loop parameters
const a0 = 6.75
const a1 = -1.95
const a2 = 2.625
const a3 = -7.44
const b3 = 0.75
const b4 = 7.5

# Rotation-specific parameters
const r0 = 1.2  # Nuclear radius parameter (fm)
const C = 0.5   # Coupling parameter for rotation effects

# Rotation-dependent coefficients for Polyakov potential
const coefficients = Dict{String, Vector{Float64}}(
    "a" => [1.0, 0.1, 0.01],    # a₀, a₁ω², a₂ω⁴
    "b" => [1.0, 0.05, 0.001],  # b₀, b₁ω², b₂ω⁴  
    "c" => [0.5, 0.02, 0.0001], # c₀, c₁ω², c₂ω⁴
    "d" => [0.0, 0.1, 0.01]     # d₀, d₁ω², d₂ω⁴
)

end  # module RotationConstants
