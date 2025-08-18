"""
Constants and parameters for the anisotropic PNJL model.

This module contains PNJL_aniso model-specific constants and parameters.
Universal physical constants are imported from PhysicalConstants.
"""
module PNJLAnisoConstants

using ...PhysicalConstants: hc, Nc, MeV_to_fm_inv

export rho0, T0, Lambda_f, G_f, K_f, m0, m0_q_f, m0_s_f
export a0, a1, a2, b3, b4

# Nuclear matter properties for PNJL_aniso model
const rho0 = 0.16  # Nuclear saturation density (fm⁻³)
const T0 = MeV_to_fm_inv(210)  # Critical temperature (fm⁻¹)

# PNJL_aniso model parameters (same as PNJL but may have different parameter sets)
const Lambda = 602.3      # Cutoff parameter (MeV)
const G_Lam2 = 1.835      # Four-fermion coupling strength
const K_Lam5 = 12.36      # Six-fermion coupling strength

const m0_q = 5.5    # Up/down quark mass (MeV)
const m0_s = 140.7  # Strange quark mass (MeV)

# Convert to natural units (fm⁻¹)
const Lambda_f = MeV_to_fm_inv(Lambda)
const G_f = G_Lam2 / Lambda_f^2
const K_f = K_Lam5 / Lambda_f^5

const m0_q_f = MeV_to_fm_inv(m0_q)
const m0_s_f = MeV_to_fm_inv(m0_s)

# Quark mass array [up, down, strange]
const m0 = [m0_q_f, m0_q_f, m0_s_f]

# Polyakov loop potential parameters
const a0 = 3.51
const a1 = -2.47
const a2 = 15.2
const b3 = -1.75
const b4 = 7.555

end  # module PNJLAnisoConstants
