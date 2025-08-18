"""
Constants and parameters for the PNJL (Polyakov-Nambu-Jona-Lasinio) model.

This module contains all physical constants and model parameters
specific to the PNJL model for QCD phase transitions.
"""
module PNJLConstants

# 直接定义常数
const π = 3.141592653589793
const hc = 197.33
MeV_to_fm_inv(energy_MeV) = energy_MeV / hc

# Nuclear matter properties
const rho0 = 0.16  # Nuclear saturation density (fm⁻³)

# Polyakov loop parameters
const T0 = 210 / hc  # Critical temperature (fm⁻¹)

# Color degrees of freedom
const Nc = 3.0

# Model parameters (in MeV, then converted to natural units)
const Lambda = 602.3      # Cutoff parameter (MeV)
const G_Lam2 = 1.835      # Four-fermion coupling strength
const K_Lam5 = 12.36      # Six-fermion coupling strength

# Current quark masses (MeV)
const m0_q = 5.5    # Up/down quark mass
const m0_s = 140.7  # Strange quark mass

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
const b4 = 7.5

export rho0, T0, Nc, Lambda_f, G_f, K_f, m0, m0_q_f, m0_s_f,
       a0, a1, a2, b3, b4

end  # module PNJLConstants
