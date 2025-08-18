"""
Common physical constants and unit conversions.

This module provides fundamental physical constants and conversion factors
used throughout the PNJL physics simulation package.
"""
module PhysicalConstants

export π, hc, MeV_to_fm_inv

# Mathematical constants
const π = 3.141592653589793

# Physical constants
const hc = 197.33  # MeV⋅fm (ℏc conversion factor)

# Unit conversion functions
"""
    MeV_to_fm_inv(energy_MeV)

Convert energy from MeV to fm⁻¹ units using ℏc = 197.33 MeV⋅fm.
"""
MeV_to_fm_inv(energy_MeV) = energy_MeV / hc

end  # module PhysicalConstants
