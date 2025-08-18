"""
Universal Physical Constants

This module provides fundamental physical constants used across all physics models.
Only includes truly universal constants - model-specific parameters are in their respective modules.
"""
module PhysicalConstants

export hc, Nc, MeV_to_fm_inv

# Universal physical constants
const hc = 197.33  # ℏc (MeV⋅fm) - conversion factor between natural and SI units
const Nc = 3       # Number of colors in QCD (universal gauge theory constant)

# Unit conversion functions
"""
    MeV_to_fm_inv(energy_MeV)

Convert energy from MeV to fm⁻¹ units using ℏc = 197.33 MeV⋅fm.
"""
MeV_to_fm_inv(energy_MeV) = energy_MeV / hc

end  # module PhysicalConstants
