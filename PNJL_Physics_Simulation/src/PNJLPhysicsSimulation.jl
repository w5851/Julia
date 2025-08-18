"""
PNJL Physics Simulation Package

A comprehensive Julia package for studying QCD phase transitions.
"""
module PNJLPhysicsSimulation

# Core modules
include("core/constants.jl")
include("core/math_utils.jl")
include("core/integration.jl")  
include("core/thermodynamics.jl")

# Model modules
include("models/gas_liquid/constants.jl")
include("models/gas_liquid/functions.jl")
include("models/pnjl/constants.jl")
include("models/pnjl/functions.jl")
include("models/pnjl_aniso/constants.jl")
include("models/pnjl_aniso/functions.jl")
include("models/rotation/constants.jl")
include("models/rotation/functions.jl")

# Export everything for now to test
export PhysicalConstants, Integration, Thermodynamics, MathUtils,
       GasLiquidConstants, GasLiquidFunctions,
       PNJLConstants, PNJLFunctions, 
       PNJLAnisoConstants, RotationConstants

end  # module PNJLPhysicsSimulation
