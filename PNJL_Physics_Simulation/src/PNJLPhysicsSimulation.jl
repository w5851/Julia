"""
PNJL Physics Simulation Package

A comprehensive Julia package for studying QCD phase transitions.
"""
module PNJLPhysicsSimulation

# Core modules
include("core/constants.jl")
include("core/math_utils.jl")
include("core/integration.jl")
include("core/integration_interface.jl")
include("core/thermodynamics.jl")
include("core/function_registry.jl")  # Function registry

# Model modules - must be included before the unified interface
include("models/gas_liquid/constants.jl")
include("models/gas_liquid/functions.jl")
include("models/pnjl/constants.jl")
include("models/pnjl/functions.jl")
include("models/pnjl_aniso/constants.jl")
include("models/pnjl_aniso/functions.jl")
include("models/rotation/constants.jl")
include("models/rotation/functions.jl")

# High-level interfaces - include after models
include("core/model_configuration.jl")
# include("core/unified_physics_interface.jl")  # 暂时注释掉避免依赖问题

# Export everything for now to test
export PhysicalConstants, Integration, IntegrationInterface, ModelConfiguration, 
       FunctionRegistry, Thermodynamics, MathUtils,
       GasLiquidConstants, GasLiquidFunctions,
       PNJLConstants, PNJLFunctions, 
       PNJLAnisoConstants, PNJLAnisoFunctions,
       RotationConstants, RotationFunctions

end  # module PNJLPhysicsSimulation
