#!/usr/bin/env julia

println("=== Testing PNJL Models ===")

try
    # Test basic module loading
    include("src/core/constants.jl")
    include("src/core/integration.jl")
    include("src/models/pnjl_aniso/constants.jl")
    include("src/models/rotation/constants.jl")
    
    println("✓ All modules loaded successfully")
    
    # Test constants
    println("PNJL Aniso hc: ", PNJLAnisoConstants.hc)
    println("Rotation hc: ", RotationConstants.hc)
    println("✓ Constants accessible")
    
    # Test integration
    nodes, weights = Integration.gauleg(0.0, 1.0, 5)
    println("✓ Integration function works")
    
    println("SUCCESS: All basic components working!")
    
catch e
    println("ERROR: ", e)
    println("Please check file paths and dependencies")
end
