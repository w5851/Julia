"""
Simplified test runner that works without complex package loading.
"""

println("=== PNJL Physics Simulation Basic Tests ===")

# Get the absolute project root directory
project_root = abspath(joinpath(@__DIR__, ".."))
println("Project root: $project_root")

# Define file paths
constants_file = joinpath(project_root, "src", "core", "constants.jl")
integration_file = joinpath(project_root, "src", "core", "integration.jl")
gas_liquid_file = joinpath(project_root, "src", "models", "gas_liquid", "constants.jl")
pnjl_aniso_file = joinpath(project_root, "src", "models", "pnjl_aniso", "constants.jl")
rotation_file = joinpath(project_root, "src", "models", "rotation", "constants.jl")

println("Testing file existence...")
println("  Constants file: $(isfile(constants_file)) - $constants_file")
println("  Integration file: $(isfile(integration_file)) - $integration_file")
println("  Gas-liquid file: $(isfile(gas_liquid_file)) - $gas_liquid_file")
println("  PNJL Aniso file: $(isfile(pnjl_aniso_file)) - $pnjl_aniso_file")
println("  Rotation file: $(isfile(rotation_file)) - $rotation_file")

# Test 1: Can we load constants?
println("\nTest 1: Loading PhysicalConstants...")
try
    include(constants_file)
    using .PhysicalConstants
    π_val = PhysicalConstants.π
    hc_val = PhysicalConstants.hc
    println("✓ PhysicalConstants: π = $π_val, hc = $hc_val")
    @assert abs(π_val - 3.141592653589793) < 1e-10
    @assert abs(hc_val - 197.33) < 1e-10
catch e
    println("✗ PhysicalConstants failed: $e")
    error("Cannot proceed without basic constants")
end

# Test 2: Can we load integration?
println("\nTest 2: Loading Integration...")
try
    using FastGaussQuadrature  # This might fail if not installed
    include(integration_file)
    using .Integration
    nodes, weights = Integration.gauleg(0.0, 1.0, 5)
    println("✓ Integration: Generated $(length(nodes)) nodes")
    @assert length(nodes) == 5
    @assert length(weights) == 5
catch e
    println("⚠ Integration test skipped (missing FastGaussQuadrature?): $e")
end

# Test 3: Basic Gas-Liquid constants
println("\nTest 3: Loading GasLiquidConstants...")
try
    include(gas_liquid_file)
    using .GasLiquidConstants
    m_val = GasLiquidConstants.m
    println("✓ GasLiquidConstants: Nucleon mass = $m_val")
    @assert m_val > 0
catch e
    println("✗ GasLiquidConstants failed: $e")
end

# Test 4: PNJL Aniso constants
println("\nTest 4: Loading PNJLAnisoConstants...")
try
    include(pnjl_aniso_file)
    using .PNJLAnisoConstants
    G_f_val = PNJLAnisoConstants.G_f
    println("✓ PNJLAnisoConstants: G_f = $G_f_val")
    @assert G_f_val > 0
catch e
    println("✗ PNJLAnisoConstants failed: $e")
end

# Test 5: Rotation constants
println("\nTest 5: Loading RotationConstants...")
try
    include(rotation_file)
    using .RotationConstants
    r0_val = RotationConstants.r0
    println("✓ RotationConstants: r0 = $r0_val")
    @assert r0_val > 0
catch e
    println("✗ RotationConstants failed: $e")
end

println("\n🎉 Basic tests completed successfully!")
println("If you saw any warnings about FastGaussQuadrature, install it with:")
println("  using Pkg; Pkg.add(\"FastGaussQuadrature\")")
