"""
Package installation and setup script.

This script installs all required dependencies and sets up the environment
for the PNJL Physics Simulation package.
"""

using Pkg

println("=== PNJL Physics Simulation Package Setup ===")
println("Installing required packages...")

# Required packages
packages = [
    "FastGaussQuadrature",
    "ForwardDiff", 
    "NLsolve",
    "FiniteDifferences",
    "BenchmarkTools",
    "StaticArrays",
    "SpecialFunctions"
]

# Install packages
for pkg in packages
    try
        println("Installing $pkg...")
        Pkg.add(pkg)
        println("  ✓ $pkg installed successfully")
    catch e
        println("  ✗ Failed to install $pkg: $e")
    end
end

# Activate the local project
println("\nActivating project environment...")
try
    Pkg.activate(".")
    Pkg.instantiate()
    println("  ✓ Project environment activated")
catch e
    println("  ✗ Failed to activate project: $e")
end

# Run tests to verify installation
println("\nRunning tests to verify installation...")
try
    include("test/runtests.jl")
    println("  ✓ All tests passed!")
catch e
    println("  ✗ Tests failed: $e")
    println("  Please check the installation and try again.")
end

println("\n=== Setup Complete ===")
println("You can now use the package with:")
println("  using PNJLPhysicsSimulation")
println("  using PNJLPhysicsSimulation.GasLiquidModel")
println("  using PNJLPhysicsSimulation.PNJLModel")
println("\nSee scripts/ directory for usage examples.")
println("See docs/ directory for detailed documentation.")
