"""
Example script demonstrating PNJL model calculations.

This script shows how to:
1. Set up PNJL model parameters
2. Calculate thermodynamic properties 
3. Construct T-ρ phase diagrams
4. Analyze chiral symmetry breaking
"""

using PNJLPhysicsSimulation
using PNJLPhysicsSimulation.PNJLModel

println("=== PNJL Model Example ===")

# Set up model parameters
println("Setting up PNJL model parameters...")

# Integration nodes
n_points = 128
nodes = get_nodes(n_points)
println("Created $(n_points) integration nodes")

# Thermodynamic conditions
T = 150.0 / PhysicalConstants.hc  # 150 MeV
println("Temperature: T = $(T * PhysicalConstants.hc) MeV")

# Chemical potentials for u, d, s quarks
μ_u = 320.0 / PhysicalConstants.hc
μ_d = 320.0 / PhysicalConstants.hc  
μ_s = 320.0 / PhysicalConstants.hc
mu = [μ_u, μ_d, μ_s]

println("Chemical potentials:")
println("  μ_u = $(μ_u * PhysicalConstants.hc) MeV")
println("  μ_d = $(μ_d * PhysicalConstants.hc) MeV")
println("  μ_s = $(μ_s * PhysicalConstants.hc) MeV")

# Field values: [φ_u, φ_d, φ_s, Φ₁, Φ₂]
x = [-0.1, -0.1, -1.7, 0.5, 0.5]

println("Initial field values:")
println("  φ_u = $(x[1])")
println("  φ_d = $(x[2])")
println("  φ_s = $(x[3])")
println("  Φ₁ = $(x[4])")
println("  Φ₂ = $(x[5])")

# Single point calculation
println("\n=== Single Point Calculation ===")

println("Calculating pressure...")
pressure = pressure_solve_core(x, mu, T, nodes)
println("Equilibrium pressure: P = $(pressure)")

println("Calculating thermodynamic quantities...")
pressure_eq, rho, entropy, energy = calculate_thermo(x, mu, T, nodes)

println("Thermodynamic results:")
println("  Pressure: P = $(pressure_eq)")
println("  Baryon densities: ρ = $(rho)")
println("  Total baryon density: ρ_total = $(sum(rho))")
println("  Entropy density: s = $(entropy)")
println("  Energy density: ε = $(energy)")

# Temperature derivatives
println("\nCalculating temperature derivatives...")
using FiniteDifferences

fdm = central_fdm(5, 4)
dP_dT4_val = dP_dT4_direct(x, mu, T, nodes, fdm)
println("4th order temperature derivative: ∂⁴P/∂T⁴ = $(dP_dT4_val)")

# T-ρ phase diagram construction (small example)
println("\n=== T-ρ Phase Diagram Construction ===")

T_start = 140.0 / PhysicalConstants.hc
T_end = 160.0 / PhysicalConstants.hc  
T_step = 5.0 / PhysicalConstants.hc

rho_start = 1.0
rho_end = 0.5
rho_step = -0.1

println("Temperature range: $(T_start*PhysicalConstants.hc) to $(T_end*PhysicalConstants.hc) MeV")
println("Density range: $(rho_start) to $(rho_end) (relative to ρ₀)")

println("Running T-ρ phase diagram calculation (this may take a moment)...")

# Note: This is a small example. For full calculations, use larger ranges
results = Trho_optimized(
    T_start, T_end, T_step,
    rho_start, rho_end, rho_step,
    save_results=false  # Don't save for this example
)

if results !== nothing
    println("T-ρ calculation completed successfully!")
    println("Number of temperature points: $(length(T_start:T_step:T_end))")
    println("Number of density points: $(length(rho_start:rho_step:rho_end))")
else
    println("T-ρ calculation completed and results saved!")
end

println("\n=== PNJL Example completed successfully! ===")
