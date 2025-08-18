"""
Example script demonstrating Gas-Liquid model calculations.

This script shows how to:
1. Set up model parameters and coupling constants
2. Calculate pressure and thermodynamic derivatives  
3. Analyze phase transition properties
4. Save results for further analysis
"""

using PNJLPhysicsSimulation
using PNJLPhysicsSimulation.GasLiquidModel

println("=== Gas-Liquid Model Example ===")

# Set up physical parameters
println("Setting up model parameters...")

# Integration nodes
n_points = 256
nodes = get_nodes(n_points)
println("Created $(n_points) integration nodes")

# Temperature and thermodynamic conditions
T = 50.0 / PhysicalConstants.hc  # 50 MeV in natural units
println("Temperature: T = $(T * PhysicalConstants.hc) MeV")

# Nuclear matter properties for coupling calculation
ρ₀ = 0.16                    # Nuclear saturation density (fm⁻³)
B_A = -16.0 / PhysicalConstants.hc  # Binding energy per nucleon 
K = 240.0 / PhysicalConstants.hc    # Incompressibility
m_ratio = 0.75               # Effective mass ratio
E_sym = 31.3 / PhysicalConstants.hc # Symmetry energy

# Calculate coupling constants
println("Calculating coupling constants...")
couplings = calculate_couplings(ρ₀, B_A, K, m_ratio, E_sym)
fσ, fω, fρ, fδ, b, c = couplings

println("Coupling constants:")
println("  fσ = $(fσ)")
println("  fω = $(fω)")  
println("  fρ = $(fρ)")
println("  fδ = $(fδ)")
println("  b = $(b)")
println("  c = $(c)")

# Single point calculation
println("\n=== Single Point Calculation ===")

μ_B = 1001.0 / PhysicalConstants.hc  # Baryon chemical potential
println("Baryon chemical potential: μ_B = $(μ_B * PhysicalConstants.hc) MeV")

# Initial guess for field values and chemical potentials
x0 = [1.25, 0.01, μ_B/2, μ_B/2]  # [gσ, gδ, μ_p, μ_n]

println("Calculating pressure and derivatives...")
pressure, dpre_dmu1, dpre_dmu2, dpre_dmu3, dpre_dmu4 = 
    calculate_pressure_derivatives_efficient(μ_B, T, x0, nodes, couplings)

println("Results:")
println("  Pressure: P = $(pressure)")
println("  ∂P/∂μ₁ = $(dpre_dmu1)")
println("  ∂²P/∂μ² = $(dpre_dmu2)")
println("  ∂³P/∂μ³ = $(dpre_dmu3)")
println("  ∂⁴P/∂μ⁴ = $(dpre_dmu4)")

# Calculate thermodynamic fluctuations
println("\nCalculating thermodynamic fluctuations...")
kappa1, kappa2, kappa3, kappa4, fluctuation_ratios = 
    calculate_thermodynamic_fluctuations(μ_B, T, x0, nodes, couplings)

println("Cumulants:")
println("  κ₁ = $(kappa1)")
println("  κ₂ = $(kappa2)")
println("  κ₃ = $(kappa3)")  
println("  κ₄ = $(kappa4)")

println("Fluctuation ratios:")
println("  κ₂/κ₁ = $(fluctuation_ratios[1])")
println("  κ₃/κ₂ = $(fluctuation_ratios[2])")
println("  κ₄/κ₂ = $(fluctuation_ratios[3])")

# Batch calculation example
println("\n=== Batch Calculation Example ===")

# Create range of chemical potentials
μ_B_start = 1001.0 / PhysicalConstants.hc
μ_B_end = 600.0 / PhysicalConstants.hc
μ_B_step = -10.0 / PhysicalConstants.hc

μ_B_range = μ_B_start:μ_B_step:μ_B_end
println("Chemical potential range: $(μ_B_start*PhysicalConstants.hc) to $(μ_B_end*PhysicalConstants.hc) MeV")
println("Number of points: $(length(μ_B_range))")

# Run batch calculation
println("Running batch calculation...")
results = calculate_derivatives_batch(
    collect(μ_B_range), T, x0, nodes, couplings,
    save_results=true, 
    output_file=joinpath(@__DIR__, "..", "output", "gas_liquid_example.dat")
)

println("Batch calculation completed!")
println("Results saved to: gas_liquid_example.dat")

println("\n=== Example completed successfully! ===")
