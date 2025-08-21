"""
Test suite for the PNJL model.

Tests basic functionality including:
- Integration node generation
- Chiral condensate calculations
- Polyakov loop potential
- Pressure calculations
"""

using Test

# Get project root and load modules directly
project_root = abspath(joinpath(@__DIR__, ".."))
println("Loading PNJL model from: $project_root")

# Load core modules
include(joinpath(project_root, "src", "core", "constants.jl"))
include(joinpath(project_root, "src", "core", "integration_interface.jl"))
include(joinpath(project_root, "src", "core", "thermodynamics.jl"))

# Load PNJL model
include(joinpath(project_root, "src", "models", "pnjl", "constants.jl"))
include(joinpath(project_root, "src", "models", "pnjl", "functions.jl"))

# Import modules
using .PhysicalConstants
using .IntegrationInterface
using .PNJLConstants
using .PNJLFunctions

@testset "PNJL Model Tests" begin
    
    @testset "Integration Nodes" begin
        nodes = get_nodes(64)
        @test length(nodes) == 4  # [p_nodes1, p_nodes2, coef1, coef2]
        
        p_nodes1, p_nodes2, coef1, coef2 = nodes
        @test length(p_nodes1) == 64
        @test length(p_nodes2) == 64
        @test length(coef1) == 64
        @test length(coef2) == 64
        
        # Check that all nodes and coefficients are positive
        @test all(p_nodes1 .> 0)
        @test all(p_nodes2 .> 0)
        @test all(coef1 .> 0)
        @test all(coef2 .> 0)
    end
    
    @testset "Chiral Condensate" begin
        # Test chiral condensate calculation
        phi = [-0.1, -0.1, -1.7]  # [φ_u, φ_d, φ_s]
        
        chi = calculate_chiral(phi)
        @test isfinite(chi)
        @test isreal(chi)
        
        # Test with zero condensates
        phi_zero = [0.0, 0.0, 0.0]
        chi_zero = calculate_chiral(phi_zero)
        @test chi_zero == 0.0
    end
    
    @testset "Polyakov Loop Potential" begin
        T = 150.0 / PhysicalConstants.hc
        Phi1 = 0.5
        Phi2 = 0.5
        
        U = calculate_U(T, Phi1, Phi2)
        @test isfinite(U)
        @test isreal(U)
        
        # Test symmetry: U(Phi1, Phi2) = U(Phi2, Phi1)
        U_swap = calculate_U(T, Phi2, Phi1)
        @test abs(U - U_swap) < 1e-12
        
        # Test with zero Polyakov loops
        U_zero = calculate_U(T, 0.0, 0.0)
        @test isfinite(U_zero)
    end
    
    @testset "Effective Masses" begin
        phi = [-0.1, -0.1, -1.7]
        masses = calculate_mass_vec(phi)
        
        @test length(masses) == 3
        @test all(isfinite.(masses))
        @test all(isreal.(masses))
        
        # Masses should be positive for physical condensates
        @test all(masses .> 0)
    end
    
    @testset "Energy Calculation" begin
        mass = 0.3  # fm^-1
        p = 1.0     # fm^-1
        
        E = calculate_energy(mass, p)
        @test isfinite(E)
        @test isreal(E)
        @test E > 0
        
        # Check relativistic energy relation
        expected = sqrt(p^2 + mass^2)
        @test abs(E - expected) < 1e-12
        
        # Test limit: E → |p| as m → 0
        E_massless = calculate_energy(0.0, p)
        @test abs(E_massless - p) < 1e-12
    end
    
    @testset "Pressure Calculation" begin
        nodes = get_nodes(64)
        T = 150.0 / PhysicalConstants.hc
        phi = [-0.1, -0.1, -1.7]
        Phi1, Phi2 = 0.5, 0.5
        mu = [320.0/PhysicalConstants.hc, 320.0/PhysicalConstants.hc, 320.0/PhysicalConstants.hc]
        
        @test_nowarn pressure = calculate_pressure(phi, Phi1, Phi2, mu, T, nodes)
        
        pressure = calculate_pressure(phi, Phi1, Phi2, mu, T, nodes)
        @test isfinite(pressure)
        @test isreal(pressure)
    end
    
    @testset "Pressure Wrapper" begin
        nodes = get_nodes(64)
        T = 150.0 / PhysicalConstants.hc
        x = [-0.1, -0.1, -1.7, 0.5, 0.5]  # [φ_u, φ_d, φ_s, Φ₁, Φ₂]
        mu = [320.0/PhysicalConstants.hc, 320.0/PhysicalConstants.hc, 320.0/PhysicalConstants.hc]
        
        @test_nowarn pressure = pressure_wrapper(x, mu, T, nodes)
        
        pressure = pressure_wrapper(x, mu, T, nodes)
        @test isfinite(pressure)
        @test isreal(pressure)
    end
    
    @testset "Gradient Calculation" begin
        nodes = get_nodes(32)  # Smaller for faster testing
        T = 150.0 / PhysicalConstants.hc
        x = [-0.1, -0.1, -1.7, 0.5, 0.5]
        mu = [320.0/PhysicalConstants.hc, 320.0/PhysicalConstants.hc, 320.0/PhysicalConstants.hc]
        
        @test_nowarn grad = calculate_core(x, mu, T, nodes)
        
        grad = calculate_core(x, mu, T, nodes)
        @test length(grad) == 5
        @test all(isfinite.(grad))
        @test all(isreal.(grad))
    end
    
    @testset "Density Calculation" begin
        nodes = get_nodes(32)
        T = 150.0 / PhysicalConstants.hc
        x = [-0.1, -0.1, -1.7, 0.5, 0.5]
        mu = [320.0/PhysicalConstants.hc, 320.0/PhysicalConstants.hc, 320.0/PhysicalConstants.hc]
        
        @test_nowarn rho = calculate_rho(x, mu, T, nodes)
        
        rho = calculate_rho(x, mu, T, nodes)
        @test length(rho) == 3
        @test all(isfinite.(rho))
        @test all(isreal.(rho))
    end
    
    @testset "Thermodynamic Quantities" begin
        nodes = get_nodes(32)
        T = 150.0 / PhysicalConstants.hc
        x = [-0.1, -0.1, -1.7, 0.5, 0.5]
        mu = [320.0/PhysicalConstants.hc, 320.0/PhysicalConstants.hc, 320.0/PhysicalConstants.hc]
        
        @test_nowarn pressure, rho, entropy, energy = calculate_thermo(x, mu, T, nodes)
        
        pressure, rho, entropy, energy = calculate_thermo(x, mu, T, nodes)
        
        @test isfinite(pressure) && isreal(pressure)
        @test length(rho) == 3 && all(isfinite.(rho)) && all(isreal.(rho))
        @test isfinite(entropy) && isreal(entropy)
        @test isfinite(energy) && isreal(energy)
        
        # Basic thermodynamic consistency check
        # For this test, we just check that the thermodynamic relation holds approximately
        energy_check = -pressure + sum(mu .* rho) + T * entropy
        @test abs(energy - energy_check) < 1e-10
    end
    
end
