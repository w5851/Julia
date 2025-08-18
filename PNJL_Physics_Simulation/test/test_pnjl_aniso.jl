"""
Test suite for the PNJL Anisotropic model.

Tests functionality including:
- Integration node generation with anisotropy
- Anisotropic energy calculations
- Pressure calculations with xi parameter
"""

using Test

# Get project root and load modules directly
project_root = abspath(joinpath(@__DIR__, ".."))
println("Loading PNJL Aniso model from: $project_root")

# Load core modules
include(joinpath(project_root, "src", "core", "constants.jl"))
include(joinpath(project_root, "src", "core", "integration.jl"))
include(joinpath(project_root, "src", "core", "thermodynamics.jl"))

# Load PNJL Aniso model
include(joinpath(project_root, "src", "models", "pnjl_aniso", "constants.jl"))
include(joinpath(project_root, "src", "models", "pnjl_aniso", "functions.jl"))

# Import modules
using .PhysicalConstants
using .Integration
using .PNJLAnisoConstants

@testset "PNJL Anisotropic Model Tests" begin
    
    @testset "Constants Loading" begin
        @test PNJLAnisoConstants.hc > 0
        @test PNJLAnisoConstants.G_f > 0
        @test PNJLAnisoConstants.Lambda_f > 0
        @test PNJLAnisoConstants.m0_q_f > 0
        @test PNJLAnisoConstants.m0_s_f > 0
        println("✓ All constants loaded correctly")
    end
    
    @testset "Integration Nodes" begin
        nodes_1, nodes_2 = get_nodes(8, 4)  # Small size for testing
        
        @test length(nodes_1) == 3  # [p_mesh, t_mesh, coefficient]
        @test length(nodes_2) == 3
        
        p_mesh1, t_mesh1, coef1 = nodes_1
        p_mesh2, t_mesh2, coef2 = nodes_2
        
        @test size(p_mesh1) == (8, 4)  # p_num × t_num
        @test size(t_mesh1) == (8, 4)
        @test size(coef1) == (8, 4)
        
        @test all(p_mesh1 .> 0)
        @test all(0 .<= t_mesh1 .<= 1)  # cosθ ∈ [0,1]
        @test all(coef1 .> 0)
        
        println("✓ Node generation works correctly")
    end
    
    @testset "Energy Calculation with Anisotropy" begin
        mass = 0.3
        p = 1.0
        xi = 0.4  # anisotropy parameter
        t = 0.5   # cosθ
        
        # Test without anisotropy (xi=0)
        E0 = calculate_energy(mass, p, 0.0, t)
        expected_E0 = sqrt(p^2 + mass^2)
        @test abs(E0 - expected_E0) < 1e-10
        
        # Test with anisotropy
        E_xi = calculate_energy(mass, p, xi, t)
        expected_E_xi = sqrt(p^2 + mass^2 + xi * (p*t)^2)
        @test abs(E_xi - expected_E_xi) < 1e-10
        
        # Anisotropic energy should be larger
        @test E_xi > E0
        
        println("✓ Energy calculation with anisotropy works correctly")
    end
    
    @testset "Mass Vector Calculation" begin
        phi = [-0.1, -0.1, -1.7]
        masses = calculate_mass_vec(phi)
        
        @test length(masses) == 3
        @test eltype(masses) == eltype(phi)
        
        # Test that masses are positive (for reasonable phi values)
        @test all(masses .> 0)
        
        println("✓ Mass vector calculation works correctly")
    end
    
    @testset "Chiral Condensate" begin
        phi = [-0.1, -0.1, -1.7]
        chi = calculate_chiral(phi)
        
        @test isa(chi, Real)
        # For negative phi values, chiral term should be positive
        @test chi > 0
        
        println("✓ Chiral condensate calculation works correctly")
    end
    
    @testset "Polyakov Loop Potential" begin
        T = 150 / PNJLAnisoConstants.hc
        Phi1, Phi2 = 0.5, 0.5
        
        U = calculate_U(T, Phi1, Phi2)
        @test isa(U, Real)
        @test isfinite(U)
        
        println("✓ Polyakov loop potential calculation works correctly")
    end
    
    @testset "Pressure Calculation" begin
        # Set up test parameters
        phi = [-0.1, -0.1, -1.7]
        Phi1, Phi2 = 0.5, 0.5
        mu = [320/PNJLAnisoConstants.hc, 320/PNJLAnisoConstants.hc, 320/PNJLAnisoConstants.hc]
        T = 150 / PNJLAnisoConstants.hc
        xi = 0.4
        
        # Generate small nodes for fast testing
        nodes_1, nodes_2 = get_nodes(4, 2)
        
        # Test pressure calculation
        pressure = calculate_pressure(phi, Phi1, Phi2, mu, T, nodes_1, nodes_2, xi)
        @test isa(pressure, Real)
        @test isfinite(pressure)
        
        # Test wrapper function
        x = [phi..., Phi1, Phi2]
        pressure_wrap = pressure_wrapper(x, mu, T, nodes_1, nodes_2, xi)
        @test abs(pressure - pressure_wrap) < 1e-10
        
        println("✓ Pressure calculation works correctly")
    end
    
    @testset "Gradient Calculation" begin
        # Set up test parameters
        phi = [-0.1, -0.1, -1.7]
        Phi1, Phi2 = 0.5, 0.5
        mu = [320/PNJLAnisoConstants.hc, 320/PNJLAnisoConstants.hc, 320/PNJLAnisoConstants.hc]
        T = 150 / PNJLAnisoConstants.hc
        xi = 0.4
        x = [phi..., Phi1, Phi2]
        
        # Generate small nodes for fast testing
        nodes_1, nodes_2 = get_nodes(4, 2)
        
        # Test gradient calculation
        grad = calculate_core(x, mu, T, nodes_1, nodes_2, xi)
        @test length(grad) == 5
        @test all(isfinite.(grad))
        
        println("✓ Gradient calculation works correctly")
    end
    
    @testset "Density Calculation" begin
        # Set up test parameters
        phi = [-0.1, -0.1, -1.7]
        Phi1, Phi2 = 0.5, 0.5
        mu = [320/PNJLAnisoConstants.hc, 320/PNJLAnisoConstants.hc, 320/PNJLAnisoConstants.hc]
        T = 150 / PNJLAnisoConstants.hc
        xi = 0.4
        x = [phi..., Phi1, Phi2]
        
        # Generate small nodes for fast testing
        nodes_1, nodes_2 = get_nodes(4, 2)
        
        # Test density calculation
        rho = calculate_rho(x, mu, T, nodes_1, nodes_2, xi)
        @test length(rho) == 3
        @test all(isfinite.(rho))
        @test all(rho .>= 0)  # Densities should be non-negative
        
        println("✓ Density calculation works correctly")
    end
end

println("🎉 All PNJL Anisotropic model tests passed!")
