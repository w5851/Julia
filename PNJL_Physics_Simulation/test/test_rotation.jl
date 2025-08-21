"""
Test suite for the Rotation model.

Tests functionality including:
- Bessel function integration
- Rotation energy calculations
- Pressure calculations with angular velocity
"""

using Test

# Get project root and load modules directly
project_root = abspath(joinpath(@__DIR__, ".."))
println("Loading Rotation model from: $project_root")

# Load core modules
include(joinpath(project_root, "src", "core", "constants.jl"))
include(joinpath(project_root, "src", "core", "integration_interface.jl"))
include(joinpath(project_root, "src", "core", "thermodynamics.jl"))

# Load Rotation model
include(joinpath(project_root, "src", "models", "rotation", "constants.jl"))
include(joinpath(project_root, "src", "models", "rotation", "functions.jl"))

# Import modules
using .PhysicalConstants
using .IntegrationInterface
using .RotationConstants

@testset "Rotation Model Tests" begin
    
    @testset "Constants Loading" begin
        @test RotationConstants.hc > 0
        @test RotationConstants.G_f > 0
        @test RotationConstants.Lambda_f > 0
        @test RotationConstants.m0_q_f > 0
        @test RotationConstants.r0 > 0
        @test haskey(RotationConstants.coefficients, "a")
        @test haskey(RotationConstants.coefficients, "b")
        println("✓ All constants loaded correctly")
    end
    
    @testset "Bessel Function Initialization" begin
        p = [1.0, 2.0]
        theta = [π/4, π/2]
        n = [0, 1]
        w = [0.5, 0.8]
        
        coef = init_bessel(p, theta, n, w)
        @test length(coef) == length(p)
        @test all(isfinite.(coef))
        
        println("✓ Bessel function initialization works correctly")
    end
    
    @testset "Integration Nodes 3D" begin
        nodes1, nodes2 = get_nodes(4, 3)  # Small size for testing
        
        @test length(nodes1) == 3  # [p_mesh, n_mesh, coefficient]
        @test length(nodes2) == 3
        
        p_mesh1, n_mesh1, coef1 = nodes1
        p_mesh2, n_mesh2, coef2 = nodes2
        
        # Check dimensions: p_num × t_num × n_num = 4 × 3 × 11
        @test size(p_mesh1) == (4, 3, 11)
        @test size(n_mesh1) == (4, 3, 11) 
        @test size(coef1) == (4, 3, 11)
        
        @test all(p_mesh1 .> 0)
        @test all(-5 .<= n_mesh1 .<= 5)  # quantum numbers from -5 to 5
        @test all(isfinite.(coef1))
        
        println("✓ 3D node generation works correctly")
    end
    
    @testset "Energy Calculation with Rotation" begin
        mass = 0.3
        p = 1.0
        n = 2      # quantum number
        omega = 0.1  # angular velocity
        
        # Test energy calculation
        E = calculate_energy(mass, p, n, omega)
        expected_E = sqrt(p^2 + mass^2) - (0.5 + n) * omega
        @test abs(E - expected_E) < 1e-10
        
        # Test different quantum numbers
        E_minus = calculate_energy(mass, p, -1, omega)
        E_plus = calculate_energy(mass, p, 1, omega)
        
        # Higher quantum number should give lower energy (due to negative term)
        @test E_plus < E_minus
        
        println("✓ Energy calculation with rotation works correctly")
    end
    
    @testset "Mass Calculation" begin
        phi = 0.1
        mass = calculate_mass(phi)
        
        expected_mass = RotationConstants.m0_q_f - 2.0 * RotationConstants.G_f * phi
        @test abs(mass - expected_mass) < 1e-10
        @test mass > 0  # Mass should be positive for reasonable phi
        
        println("✓ Mass calculation works correctly")
    end
    
    @testset "Chiral Condensate (Rotation)" begin
        phi = 0.1
        chi = calculate_chiral(phi)
        
        expected_chi = RotationConstants.G_f * phi^2
        @test abs(chi - expected_chi) < 1e-10
        @test chi >= 0  # Should be non-negative for phi^2
        
        println("✓ Chiral condensate calculation works correctly")
    end
    
    @testset "Standard Polyakov Potential" begin
        T = 150 / RotationConstants.hc
        Phi1, Phi2 = 0.2, 0.3
        
        U = calculate_U(T, Phi1, Phi2)
        @test isa(U, Real)
        @test isfinite(U)
        
        println("✓ Standard Polyakov potential calculation works correctly")
    end
    
    @testset "Rotation-dependent Factors" begin
        T = 150 / RotationConstants.hc
        omega = 100 / RotationConstants.hc
        
        f, f_inv = calc_factors(T, omega)
        @test isa(f, Real)
        @test isa(f_inv, Real)
        @test isfinite(f)
        @test isfinite(f_inv)
        @test f * f_inv ≈ 1.0
        
        println("✓ Rotation-dependent factor calculation works correctly")
    end
    
    @testset "Rotation-dependent Polyakov Potential" begin
        T = 150 / RotationConstants.hc
        Phi1, Phi2 = 0.2, 0.3
        omega = 100 / RotationConstants.hc
        
        U_rot = calc_U(T, Phi1, Phi2, omega)
        @test isa(U_rot, Real)
        @test isfinite(U_rot)
        
        # Compare with omega = 0
        U_no_rot = calc_U(T, Phi1, Phi2, 0.0)
        @test U_rot != U_no_rot  # Should be different
        
        println("✓ Rotation-dependent Polyakov potential works correctly")
    end
    
    @testset "Pressure Calculation (Rotation)" begin
        # Set up test parameters
        phi = 0.1
        Phi1, Phi2 = 0.2, 0.3
        mu = 100 / RotationConstants.hc
        T = 150 / RotationConstants.hc
        omega = 100 / RotationConstants.hc
        
        # Generate small nodes for fast testing
        nodes1, nodes2 = get_nodes(2, 2)  # Very small for speed
        
        # Test pressure calculation
        pressure = calculate_pressure(phi, Phi1, Phi2, mu, T, nodes1, omega)
        @test isa(pressure, Real)
        @test isfinite(pressure)
        
        # Test wrapper function
        x = [phi, Phi1, Phi2]
        pressure_wrap = pressure_wrapper(x, mu, T, nodes1, omega)
        @test abs(pressure - pressure_wrap) < 1e-10
        
        println("✓ Pressure calculation works correctly")
    end
    
    @testset "Gradient Calculation (Rotation)" begin
        # Set up test parameters
        phi = 0.1
        Phi1, Phi2 = 0.2, 0.3
        mu = 100 / RotationConstants.hc
        T = 150 / RotationConstants.hc
        omega = 100 / RotationConstants.hc
        x = [phi, Phi1, Phi2]
        
        # Generate small nodes for fast testing
        nodes1, nodes2 = get_nodes(2, 2)
        
        # Test gradient calculation
        grad = calculate_core(x, mu, T, nodes1, omega)
        @test length(grad) == 3
        @test all(isfinite.(grad))
        
        println("✓ Gradient calculation works correctly")
    end
    
    @testset "Density Calculation (Rotation)" begin
        # Set up test parameters
        phi = 0.1
        Phi1, Phi2 = 0.2, 0.3
        mu = 100 / RotationConstants.hc
        T = 150 / RotationConstants.hc
        omega = 100 / RotationConstants.hc
        x = [phi, Phi1, Phi2]
        
        # Generate small nodes for fast testing
        nodes1, nodes2 = get_nodes(2, 2)
        
        # Test density calculation
        rho = calculate_rho(x, mu, T, nodes1, omega)
        @test isa(rho, Real)
        @test isfinite(rho)
        @test rho >= 0  # Density should be non-negative
        
        println("✓ Density calculation works correctly")
    end
    
    @testset "Thermodynamic Quantities (Rotation)" begin
        # Set up test parameters
        phi = 0.1
        Phi1, Phi2 = 0.2, 0.3
        mu = 100 / RotationConstants.hc
        T = 150 / RotationConstants.hc
        omega = 100 / RotationConstants.hc
        x = [phi, Phi1, Phi2]
        
        # Generate small nodes for fast testing
        nodes1, nodes2 = get_nodes(2, 2)
        
        # Test thermodynamic calculation
        pressure, rho, entropy, energy = calculate_thermo(x, mu, T, nodes1, omega)
        @test all(isfinite.([pressure, rho, entropy, energy]))
        @test rho >= 0
        
        println("✓ Thermodynamic quantities calculation works correctly")
    end
end

println("🎉 All Rotation model tests passed!")
