"""
Test suite for the Gas-Liquid model.

Tests basic functionality including:
- Integration node generation
- Coupling constant calculation  
- Pressure calculations
- Derivative calculations
"""

using Test

# Get project root and load modules directly
project_root = abspath(joinpath(@__DIR__, ".."))
println("Loading Gas-Liquid model from: $project_root")

# Load core modules
include(joinpath(project_root, "src", "core", "constants.jl"))
include(joinpath(project_root, "src", "core", "integration.jl"))
include(joinpath(project_root, "src", "core", "thermodynamics.jl"))

# Load Gas-Liquid model
include(joinpath(project_root, "src", "models", "gas_liquid", "constants.jl"))
include(joinpath(project_root, "src", "models", "gas_liquid", "functions.jl"))

# Import modules
using .PhysicalConstants
using .Integration
using .GasLiquidConstants
using .GasLiquidFunctions

@testset "Gas-Liquid Model Tests" begin
    
    @testset "Integration Nodes" begin
        nodes = GasLiquidFunctions.get_nodes(64)
        p_nodes, coefficients = nodes
        
        @test length(p_nodes) == 64
        @test length(coefficients) == 64
        @test all(p_nodes .> 0)
        @test all(coefficients .> 0)
        
        # Test integration accuracy with simple function
        # ∫₀²⁰ p² dp = 8000/3
        exact = 8000.0/3
        numerical = sum(p_nodes.^2 .* coefficients) * π^2
        @test abs(numerical - exact) / exact < 0.01  # 1% accuracy
    end
    
    @testset "Coupling Constants" begin
        # Standard nuclear matter parameters
        ρ₀ = 0.16
        B_A = -16.0 / PhysicalConstants.hc
        K = 240.0 / PhysicalConstants.hc
        m_ratio = 0.75
        E_sym = 31.3 / PhysicalConstants.hc
        
        couplings = GasLiquidFunctions.calculate_couplings(ρ₀, B_A, K, m_ratio, E_sym)
        fσ, fω, fρ, fδ, b, c = couplings
        
        # Check that couplings are physical (finite and real)
        @test isfinite(fσ) && isreal(fσ)
        @test isfinite(fω) && isreal(fω)
        @test isfinite(fρ) && isreal(fρ)
        @test isfinite(fδ) && isreal(fδ)
        @test isfinite(b) && isreal(b)
        @test isfinite(c) && isreal(c)
        
        # Check signs (these should be consistent with model physics)
        @test fσ > 0  # Attractive scalar interaction
        @test fω > 0  # Repulsive vector interaction
    end
    
    @testset "Pressure Calculation" begin
        nodes = GasLiquidFunctions.get_nodes(128)
        T = 50.0 / PhysicalConstants.hc
        couplings = GasLiquidFunctions.calculate_couplings(0.16, -16.0/PhysicalConstants.hc, 
                                       240.0/PhysicalConstants.hc, 0.75, 31.3/PhysicalConstants.hc)
        
        μ_B = 500.0 / PhysicalConstants.hc
        x0 = [1.0, 0.0, μ_B/2, μ_B/2]
        
        # Test that pressure calculation doesn't throw errors
        @test_nowarn pressure = calculate_pressure_solved(μ_B, T, x0, nodes, couplings)
        
        pressure = calculate_pressure_solved(μ_B, T, x0, nodes, couplings)
        @test isfinite(pressure)
        @test isreal(pressure)
    end
    
    @testset "Pressure Derivatives" begin
        nodes = GasLiquidFunctions.get_nodes(64)  # Smaller for faster testing
        T = 50.0 / PhysicalConstants.hc
        couplings = GasLiquidFunctions.calculate_couplings(0.16, -16.0/PhysicalConstants.hc, 
                                       240.0/PhysicalConstants.hc, 0.75, 31.3/PhysicalConstants.hc)
        
        μ_B = 500.0 / PhysicalConstants.hc
        x0 = [1.0, 0.0, μ_B/2, μ_B/2]
        
        # Test derivative calculation
        @test_nowarn pressure, dpre_dmu1, dpre_dmu2, dpre_dmu3, dpre_dmu4 = 
            GasLiquidFunctions.calculate_pressure_derivatives_efficient(μ_B, T, x0, nodes, couplings)
            
        pressure, dpre_dmu1, dpre_dmu2, dpre_dmu3, dpre_dmu4 = 
            GasLiquidFunctions.calculate_pressure_derivatives_efficient(μ_B, T, x0, nodes, couplings)
            
        # Check that all derivatives are finite
        @test isfinite(pressure)
        @test isfinite(dpre_dmu1) 
        @test isfinite(dpre_dmu2)
        @test isfinite(dpre_dmu3)
        @test isfinite(dpre_dmu4)
        
        # Physical check: first derivative should be positive (∂P/∂μ ~ ρ > 0)
        @test dpre_dmu1 > 0
    end
    
    @testset "Thermodynamic Fluctuations" begin
        nodes = GasLiquidFunctions.get_nodes(64)
        T = 50.0 / PhysicalConstants.hc
        couplings = GasLiquidFunctions.calculate_couplings(0.16, -16.0/PhysicalConstants.hc, 
                                       240.0/PhysicalConstants.hc, 0.75, 31.3/PhysicalConstants.hc)
        
        μ_B = 500.0 / PhysicalConstants.hc
        x0 = [1.0, 0.0, μ_B/2, μ_B/2]
        
        @test_nowarn kappa1, kappa2, kappa3, kappa4, fluctuation_ratios = 
            GasLiquidFunctions.calculate_thermodynamic_fluctuations(μ_B, T, x0, nodes, couplings)
            
        kappa1, kappa2, kappa3, kappa4, fluctuation_ratios = 
            GasLiquidFunctions.calculate_thermodynamic_fluctuations(μ_B, T, x0, nodes, couplings)
            
        # Check that cumulants are finite
        @test isfinite(kappa1)
        @test isfinite(kappa2)
        @test isfinite(kappa3)
        @test isfinite(kappa4)
        
        # Check fluctuation ratios
        @test length(fluctuation_ratios) == 3
        @test all(isfinite.(fluctuation_ratios))
    end
    
end
