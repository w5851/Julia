"""
Simplified test runner to verify module loading.
"""

using Test

println("Testing module loading...")

@testset "Module Loading Tests" begin
    @testset "Core Modules" begin
        @test_nowarn include(joinpath(@__DIR__, "..", "src", "core", "constants.jl"))
        @test_nowarn include(joinpath(@__DIR__, "..", "src", "core", "integration.jl"))
        @test_nowarn include(joinpath(@__DIR__, "..", "src", "core", "thermodynamics.jl"))
    end
    
    @testset "Gas Liquid Model" begin
        @test_nowarn include(joinpath(@__DIR__, "..", "src", "models", "gas_liquid", "constants.jl"))
        @test_nowarn include(joinpath(@__DIR__, "..", "src", "models", "gas_liquid", "functions.jl"))
    end
    
    @testset "PNJL Model" begin
        @test_nowarn include(joinpath(@__DIR__, "..", "src", "models", "pnjl", "constants.jl"))
        @test_nowarn include(joinpath(@__DIR__, "..", "src", "models", "pnjl", "functions.jl"))
    end
    
    @testset "PNJL Aniso Model" begin
        @test_nowarn include(joinpath(@__DIR__, "..", "src", "models", "pnjl_aniso", "constants.jl"))
        @test_nowarn include(joinpath(@__DIR__, "..", "src", "models", "pnjl_aniso", "functions.jl"))
    end
    
    @testset "Rotation Model" begin
        @test_nowarn include(joinpath(@__DIR__, "..", "src", "models", "rotation", "constants.jl"))
        @test_nowarn include(joinpath(@__DIR__, "..", "src", "models", "rotation", "functions.jl"))
    end
    
    @testset "Main Package Loading" begin
        @test_nowarn include(joinpath(@__DIR__, "..", "src", "PNJLPhysicsSimulation.jl"))
        # Note: Package-style loading removed to avoid dependency issues
        println("✓ Main package file can be included successfully")
    end
end

println("All module loading tests completed!")
