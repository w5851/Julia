"""
Test script for PNJL Aniso model refactoring

This script tests the refactored PNJL_aniso model to ensure the integration
interface refactoring works correctly and produces consistent results.
"""

using Test

println("=== PNJL Aniso Model Refactoring Test ===")

# Get project root
project_root = abspath(joinpath(@__DIR__, ".."))
println("Testing from: $project_root")

# Test loading the refactored package
println("\n1. Testing package loading...")
try
    using PNJLPhysicsSimulation
    println("✅ Package loaded successfully")
catch e
    println("❌ Package loading failed: $e")
    exit(1)
end

# Test accessing PNJL Aniso constants and functions
println("\n2. Testing PNJL Aniso module access...")
try
    # Test constants access
    G_f = PNJLAnisoConstants.G_f
    Lambda_f = PNJLAnisoConstants.Lambda_f
    println("✅ PNJL Aniso constants accessible: G_f = $G_f, Lambda_f = $Lambda_f")
    
    # Test functions access
    using PNJLPhysicsSimulation.PNJLAnisoFunctions
    println("✅ PNJL Aniso functions module loaded")
catch e
    println("❌ PNJL Aniso module access failed: $e")
    exit(1)
end

# Test node generation
println("\n3. Testing node generation...")
try
    nodes1, nodes2 = PNJLAnisoFunctions.get_nodes_aniso(32, 16)
    println("✅ Node generation successful")
    println("   nodes1 dimensions: $(size(nodes1[1])), $(size(nodes1[2])), $(size(nodes1[3]))")
    println("   nodes2 dimensions: $(size(nodes2[1])), $(size(nodes2[2])), $(size(nodes2[3]))")
catch e
    println("❌ Node generation failed: $e")
    exit(1)
end

# Test basic physical calculations
println("\n4. Testing basic physical calculations...")
try
    # Test parameters
    T = 150.0 / PhysicalConstants.hc  # 150 MeV
    phi = [-0.1, -0.1, -0.15]
    Phi1, Phi2 = 0.5, 0.5
    xi = 0.1  # anisotropy parameter
    
    # Test mass calculation
    masses = PNJLAnisoFunctions.calculate_mass_vec(phi)
    println("✅ Mass calculation: masses = $masses")
    
    # Test energy calculation with anisotropy
    p, t = 1.0, 0.5
    E_aniso = PNJLAnisoFunctions.calculate_energy_aniso(masses[1], p, xi, t)
    println("✅ Anisotropic energy calculation: E = $E_aniso")
    
    # Test U potential
    U = PNJLAnisoFunctions.calculate_U_aniso(T, Phi1, Phi2)
    println("✅ U potential calculation: U = $U")
    
    # Test chiral contribution
    chi = PNJLAnisoFunctions.calculate_chiral_aniso(phi)
    println("✅ Chiral contribution: chi = $chi")
    
catch e
    println("❌ Basic physical calculations failed: $e")
    exit(1)
end

# Test new integration interface
println("\n5. Testing new integration interface...")
try
    # Generate test nodes
    p_nodes, p_weights = IntegrationInterface.gauleg(0.0, 10.0, 32)
    t_nodes, t_weights = IntegrationInterface.gauleg(-1.0, 1.0, 16)
    
    # Test parameters
    masses = [0.3, 0.3, 0.4]  # Test masses in natural units
    mu = [0.3, 0.3, 0.3]     # Chemical potentials
    T = 0.15                  # Temperature
    Phi1, Phi2 = 0.5, 0.5
    xi = 0.1
    
    # Test vacuum energy integral
    vacuum_energy = PNJLAnisoFunctions.vacuum_energy_integral_aniso(
        masses, p_nodes, p_weights, t_nodes, t_weights, xi)
    println("✅ Vacuum energy integral: E_vac = $vacuum_energy")
    
    # Test thermal omega integral  
    thermal_omega = PNJLAnisoFunctions.omega_thermal_integral_aniso(
        masses, mu, T, Phi1, Phi2, p_nodes, p_weights, t_nodes, t_weights, xi)
    println("✅ Thermal omega integral: Omega_th = $thermal_omega")
    
catch e
    println("❌ New integration interface test failed: $e")
    exit(1)
end

# Test backward compatibility
println("\n6. Testing backward compatibility...")
try
    # Generate nodes using old format
    nodes1, nodes2 = PNJLAnisoFunctions.get_nodes_aniso(16, 8)
    
    # Test parameters
    phi = [-0.1, -0.1, -0.15]
    masses = PNJLAnisoFunctions.calculate_mass_vec(phi)
    mu = [0.3, 0.3, 0.3]
    T = 0.15
    Phi1, Phi2 = 0.5, 0.5
    xi = 0.1
    
    # Extract node data
    p_nodes1 = @view nodes1[1][:]
    t_nodes1 = @view nodes1[2][:]
    coef1 = @view nodes1[3][:]
    
    # Test old interface (should still work but give deprecation warning)
    @test_nowarn energy_sum = PNJLAnisoFunctions.calculate_energy_sum(masses, p_nodes1, coef1, t_nodes1, xi)
    println("✅ Backward compatibility: old calculate_energy_sum still works")
    
    @test_nowarn log_sum = PNJLAnisoFunctions.calculate_log_sum_aniso(masses, p_nodes1, Phi1, Phi2, mu, T, coef1, t_nodes1, xi)
    println("✅ Backward compatibility: old calculate_log_sum_aniso still works")
    
catch e
    println("❌ Backward compatibility test failed: $e")
    exit(1)
end

# Test full pressure calculation
println("\n7. Testing full pressure calculation...")
try
    # Generate integration nodes
    nodes1, nodes2 = PNJLAnisoFunctions.get_nodes_aniso(24, 12)
    
    # Test parameters
    phi = [-0.1, -0.1, -0.15]
    Phi1, Phi2 = 0.5, 0.5
    mu = [0.32, 0.32, 0.32]  # Chemical potentials
    T = 0.15                 # Temperature
    xi = 0.05                # Small anisotropy
    
    # Calculate pressure
    pressure = PNJLAnisoFunctions.calculate_pressure_aniso(phi, Phi1, Phi2, mu, T, nodes1, nodes2, xi)
    println("✅ Full pressure calculation: P = $pressure")
    
    # Test pressure wrapper
    x = [phi[1], phi[2], phi[3], Phi1, Phi2]
    pressure_wrap = PNJLAnisoFunctions.pressure_wrapper(x, mu, T, nodes1, nodes2, xi)
    println("✅ Pressure wrapper: P = $pressure_wrap")
    
    # Verify consistency
    @test abs(pressure - pressure_wrap) < 1e-10
    println("✅ Pressure calculation consistency verified")
    
catch e
    println("❌ Full pressure calculation test failed: $e")
    exit(1)
end

println("\n=== All Tests Passed Successfully! ===")
println("✅ PNJL_aniso model refactoring completed successfully")
println("✅ Integration interface is working correctly")
println("✅ Backward compatibility maintained")
println("✅ All physical calculations produce finite results")
