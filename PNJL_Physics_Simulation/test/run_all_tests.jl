"""
Master test runner for all PNJL Physics Simulation tests.

This script runs all available tests and provides a summary.
"""

println("="^70)
println("PNJL Physics Simulation - Complete Test Suite")
println("="^70)

# Store test results
test_results = Dict{String, Bool}()

function run_test_file(test_name, test_file)
    """Run a test file and capture the result"""
    println("\n" * "="^50)
    println("Running: $test_name")
    println("="^50)
    
    try
        include(test_file)
        test_results[test_name] = true
        println("✅ $test_name: PASSED")
        return true
    catch e
        test_results[test_name] = false
        println("❌ $test_name: FAILED")
        println("Error: $e")
        # Print a few lines of traceback for debugging
        bt = catch_backtrace()
        for (i, frame) in enumerate(stacktrace(bt))
            if i > 3  # Limit traceback depth
                break
            end
            println("  at $frame")
        end
        return false
    end
end

# Get test directory
test_dir = @__DIR__

# List of tests to run
tests = [
    ("Basic Module Loading", "runtests.jl"),
    ("Module Loading Tests", "test_modules.jl"),
    ("Gas-Liquid Model", "test_gas_liquid.jl"),
    ("PNJL Model", "test_pnjl.jl"),
    ("PNJL Anisotropic Model", "test_pnjl_aniso.jl"),
    ("Rotation Model", "test_rotation.jl"),
    ("Debug Test", "debug_test.jl")
]

# Run all tests
total_tests = length(tests)
passed_tests = 0

for (test_name, test_file) in tests
    test_path = joinpath(test_dir, test_file)
    if isfile(test_path)
        if run_test_file(test_name, test_path)
            passed_tests += 1
        end
    else
        println("⚠️  Test file not found: $test_file")
        test_results[test_name] = false
    end
end

# Print summary
println("\n" * "="^70)
println("TEST SUMMARY")
println("="^70)

for (test_name, passed) in test_results
    status = passed ? "✅ PASSED" : "❌ FAILED"
    println("  $status - $test_name")
end

println("\nOverall Results:")
println("  Tests Passed: $passed_tests/$total_tests")
println("  Success Rate: $(round(passed_tests/total_tests*100, digits=1))%")

if passed_tests == total_tests
    println("\n🎉 All tests passed! The PNJL Physics Simulation package is working correctly.")
else
    println("\n⚠️  Some tests failed. Please check the errors above.")
    if passed_tests > 0
        println("Note: $(passed_tests) test(s) did pass, so basic functionality is working.")
    end
end

println("\n" * "="^70)
println("Available Models:")
println("  ✓ Core utilities (constants, integration, thermodynamics)")
println("  ✓ Gas-Liquid phase transition model")
println("  ✓ Standard PNJL model")
println("  ✓ PNJL Anisotropic model (momentum-dependent anisotropy)")
println("  ✓ Rotation model (Bessel function integration)")
println("="^70)
