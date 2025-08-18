println("Debug: Starting test...")

# Print current directory
println("Current directory: ", pwd())

# Print test directory
test_dir = @__DIR__
println("Test directory: ", test_dir)

# Print project root
project_root = dirname(test_dir)
println("Project root: ", project_root)

# Check if constants file exists
constants_path = joinpath(project_root, "src", "core", "constants.jl")
println("Constants path: ", constants_path)
println("Constants file exists: ", isfile(constants_path))

# Try to load it
if isfile(constants_path)
    try
        include(constants_path)
        println("✓ Successfully included constants.jl")
        
        # Try to use the module
        using .PhysicalConstants
        println("✓ Successfully imported PhysicalConstants")
        println("  π = ", PhysicalConstants.π)
        println("  hc = ", PhysicalConstants.hc)
    catch e
        println("✗ Error: ", e)
    end
else
    println("✗ Constants file not found!")
end
