# PNJL Physics Simulation - AI Copilot Instructions

## Project Overview
This is a Julia package for quantum chromodynamics (QCD) phase transition simulations using multiple effective physics models (Gas-Liquid, PNJL, PNJL-Aniso, Rotation). The codebase implements sophisticated numerical methods for calculating thermodynamic properties and phase diagrams.

## Architecture & Design Patterns

### Module Hierarchy (Read Multiple Files to Understand)
```
Core Layer: src/core/ (constants, integration, math_utils, thermodynamics)
↓
Models Layer: src/models/{gas_liquid,pnjl,pnjl_aniso,rotation}/
↓
Physics Layer: src/physics/ (high-level physics calculations)
```

**Critical Pattern**: Each model has separate `constants.jl` + `functions.jl` files. Constants are model-specific and isolated to avoid conflicts (eliminated duplicate physical constants like π, hc, Nc across models).

### Integration Interface Architecture
The project recently implemented a major refactor around `src/core/integration_interface.jl`:
- **Before**: `calculate_log_sum` functions duplicated across 3 models
- **After**: Unified `omega_thermal_integral()` and `vacuum_energy_integral()` functions
- **Pattern**: Use `IntegrationMethod` and `IntegrationGrid` abstractions for numerical integration

### Safety-First Numerical Computing
All mathematical functions use "safe" wrappers to prevent numerical instability:
```julia
# Use safe_log(x) instead of log(x) to handle negative/zero values
# Use safe_exp(x) instead of exp(x) to prevent overflow
# Pattern appears in src/core/math_utils.jl
```

## Developer Workflows

### Testing & Validation Commands
```bash
# Run comprehensive test suite
julia test/run_all_tests.jl

# Quick integration interface tests (36 tests)  
julia --project=. -e "include(\"test/test_integration_interface.jl\")"

# Performance validation
julia scripts/setup.jl  # Install dependencies first
julia scripts/gas_liquid_example.jl  # Example calculation
```

### Project Activation Pattern
Always use this sequence for loading the project:
```julia
cd("path/to/PNJL_Physics_Simulation")
push!(LOAD_PATH, "src")
include("src/PNJLPhysicsSimulation.jl")
using .PNJLPhysicsSimulation
```

### Documentation-Driven Development
**Critical**: Before making changes, read in order:
1. `agent/requirements.md` - Current status and TODOs
2. `agent/architecture.md` - Design principles  
3. `agent/api_reference.md` - Interface contracts
4. `USAGE_GUIDE.md` - Working examples

After changes, update requirements status and `agent/changelog.md`.

## Project-Specific Conventions

### Physical Units Convention
- All calculations use natural units with `hc = 197.33` MeV⋅fm conversion
- Temperature example: `T = 50.0 / PhysicalConstants.hc` (50 MeV)
- Chemical potential: `μ_B = 1001.0 / PhysicalConstants.hc`

### Function Naming Patterns
- `calculate_*`: Core physics calculations
- `get_*`: Utility functions (e.g., `get_nodes(256)`)  
- `safe_*`: Numerically stable math functions
- `*_integral`: Integration-specific functions using the new interface

### Integration Grid System
Models use consistent grid generation:
```julia
nodes = get_nodes(256)  # 256-point Gaussian quadrature
grid = create_momentum_grid(n_points, cutoff)
result = omega_thermal_integral(masses, mu, T, Phi1, Phi2, grid, method)
```

### Error Handling Strategy
- Prefer safe math functions over try/catch blocks
- Functions return meaningful error messages with physics context
- Boundary conditions handled at the mathematical level, not application level

## Key Integration Points

### Cross-Model Dependencies
Models share integration methods but maintain isolated constants. The `IntegrationInterface` module provides the bridge between physics models and numerical methods.

### External Dependencies Patterns
- `FastGaussQuadrature.jl`: Always use via `gauleg()` wrapper in `Integration` module
- `ForwardDiff.jl`: Used for automatic differentiation in pressure calculations
- `NLsolve.jl`: Nonlinear equation solving in phase transition calculations

## AI Assistant Integration

### Context Files Priority
1. **Always read first**: `agent/requirements.md` - Current development status and priorities
2. **Architecture guidance**: `agent/architecture.md` - Design patterns and module structure
3. **API contracts**: `agent/api_reference.md` - Interface definitions and usage examples
4. **Development workflow**: `agent/prompt.md` - Coding standards and best practices

### Multi-Platform AI Configuration
- **GitHub Copilot**: Uses this file for contextual suggestions
- **Cursor**: Reads `.cursorrules` for editor-specific guidance  
- **Windsurf**: Uses `.windsurfrules` for development workflow
- **Generic LLM**: Can reference `agent/prompt.md` for comprehensive guidance

### Task Management Integration
Before implementing features, check `agent/requirements.md` for:
- Current task status (✅ completed, 🔄 in progress, 📋 planned)
- Dependencies between tasks
- Priority levels and expected completion times
- Technical debt and known issues

## Current Development Status (2025-08-18)
- ✅ Constants separation refactor completed
- ✅ Integration interface implemented with full test coverage (36/36 tests pass)
- 🔄 Model refactoring to use new integration interface (in progress)
- 📋 Next: PNJL, PNJL-Aniso, Rotation model integration calls need updating
