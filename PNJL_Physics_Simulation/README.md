# PNJL Physics Simulation

A comprehensive Julia package for studying QCD phase transitions using various effective models.

## Overview

This package provides implementations of multiple physics models for studying the QCD phase diagram:

- **Gas-Liquid Model**: Nuclear liquid-gas phase transition
- **PNJL Model**: Polyakov-Nambu-Jona-Lasinio model for chiral symmetry breaking
- **PNJL Anisotropic Model**: PNJL with anisotropic momentum effects  
- **Rotation Model**: PNJL with rotation effects

## Project Structure

```
PNJL_Physics_Simulation/
├── src/
│   ├── PNJLPhysicsSimulation.jl    # Main module
│   ├── core/                        # Shared utilities
│   │   ├── constants.jl            # Physical constants
│   │   ├── integration.jl          # Numerical integration
│   │   └── thermodynamics.jl       # Thermodynamic functions
│   ├── models/                      # Physics models
│   │   ├── gas_liquid/             # Gas-liquid transition
│   │   ├── pnjl/                   # PNJL model
│   │   ├── pnjl_aniso/             # Anisotropic PNJL
│   │   └── rotation/               # Rotation effects
│   └── physics/                     # High-level physics calculations
├── test/                           # Test suite
├── docs/                           # Documentation
├── output/                         # Calculation results
├── scripts/                        # Example scripts
└── Project.toml                    # Package dependencies
```

## Features

### Modular Design
- Each physics model is implemented as an independent module
- Shared utilities for common calculations
- Clean separation of constants and functions

### Professional Code Organization
- Clear module boundaries with well-defined interfaces
- Comprehensive documentation and type annotations
- Optimized numerical algorithms with minimal allocations

### Extensible Architecture
- Easy to add new physics models
- Consistent API across different models
- Support for automatic differentiation and optimization

## Quick Start

```julia
using PNJLPhysicsSimulation

# Gas-liquid model example
using PNJLPhysicsSimulation.GasLiquidModel

# Set up calculation parameters
nodes = get_nodes(256)
T = 50.0 / PhysicalConstants.hc
couplings = calculate_couplings(0.16, -16.0/PhysicalConstants.hc, 
                               240.0/PhysicalConstants.hc, 0.75, 31.3/PhysicalConstants.hc)

# Calculate pressure and derivatives
μ_B = 1001.0 / PhysicalConstants.hc
x0 = [1.25, 0.01, μ_B/2, μ_B/2]
pressure, dpre_dmu1, dpre_dmu2, dpre_dmu3, dpre_dmu4 = 
    calculate_pressure_derivatives_efficient(μ_B, T, x0, nodes, couplings)

# PNJL model example  
using PNJLPhysicsSimulation.PNJLModel

nodes = get_nodes(128)
T = 150.0 / PhysicalConstants.hc
x = [-0.1, -0.1, -1.7, 0.5, 0.5]
mu = [320.0/PhysicalConstants.hc, 320.0/PhysicalConstants.hc, 320.0/PhysicalConstants.hc]

pressure = pressure_solve_core(x, mu, T, nodes)
```

## Dependencies

- Julia 1.6+
- FastGaussQuadrature.jl: Gaussian quadrature integration
- ForwardDiff.jl: Automatic differentiation
- NLsolve.jl: Nonlinear equation solving
- FiniteDifferences.jl: Numerical derivatives
- BenchmarkTools.jl: Performance benchmarking
- StaticArrays.jl: Fixed-size arrays
- SpecialFunctions.jl: Mathematical special functions

## Installation

```julia
using Pkg
Pkg.develop(path="path/to/PNJL_Physics_Simulation")
```

## Development Status

- ✅ Core infrastructure
- ✅ Gas-Liquid model (complete)  
- ✅ PNJL model (complete)
- 🚧 PNJL Anisotropic model (in progress)
- 🚧 Rotation model (in progress)
- 📋 Advanced physics calculations (planned)

## Contributing

This package follows professional Julia development practices:

1. All modules are self-contained with clear dependencies
2. Functions include comprehensive documentation
3. Performance-critical code is optimized for minimal allocations
4. Consistent naming conventions across all modules
5. Extensive test coverage for numerical accuracy

## License

[Add appropriate license information]
