"""
Core integration utilities for physics simulations.

This module provides common integration functions used across different
physics models in the PNJL simulation package.
"""
module Integration

using FastGaussQuadrature

export gauleg

"""
    gauleg(a, b, n) -> (nodes, weights)

Gauss-Legendre quadrature integration nodes and weights.

# Arguments
- `a::Real`: Lower integration bound
- `b::Real`: Upper integration bound  
- `n::Int`: Number of integration points

# Returns
- `nodes::Vector{Float64}`: Integration nodes
- `weights::Vector{Float64}`: Integration weights

# Example
```julia
nodes, weights = gauleg(0.0, 1.0, 10)
```
"""
function gauleg(a, b, n)
    t_nodes, t_weights = gausslegendre(n)
    nodes   = @. (b - a)/2 * t_nodes + (a + b)/2
    weights = @. (b - a)/2 * t_weights
    return nodes, weights
end

end  # module Integration
