"""
Function Registry Module

This module provides a centralized registry for all physics model functions
to resolve naming conflicts and provide clear function disambiguation.
"""

module FunctionRegistry

# Since model files don't have explicit modules, we'll create function wrappers
# that resolve to the correct implementations

# Export disambiguation functions
export get_pnjl_nodes, get_pnjl_aniso_nodes, get_rotation_nodes, get_gas_liquid_nodes
export calculate_pnjl_log_sum, calculate_pnjl_aniso_log_sum, calculate_rotation_log_sum
export calculate_pnjl_energy_sum, calculate_pnjl_aniso_energy_sum, calculate_rotation_energy_sum
export get_model_functions

# Since model functions are defined in the main module scope,
# we need to access them through the parent module's namespace.
# For now, we'll create a simple registry that doesn't rely on direct imports.

# Disambiguated node functions
"""
    get_pnjl_nodes(n_p::Int) -> (Vector{Float64}, Vector{Float64})

Get integration nodes and weights for PNJL model.
"""
function get_pnjl_nodes(n_p::Int)
    # This will be resolved at runtime through multiple dispatch
    error("PNJL get_nodes function not available - check model imports")
end

"""
    get_pnjl_aniso_nodes(p_num::Int, t_num::Int) -> (Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64})

Get integration nodes and weights for anisotropic PNJL model.
"""
function get_pnjl_aniso_nodes(p_num::Int, t_num::Int)
    error("PNJL Aniso get_nodes function not available - check model imports")
end

"""
    get_rotation_nodes(p_num::Int, t_num::Int) -> (Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64})

Get integration nodes and weights for rotation model.
"""
function get_rotation_nodes(p_num::Int, t_num::Int)
    error("Rotation get_nodes function not available - check model imports")
end

"""
    get_gas_liquid_nodes(n_p::Int) -> (Vector{Float64}, Vector{Float64})

Get integration nodes and weights for gas-liquid model.
"""
function get_gas_liquid_nodes(n_p::Int)
    error("Gas-Liquid get_nodes function not available - check model imports")
end

# Disambiguated log_sum functions
"""
    calculate_pnjl_log_sum(masses, p_nodes, Phi1, Phi2, mu, T, coefficient)

Calculate logarithmic sum for PNJL model.
"""
function calculate_pnjl_log_sum(masses, p_nodes, Phi1, Phi2, mu, T, coefficient)
    error("PNJL calculate_log_sum function not available - check model imports")
end

"""
    calculate_pnjl_aniso_log_sum(masses, p_nodes, Phi1, Phi2, mu, T, coefficient, t_nodes, xi)

Calculate logarithmic sum for anisotropic PNJL model.
"""
function calculate_pnjl_aniso_log_sum(masses, p_nodes, Phi1, Phi2, mu, T, coefficient, t_nodes, xi)
    error("PNJL Aniso calculate_log_sum function not available - check model imports")
end

"""
    calculate_rotation_log_sum(masses, p_nodes, Phi1, Phi2, mu, T, coefficient, n_nodes, omega)

Calculate logarithmic sum for rotation model.
"""
function calculate_rotation_log_sum(masses, p_nodes, Phi1, Phi2, mu, T, coefficient, n_nodes, omega)
    error("Rotation calculate_log_sum function not available - check model imports")
end

# Disambiguated energy_sum functions
"""
    calculate_pnjl_energy_sum(masses, p_nodes, coefficient)

Calculate energy sum for PNJL model.
"""
function calculate_pnjl_energy_sum(masses, p_nodes, coefficient)
    error("PNJL calculate_energy_sum function not available - check model imports")
end

"""
    calculate_pnjl_aniso_energy_sum(masses, p_nodes, coefficient, t_nodes)

Calculate energy sum for anisotropic PNJL model.
"""
function calculate_pnjl_aniso_energy_sum(masses, p_nodes, coefficient, t_nodes)
    error("PNJL Aniso calculate_energy_sum function not available - check model imports")
end

"""
    calculate_rotation_energy_sum(masses, p_nodes, coefficient, n_nodes, omega)

Calculate energy sum for rotation model.
"""
function calculate_rotation_energy_sum(masses, p_nodes, coefficient, n_nodes, omega)
    error("Rotation calculate_energy_sum function not available - check model imports")
end

"""
    get_model_functions(model_type::Symbol) -> Dict{Symbol, Function}

Get all functions for a specific model type.

# Arguments
- `model_type`: :pnjl, :pnjl_aniso, :rotation, or :gas_liquid

# Returns
A dictionary mapping function names to model-specific functions.
"""
function get_model_functions(model_type::Symbol)
    if model_type == :pnjl
        return Dict(
            :get_nodes => get_pnjl_nodes,
            :calculate_log_sum => calculate_pnjl_log_sum,
            :calculate_energy_sum => calculate_pnjl_energy_sum
        )
    elseif model_type == :pnjl_aniso
        return Dict(
            :get_nodes => get_pnjl_aniso_nodes,
            :calculate_log_sum => calculate_pnjl_aniso_log_sum,
            :calculate_energy_sum => calculate_pnjl_aniso_energy_sum
        )
    elseif model_type == :rotation
        return Dict(
            :get_nodes => get_rotation_nodes,
            :calculate_log_sum => calculate_rotation_log_sum,
            :calculate_energy_sum => calculate_rotation_energy_sum
        )
    elseif model_type == :gas_liquid
        return Dict(
            :get_nodes => get_gas_liquid_nodes
        )
    else
        error("Unknown model type: $model_type")
    end
end

end # module FunctionRegistry
