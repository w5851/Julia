"""
Rotation Model Functions

This module implements the PNJL model with rotation effects, including
Bessel function integration and angular momentum quantization.
"""

using SpecialFunctions: log, exp, besselj
using ForwardDiff
using NLsolve
using StaticArrays
using FiniteDifferences
using FastGaussQuadrature
using ..IntegrationInterface: GaussLegendreIntegration, ProductGrid, 
                             integrate_2d, MomentumGrid, AngleGrid,
                             create_angular_momentum_grid
using ..UnifiedConstants: physics_pi, hc, rho0, T0, Nc, Lambda, G_Lam2, K_Lam5, 
                         m0_q, m0_s, m0_q_f, m0_s_f, Lambda_f, G_f, K_f

# Use unified constants (avoid redefinition)
const π = physics_pi
const a0 = 6.75
const a1 = -1.95
const a2 = 2.625
const a3 = -7.44
const b3 = 0.75
const b4 = 7.5
const r0 = 1.2
const C = 0.5

# Rotation-dependent coefficients
const coefficients = Dict{String, Vector{Float64}}(
    "a" => [1.0, 0.1, 0.01],
    "b" => [1.0, 0.05, 0.001],
    "c" => [0.5, 0.02, 0.0001],
    "d" => [0.0, 0.1, 0.01]
)

function init_bessel(p, theta, n, w)
    """
    Initialize Bessel function coefficients for rotation integration
    
    Args:
        p: Momentum values
        theta: Angular values
        n: Quantum numbers
        w: Weight values
        
    Returns:
        Normalized coefficients including Bessel functions
    """
    r0 = RotationConstants.r0
    t = sin.(theta)
    ptr2 = (p .* t .* r0).^2

    # Bessel function terms with array broadcasting
    bessel_term = besselj.(n .+ 1, ptr2) .+ besselj.(n, ptr2)

    coefficient = w .* p.^2 .* t .* bessel_term

    return coefficient ./ (4.0 .* π^2)
end

function get_nodes_rotation(p_num::Int, t_num::Int)
    """
    Generate integration nodes for rotation model
    
    Args:
        p_num: Number of momentum nodes
        t_num: Number of angular nodes
        
    Returns:
        Tuple of (nodes1, nodes2) with 3D meshgrids for p, t, n integration
    """
    # Momentum nodes and weights
    p_nodes, p_weights = gauleg(0.0, Lambda_f, p_num)
    p2_nodes, p2_weights = gauleg(0.0, 20.0, p_num)

    # Angular nodes and weights (polar angle θ ∈ [0, π])
    t_nodes, t_weights = gauleg(0.0, π, t_num)

    # Discrete angular momentum nodes and weights
    n_nodes = collect(-5:1:5)              # -5 to 5, total 11 integers
    n_weights = ones(length(n_nodes))      # All weights equal to 1

    # Create 3D joint meshgrid (p, t, n)
    p_mesh  = reshape(p_nodes,  p_num, 1, 1) .* ones(1, t_num, length(n_nodes))
    t_mesh  = reshape(t_nodes,  1, t_num, 1) .* ones(p_num, 1, length(n_nodes))
    n_mesh  = reshape(n_nodes,  1, 1, length(n_nodes)) .* ones(p_num, t_num, 1)
    w_mesh  = reshape(p_weights, p_num, 1, 1) .* reshape(t_weights, 1, t_num, 1) .* reshape(n_weights, 1, 1, length(n_nodes))
    coefficient1 = init_bessel(p_mesh, t_mesh, n_mesh, w_mesh)

    p2_mesh = reshape(p2_nodes, p_num, 1, 1) .* ones(1, t_num, length(n_nodes))
    t2_mesh = t_mesh
    n2_mesh = n_mesh
    w2_mesh = reshape(p2_weights, p_num, 1, 1) .* reshape(t_weights, 1, t_num, 1) .* reshape(n_weights, 1, 1, length(n_nodes))
    coefficient2 = init_bessel(p2_mesh, t2_mesh, n2_mesh, w2_mesh)

    # Return structure consistent with Python np.stack
    nodes1 = [p_mesh, n_mesh, coefficient1]
    nodes2 = [p2_mesh, n2_mesh, coefficient2]
    return nodes1, nodes2
end

@inline function calculate_chiral_rotation(T)
    """Calculate chiral condensate contribution for rotation model"""
    term1 = G_f * mean(phi.^2)
    return term1
end

@inline function calculate_U_rotation(T, Phi1, Phi2)
    """Calculate standard Polyakov-loop potential (non-rotation dependent)"""
    T_ratio = T0 / T
    T_b2 = a0 + a1 * T_ratio + a2 * T_ratio^2 + a3 * T_ratio^3
    poly = -0.5 * T_b2 * Phi1 * Phi2 - b3/6 * (Phi1^3 + Phi2^3) + b4/4 * (Phi1^2 * Phi2^2)
    return T^4 * poly
end

@inline function calculate_mass_rotation(phi)
    """Calculate effective mass for rotation model (single flavor)"""
    return m0_q_f - 2.0 * G_f * phi
end

@inline function calculate_energy_rotation(mass, p, n, omega)
    """Calculate energy with rotation effects"""
    p2 = p^2
    mass2 = mass^2
    return sqrt(p2 + mass2) - (0.5 + n) * omega
end

@inline function AA(x, T, Phi1, Phi2)
    """Auxiliary function for thermal distribution"""
    exp1 = exp(-x / T)
    exp2 = exp1 * exp1
    exp3 = exp1 * exp2
    f1 = 1.0 + 3.0 * Phi1 * exp1 + 3.0 * Phi2 * exp2 + exp3
    return f1
end

@inline function AAbar(x, T, Phi1, Phi2)
    """Auxiliary function for antiparticle thermal distribution"""
    exp1 = exp(-x / T)
    exp2 = exp1 * exp1
    exp3 = exp1 * exp2
    f2 = 1.0 + 3.0 * Phi2 * exp1 + 3.0 * Phi1 * exp2 + exp3
    return f2
end

@inline function calculate_log_term_rotation(E, mu, T, Phi1, Phi2)
    """Calculate logarithmic contribution for rotation model"""
    f1 = AA(E_i - mu_i, T, Phi1, Phi2)
    f2 = AA(-E_i - mu_i, T, Phi1, Phi2)
    f3 = AAbar(-E_i + mu_i, T, Phi1, Phi2)
    f4 = AAbar(E_i + mu_i, T, Phi1, Phi2)
    return log(f1) + log(f2) + log(f3) + log(f4)
end

@inline function calculate_log_sum_rotation(masses, p_nodes, Phi1, Phi2, mu, T, coefficient, n_nodes, omega)
    """Calculate logarithmic sum for rotation model
    
    **DEPRECATED**: This function is deprecated. Use new integration interface instead.
    """
    # 创建二维网格 (动量 x 角动量量子数)
    p_domain = (0.0, maximum(p_nodes))
    p_grid = MomentumGrid(p_nodes, coefficient, p_domain, p_domain[2])
    
    # 对于角动量量子数，使用离散网格
    n_domain = (minimum(n_nodes), maximum(n_nodes))
    n_grid = AngleGrid(n_nodes, ones(length(n_nodes)), n_domain)
    
    method = GaussLegendreIntegration()
    total = 0.0
    
    # 对每个夸克味道计算贡献
    @inbounds for (i, mass_i) in enumerate(masses)
        mu_i = mu[i]
        
        # 定义被积函数
        integrand = function(p, n)
            E_i = calculate_energy_rotation(mass_i, p, n, omega)
            log_term = calculate_log_term_rotation(E_i, mu_i, T, Phi1, Phi2)
            return log_term
        end
        
        # 执行二维积分
        contribution = integrate_2d(method, p_grid, n_grid, integrand)
        total += contribution
    end
    
    return total * (-T)
end

function calc_factors(T, omega)
    """
    Calculate temperature and rotation dependent factors
    
    Args:
        T: Temperature
        omega: Angular velocity
        
    Returns:
        Tuple of (f, f_inv) factors
    """
    a = coefficients["a"][1] + coefficients["a"][2] * omega^2 + coefficients["a"][3] * omega^4
    b = coefficients["b"][1] + coefficients["b"][2] * omega^2 + coefficients["b"][3] * omega^4
    c_ = coefficients["c"][1] + coefficients["c"][2] * omega^2 + coefficients["c"][3] * omega^4
    d = coefficients["d"][1] + coefficients["d"][2] * omega^2 + coefficients["d"][3] * omega^4
    f = a .* tanh.(b * (T / T0 - c_)) .+ d
    return f, 1 ./ f
end

function calc_U(T, Phi1, Phi2, omega)
    """
    Calculate rotation-dependent Polyakov potential
    
    Args:
        T: Temperature
        Phi1, Phi2: Polyakov loop variables
        omega: Angular velocity
        
    Returns:
        Rotation-dependent potential energy
    """
    f, f_inv = calc_factors(T, omega)
    
    term = -C * f * (T/T0)^2 * Phi1 * Phi2 - (Phi1^3 + Phi2^3)/3 + (1/C) * f_inv * (T0/T)^2 * (Phi1*Phi2)^2
    return T^4 * term
end

function calculate_pressure_rotation(phi, Phi1, Phi2, mu, T, nodes1, omega)
    """
    Calculate pressure for rotation model
    
    Args:
        phi: Chiral condensate (scalar for single flavor)
        Phi1, Phi2: Polyakov loop variables
        mu: Chemical potential
        T: Temperature
        nodes1: Integration nodes
        omega: Angular velocity
        
    Returns:
        Pressure value
    """
    # Unpack nodes
    p_nodes2 = @view nodes1[1][:]
    n_nodes2 = @view nodes1[2][:]
    coef2 = @view nodes1[3][:]

    chi = calculate_chiral(phi)
    U = calculate_U_rotation(T, Phi1, Phi2)
    
    masses = [calculate_mass_rotation(phi)]  # Single flavor as array
    
    # Calculate log part
    log_sum = calculate_log_sum_rotation(masses, p_nodes2, Phi1, Phi2, [mu], T, coef2, n_nodes2, omega)
    
    return -(chi + U + log_sum)
end

@inline function pressure_wrapper(x, mu, T, nodes1, omega)
    """Wrapper function for pressure calculation"""
    phi = x[1]
    Phi1 = x[2]
    Phi2 = x[3]
    return calculate_pressure(phi, Phi1, Phi2, mu, T, nodes1, omega)
end

function calculate_core(x, mu, T, nodes1, omega)
    """Core calculation function for gradient"""
    f = x -> pressure_wrapper(x, mu, T, nodes1, omega)
    return ForwardDiff.gradient(f, x)
end

@inline function calculate_rho(x, mu, T, nodes1, omega)
    """Calculate density for rotation model"""
    f_mu = mu -> pressure_wrapper(x, mu, T, nodes1, omega)
    rho = ForwardDiff.derivative(f_mu, mu)
    return rho
end

@inline function calculate_thermo(x, mu, T, nodes1, omega)
    """Calculate thermodynamic quantities for rotation model"""
    rho = calculate_rho(x, mu, T, nodes1, omega)

    f_T = T -> pressure_wrapper(x, mu, T, nodes1, omega)
    entropy = ForwardDiff.derivative(f_T, T)

    pressure = pressure_wrapper(x, mu, T, nodes1, omega)
    energy = -pressure + mu * rho + T * entropy

    return pressure, rho, entropy, energy
end

# Additional utility functions for rotation model analysis

function calculate_moment_of_inertia(phi, T, omega, nodes1)
    """
    Calculate moment of inertia for rotating quark matter
    
    Args:
        phi: Chiral condensate
        T: Temperature  
        omega: Angular velocity
        nodes1: Integration nodes
        
    Returns:
        Moment of inertia value
    """
    # Implementation would depend on specific rotation model details
    # This is a placeholder for the complete implementation
    return 0.0  # Placeholder
end

function calculate_angular_momentum(phi, T, omega, nodes1)
    """
    Calculate total angular momentum
    
    Args:
        phi: Chiral condensate
        T: Temperature
        omega: Angular velocity  
        nodes1: Integration nodes
        
    Returns:
        Angular momentum value
    """
    # Implementation would depend on specific rotation model details
    # This is a placeholder for the complete implementation
    return 0.0  # Placeholder
end
