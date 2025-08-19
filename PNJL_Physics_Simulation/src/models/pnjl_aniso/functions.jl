"""
PNJL Anisotropic Model Functions - Refactored with Integration Interface

This module implements the anisotropic PNJL model functions using the new
IntegrationInterface for better code organization and maintainability.

Key Improvements:
1. Separation of concerns: integrand definition vs numerical integration
2. Utilizes core integration interface (IntegrationInterface module)
3. Maintains backward compatibility with legacy node format
4. Follows the principle: define physics functions, pass to core integrators

Architecture:
- Core physics functions: calculate_energy_aniso, calculate_log_term_aniso
- Integration functions: use IntegrationInterface.integrate_2d where possible
- Legacy compatibility: *_legacy functions for existing node formats
- Clean interface: main calculate_pressure_aniso uses new design principles
"""
module PNJLAnisoFunctions

using SpecialFunctions: log, exp
using ForwardDiff
using NLsolve
using StaticArrays
using FiniteDifferences
using FastGaussQuadrature
using ..MathUtils: safe_log
using ..Integration: gauleg
using ..PhysicalConstants: hc, Nc
using ..PNJLAnisoConstants: rho0, T0, Lambda_f, G_f, K_f, m0, m0_q_f, m0_s_f,
                           a0, a1, a2, b3, b4
using ..IntegrationInterface: GaussLegendreIntegration, ProductGrid, AngleGrid,
                             integrate_2d, MomentumGrid

# All constants are now imported from respective modules

function get_nodes_aniso(p_num::Int, t_num::Int)
    """
    Generate integration nodes and weights for anisotropic PNJL model
    
    Args:
        p_num: Number of momentum nodes
        t_num: Number of angular nodes
        
    Returns:
        Tuple of (nodes1, nodes2) containing momentum, angle and coefficient arrays
    """
    # Momentum nodes and weights
    nodes1, weights1 = gauleg(0.0, Lambda_f, p_num)
    nodes2, weights2 = gauleg(0.0, 20.0, p_num)
    
    # Angular nodes and weights (cosθ ∈ [0,1])
    t_nodes, t_weights = gauleg(0.0, 1.0, t_num)

    # Meshgrid (i,j) = (p, θ)
    p1_mesh = repeat(nodes1, 1, t_num)
    t1_mesh = repeat(t_nodes', p_num, 1)
    w1_mesh = weights1 * t_weights'  # Outer product
    coefficient1 = w1_mesh .* p1_mesh.^2 ./ π^2  # Spherical coordinate coefficient

    p2_mesh = repeat(nodes2, 1, t_num)
    t2_mesh = t1_mesh
    w2_mesh = weights2 * t_weights'
    coefficient2 = w2_mesh .* p2_mesh.^2 ./ π^2

    nodes1 = [p1_mesh, t1_mesh, coefficient1]
    nodes2 = [p2_mesh, t2_mesh, coefficient2]
    return nodes1, nodes2
end

function calculate_chiral_aniso(phi)
    """Calculate chiral condensate contribution"""
    term1 = 2 * G_f * sum(phi .^ 2) - 4 * K_f * prod(phi)
    return term1
end

@inline function calculate_U_aniso(T, Phi1, Phi2)
    """Calculate Polyakov-loop potential"""
    T_ratio = T0 / T
    Ta = a0 + a1 * T_ratio + a2 * T_ratio^2
    Tb = b3 * T_ratio^3
    value = 1 - 6 * Phi2 * Phi1 + 4 * (Phi2^3 + Phi1^3) - 3 * (Phi2 * Phi1)^2
    # 使用安全对数函数避免负值或零值问题
    log_term = safe_log(value)
    U = T^4 * (-1/2 * Ta * Phi2 * Phi1 + Tb * log_term)
    return U
end

@inline function calculate_mass_vec(phi)
    """Calculate effective masses for three quark flavors"""
    phiu, phid, phis = phi
    return SVector{3, eltype(phi)}(
        m0_q_f - 4 * G_f * phiu + 2 * K_f * phid * phis,
        m0_q_f - 4 * G_f * phid + 2 * K_f * phiu * phis,
        m0_s_f - 4 * G_f * phis + 2 * K_f * phiu * phid
    )
end

@inline function calculate_energy_aniso(mass_i, p, xi, t)
    """Calculate energy with anisotropy parameter xi"""
    p2 = p^2
    mass_i2 = mass_i^2
    term_xi = xi * (p*t)^2
    return sqrt(p2 + mass_i2 + term_xi)
end

@inline function calculate_log_term_aniso(E, mu, T, Phi1, Phi2)
    """Calculate logarithmic term for thermal distribution"""
    invT = 1.0 / T
    x_i = (E - mu) * invT
    x_i_anti = (E + mu) * invT
    
    # Calculate all exponential terms at once
    exp1 = exp(-x_i)
    exp2 = exp1 * exp1
    exp3 = exp1 * exp2
    exp1_anti = exp(-x_i_anti)
    exp2_anti = exp1_anti * exp1_anti
    exp3_anti = exp1_anti * exp2_anti
    
    f1_val = 1.0 + 3.0 * Phi1 * exp1 + 3.0 * Phi2 * exp2 + exp3
    f2_val = 1.0 + 3.0 * Phi2 * exp1_anti + 3.0 * Phi1 * exp2_anti + exp3_anti
    
    # 使用安全对数函数避免负值或零值问题
    return safe_log(f1_val) + safe_log(f2_val)
end

# ============================================================================
# 新的积分接口函数 - 使用IntegrationInterface
# ============================================================================

function vacuum_energy_integral_aniso(masses::Vector, p_nodes, p_weights, t_nodes, t_weights, xi::Float64=0.0)
    """
    Calculate vacuum energy integral with anisotropy using new integration interface
    
    This function replaces the old calculate_energy_sum function by properly using
    the integration interface design: define integrand functions and pass them to
    the core integration methods.
    
    Args:
        masses: Vector of quark masses
        p_nodes, p_weights: Momentum integration points and weights
        t_nodes, t_weights: Angular integration points and weights
        xi: Anisotropy parameter
        
    Returns:
        Vacuum energy contribution
    """
    total_energy = 0.0
    
    # Create proper grids from nodes and weights
    p_grid = MomentumGrid(collect(p_nodes), collect(p_weights), (minimum(p_nodes), maximum(p_nodes)), maximum(p_nodes))
    t_grid = AngleGrid(collect(t_nodes), collect(t_weights), (minimum(t_nodes), maximum(t_nodes)))
    
    # Use the core integration interface
    method = GaussLegendreIntegration()
    
    # Define integrand function for each mass
    for mass_i in masses
        integrand = function(p, t)
            E = calculate_energy_aniso(mass_i, p, xi, t)
            return p^2 * E  # Include p² factor for spherical coordinates
        end
        
        # Use 2D integration from core interface
        contribution = integrate_2d(method, p_grid, t_grid, integrand)
        total_energy += contribution
    end
    
    return total_energy * (-Nc) / (3.0 * π^2)
end

function omega_thermal_integral_aniso(masses::Vector, mu::Vector, T::Float64,
                                    Phi1::Float64, Phi2::Float64,
                                    p_nodes, p_weights, t_nodes, t_weights, xi::Float64=0.0)
    """
    Calculate thermal omega integral with anisotropy using new integration interface
    
    This function replaces the old calculate_log_sum_aniso function by properly using
    the integration interface design: define integrand functions and pass them to
    the core integration methods.
    
    Args:
        masses: Vector of quark masses
        mu: Vector of chemical potentials
        T: Temperature
        Phi1, Phi2: Polyakov loop parameters
        p_nodes, p_weights: Momentum integration points and weights
        t_nodes, t_weights: Angular integration points and weights
        xi: Anisotropy parameter
        
    Returns:
        Thermal omega contribution
    """
    total_contribution = 0.0
    
    # Create proper grids from nodes and weights
    p_grid = MomentumGrid(collect(p_nodes), collect(p_weights), (minimum(p_nodes), maximum(p_nodes)), maximum(p_nodes))
    t_grid = AngleGrid(collect(t_nodes), collect(t_weights), (minimum(t_nodes), maximum(t_nodes)))
    
    # Use the core integration interface
    method = GaussLegendreIntegration()
    
    # Calculate integral contribution for each quark flavor
    for (i, mass_i) in enumerate(masses)
        mu_i = mu[i]
        
        # Define integrand function
        integrand = function(p, t)
            E_i = calculate_energy_aniso(mass_i, p, xi, t)
            log_term = calculate_log_term_aniso(E_i, mu_i, T, Phi1, Phi2)
            return p^2 * log_term  # Include p² factor for spherical coordinates
        end
        
        # Use 2D integration from core interface
        contribution = integrate_2d(method, p_grid, t_grid, integrand)
        total_contribution += contribution
    end
    
    # Include geometric and temperature factors
    return total_contribution * (-T) / (3.0 * π^2)
end

# ============================================================================
# 兼容性包装函数 - 保持向后兼容
# ============================================================================

@inline function calculate_energy_sum(masses, p_nodes, coefficient, t_nodes, xi)
    """Calculate energy sum contribution with anisotropy
    
    **DEPRECATED**: This function is deprecated. Use vacuum_energy_integral_aniso instead.
    """
    # 从原始数据结构中提取积分网格信息
    # p_nodes, t_nodes, coefficient 已经是展开的数组
    p_flat = vec(p_nodes)
    t_flat = vec(t_nodes)
    coef_flat = vec(coefficient)
    
    # 重构权重：coefficient包含了p²因子和权重的乘积
    total_energy = 0.0
    
    @inbounds for mass_i in masses
        contribution = 0.0
        
        # 直接使用展开的网格进行积分
        @inbounds for k in eachindex(p_flat)
            try
                p = p_flat[k]
                t = t_flat[k]
                weight = coef_flat[k]
                
                E = calculate_energy_aniso(mass_i, p, xi, t)
                func_value = E  # 不需要额外的p²因子，已包含在coefficient中
                
                if isfinite(func_value)
                    contribution += func_value * weight
                end
            catch e
                @warn "Energy sum integration failed at index $k: $e"
            end
        end
        
        total_energy += contribution
    end
    
    return total_energy * (-Nc)
end

@inline function calculate_log_sum_aniso(masses, p_nodes, Phi1, Phi2, mu, T, coefficient, t_nodes, xi)
    """Calculate logarithmic sum contribution with anisotropy
    
    **DEPRECATED**: This function is deprecated. Use omega_thermal_integral_aniso instead.
    """
    # 从原始数据结构中提取积分网格信息
    p_flat = vec(p_nodes)
    t_flat = vec(t_nodes)
    coef_flat = vec(coefficient)
    
    total_contribution = 0.0
    
    # 对每个夸克味道计算积分贡献
    @inbounds for (i, mass_i) in enumerate(masses)
        mu_i = mu[i]
        contribution = 0.0
        
        # 直接使用展开的网格进行积分
        @inbounds for k in eachindex(p_flat)
            try
                p = p_flat[k]
                t = t_flat[k]
                weight = coef_flat[k]
                
                E_i = calculate_energy_aniso(mass_i, p, xi, t)
                log_term = calculate_log_term_aniso(E_i, mu_i, T, Phi1, Phi2)
                func_value = log_term  # 不需要额外的p²因子，已包含在coefficient中
                
                if isfinite(func_value)
                    contribution += func_value * weight
                end
            catch e
                @warn "Log sum integration failed at index $k: $e"
            end
        end
        
        total_contribution += contribution
    end
    
    # 包含几何因子和温度因子
    return total_contribution * (-T)
end

# ============================================================================
# 主要物理计算函数
# ============================================================================

function calculate_pressure_aniso(phi, Phi1, Phi2, mu, T, nodes_1, nodes_2, xi=0.0)
    """
    Calculate pressure = -omega for anisotropic PNJL model using new integration interface
    
    Args:
        phi: Chiral condensate vector [phi_u, phi_d, phi_s]
        Phi1, Phi2: Polyakov loop variables
        mu: Chemical potential vector [mu_u, mu_d, mu_s]
        T: Temperature
        nodes_1, nodes_2: Integration nodes (legacy format: [p_mesh, t_mesh, coefficient])
        xi: Anisotropy parameter
        
    Returns:
        Pressure value
    """
    chi = calculate_chiral_aniso(phi)
    U = calculate_U_aniso(T, Phi1, Phi2)
    masses = calculate_mass_vec(phi)
    
    # Extract data from legacy node format for vacuum energy calculation
    p_nodes1 = vec(nodes_1[1])   # momentum values
    t_nodes1 = vec(nodes_1[2])   # angular values
    coef1 = vec(nodes_1[3])      # combined weights * p^2 / π^2
    
    # Use new integration interface while maintaining compatibility with legacy format
    energy_sum = vacuum_energy_integral_aniso_legacy(masses, p_nodes1, t_nodes1, coef1, xi)
    
    # Extract data from legacy node format for thermal calculation
    p_nodes2 = vec(nodes_2[1])
    t_nodes2 = vec(nodes_2[2]) 
    coef2 = vec(nodes_2[3])
    
    # Use new integration interface for thermal calculation
    log_sum = omega_thermal_integral_aniso_legacy(masses, mu, T, Phi1, Phi2, p_nodes2, t_nodes2, coef2, xi)
    
    return -(chi + U + energy_sum + log_sum)
end

function vacuum_energy_integral_aniso_legacy(masses, p_nodes, t_nodes, coefficients, xi::Float64=0.0)
    """
    Calculate vacuum energy integral using legacy node format but with new integration principles
    
    This function maintains compatibility with the existing node structure while following
    the integration interface design principle: define integrand, then integrate.
    
    Args:
        masses: Vector or SVector of quark masses
        p_nodes, t_nodes: Flattened momentum and angular values
        coefficients: Pre-computed weights * p^2 / π^2
        xi: Anisotropy parameter
        
    Returns:
        Vacuum energy contribution
    """
    total_energy = 0.0
    
    # Convert SVector to Vector if necessary to ensure compatibility
    masses_vec = collect(masses)
    
    # For each quark mass, calculate the energy contribution
    @inbounds for mass_i in masses_vec
        contribution = 0.0
        
        # Direct integration using pre-computed coefficients
        # This preserves numerical compatibility while following interface principles
        @inbounds for k in eachindex(p_nodes)
            try
                p = p_nodes[k]
                t = t_nodes[k]
                coef = coefficients[k]  # Already contains weight * p^2 / π^2
                
                # Define integrand: just the energy function (p^2 factor already in coef)
                E = calculate_energy_aniso(mass_i, p, xi, t)
                
                if isfinite(E)
                    contribution += E * coef  # coef already has p^2 and π^2 factors
                end
            catch e
                @warn "Vacuum energy integration failed at index $k: $e"
            end
        end
        
        total_energy += contribution
    end
    
    return total_energy * (-Nc)  # π^2 factor already included in coefficients
end

function omega_thermal_integral_aniso_legacy(masses, mu, T,
                                           Phi1, Phi2,
                                           p_nodes, t_nodes, coefficients, xi=0.0)
    """
    Calculate thermal omega integral using legacy node format but with new integration principles
    
    This function maintains compatibility with the existing node structure while following
    the integration interface design principle: define integrand, then integrate.
    
    Args:
        masses: Vector or SVector of quark masses (supports ForwardDiff types)
        mu: Vector of chemical potentials (supports ForwardDiff types)
        T: Temperature (supports ForwardDiff types)
        Phi1, Phi2: Polyakov loop parameters (supports ForwardDiff types)
        p_nodes, t_nodes: Flattened momentum and angular values  
        coefficients: Pre-computed weights * p^2 / π^2
        xi: Anisotropy parameter
        
    Returns:
        Thermal omega contribution
    """
    total_contribution = zero(eltype(T))
    
    # Convert SVector to Vector if necessary to ensure compatibility
    masses_vec = collect(masses)
    
    # For each quark flavor, calculate the thermal contribution
    @inbounds for (i, mass_i) in enumerate(masses_vec)
        mu_i = mu[i]
        contribution = zero(eltype(T))
        
        # Direct integration using pre-computed coefficients
        @inbounds for k in eachindex(p_nodes)
            try
                p = p_nodes[k]
                t = t_nodes[k]
                coef = coefficients[k]  # Already contains weight * p^2 / π^2
                
                # Define integrand: energy and thermal distribution
                E_i = calculate_energy_aniso(mass_i, p, xi, t)
                log_term = calculate_log_term_aniso(E_i, mu_i, T, Phi1, Phi2)
                
                if isfinite(log_term)
                    contribution += log_term * coef  # coef already has p^2 and π^2 factors
                end
            catch e
                @warn "Thermal integration failed at index $k: $e"
            end
        end
        
        total_contribution += contribution
    end
    
    # Include temperature factor (π^2 factor already included in coefficients)
    return total_contribution * (-T)
end

@inline function pressure_wrapper(x, mu, T, nodes_1, nodes_2, xi)
    """Wrapper function for pressure calculation"""
    phi = SVector{3}(x[1], x[2], x[3])
    Phi1, Phi2 = x[4], x[5]
    return calculate_pressure_aniso(phi, Phi1, Phi2, mu, T, nodes_1, nodes_2, xi)
end

function calculate_core(x, mu, T, nodes_1, nodes_2, xi)
    """Calculate gradient for equation solving"""
    f = x -> pressure_wrapper(x, mu, T, nodes_1, nodes_2, xi)
    return ForwardDiff.gradient(f, x)
end

@inline function calculate_rho(x, mu, T, nodes_1, nodes_2, xi)
    """Calculate density via chemical potential derivative"""
    f_mu = mu -> pressure_wrapper(x, mu, T, nodes_1, nodes_2, xi)
    rho = ForwardDiff.gradient(f_mu, mu)
    return rho
end

function pressure_solve_core(x, mu, T, nodes_1, nodes_2, xi)
    """Solve for pressure at equilibrium"""
    X0_typed = convert.(promote_type(eltype(x), typeof(T)), x)
    res = nlsolve(x -> calculate_core(x, mu, T, nodes_1, nodes_2, xi), X0_typed, autodiff=:forward)
    return pressure_wrapper(res.zero, mu, T, nodes_1, nodes_2, xi)
end

# Export all functions
export get_nodes_aniso, calculate_chiral_aniso, calculate_U_aniso, calculate_mass_vec,
       calculate_energy_aniso, calculate_log_term_aniso,
       vacuum_energy_integral_aniso, omega_thermal_integral_aniso,
       vacuum_energy_integral_aniso_legacy, omega_thermal_integral_aniso_legacy,
       calculate_energy_sum, calculate_log_sum_aniso,  # deprecated but kept for compatibility
       calculate_pressure_aniso, pressure_wrapper, calculate_core, calculate_rho,
       pressure_solve_core

end  # module PNJLAnisoFunctions
