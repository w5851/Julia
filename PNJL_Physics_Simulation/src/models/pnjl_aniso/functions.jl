"""
PNJL Anisotropic Model Functions

This module implements the anisotropic PNJL model functions based on the original
Function_PNJL_aniso.jl with momentum-dependent anisotropy effects.
"""

using SpecialFunctions: log, exp
using ForwardDiff
using NLsolve
using StaticArrays
using FiniteDifferences
using FastGaussQuadrature
using ..MathUtils: safe_log
using ..IntegrationInterface: GaussLegendreIntegration, ProductGrid, AngleGrid,
                             integrate_2d, MomentumGrid

# 直接定义需要的常数和函数
const π = 3.141592653589793
const hc = 197.33

# 从 PNJLAnisoConstants 导入常数
const rho0 = 0.16
const a0 = 6.75
const a1 = -1.95
const a2 = 2.625
const b3 = 0.75
const b4 = 7.5
const T0 = 210 / hc
const Nc = 3.0
const Lambda = 602.3
const G_Lam2 = 1.835
const K_Lam5 = 12.36
const m0_q = 5.5
const m0_s = 140.7
const Lambda_f = Lambda / hc
const G_f = G_Lam2 / Lambda_f^2
const K_f = K_Lam5 / Lambda_f^5
const m0_q_f = m0_q / hc
const m0_s_f = m0_s / hc

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

function calculate_chiral_aniso(T)
    """Calculate chiral condensate contribution"""
    term1 = 2 * G_f * sum(phi .^ 2) - 4 * K_f * prod(phi)
    return term1
end

@inline function calculate_U(T, Phi1, Phi2)
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

@inline function calculate_energy(mass_i, p, xi, t)
    """Calculate energy with anisotropy parameter xi"""
    p2 = p^2
    mass_i2 = mass_i^2
    term_xi = xi * (p*t)^2
    return sqrt(p2 + mass_i2 + term_xi)
end

@inline function calculate_log_term(E_i, mu_i, T, Phi1, Phi2)
    """Calculate logarithmic term for thermal distribution"""
    invT = 1.0 / T
    x_i = (E_i - mu_i) * invT
    x_i_anti = (E_i + mu_i) * invT
    
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

@inline function calculate_energy_sum(masses, p_nodes, coefficient, t_nodes, xi)
    """Calculate energy sum contribution with anisotropy
    
    **DEPRECATED**: This function is deprecated. Use new integration interface instead.
    """
    # 创建二维网格
    p_domain = (0.0, maximum(p_nodes))
    p_grid = MomentumGrid(p_nodes, coefficient, p_domain, p_domain[2])
    t_domain = (-1.0, 1.0)
    t_grid = AngleGrid(t_nodes, ones(length(t_nodes)), t_domain)
    
    method = GaussLegendreIntegration()
    total = 0.0
    
    # 对每个质量计算贡献
    @inbounds for mass_i in masses
        # 定义被积函数
        integrand = function(p, t)
            E = calculate_energy(mass_i, p, xi, t)
            return E
        end
        
        # 执行二维积分
        contribution = integrate_2d(method, p_grid, t_grid, integrand)
        total += contribution
    end
    
    return total * (-Nc)
end

@inline function calculate_log_sum(masses, p_nodes, Phi1, Phi2, mu, T, coefficient, t_nodes, xi)
    """Calculate logarithmic sum contribution with anisotropy
    
    **DEPRECATED**: This function is deprecated. Use new integration interface instead.
    """
    # 创建二维网格
    p_domain = (0.0, maximum(p_nodes))
    p_grid = MomentumGrid(p_nodes, coefficient, p_domain, p_domain[2])
    t_domain = (-1.0, 1.0)
    t_grid = AngleGrid(t_nodes, ones(length(t_nodes)), t_domain)
    
    method = GaussLegendreIntegration()
    total = 0.0
    
    # 对每个夸克味道计算贡献
    @inbounds for (i, mass_i) in enumerate(masses)
        mu_i = mu[i]
        
        # 定义被积函数
        integrand = function(p, t)
            E_i = calculate_energy(mass_i, p, xi, t)
            log_term = calculate_log_term(E_i, mu_i, T, Phi1, Phi2)
            return log_term
        end
        
        # 执行二维积分
        contribution = integrate_2d(method, p_grid, t_grid, integrand)
        total += contribution
    end
    
    return total * (-T)
end

function calculate_pressure(phi, Phi1, Phi2, mu, T, nodes_1, nodes_2, xi=0.0)
    """
    Calculate pressure = -omega for anisotropic PNJL model
    
    Args:
        phi: Chiral condensate vector [phi_u, phi_d, phi_s]
        Phi1, Phi2: Polyakov loop variables
        mu: Chemical potential vector [mu_u, mu_d, mu_s]
        T: Temperature
        nodes_1, nodes_2: Integration nodes
        xi: Anisotropy parameter
        
    Returns:
        Pressure value
    """
    # Unpack nodes
    p_nodes1 = @view nodes_1[1][:]
    t_nodes1 = @view nodes_1[2][:]
    coef1 = @view nodes_1[3][:]
    p_nodes2 = @view nodes_2[1][:]
    t_nodes2 = @view nodes_2[2][:]
    coef2 = @view nodes_2[3][:]

    chi = calculate_chiral(phi)
    U = calculate_U(T, Phi1, Phi2)
    
    masses = calculate_mass_vec(phi)
    energy_sum = calculate_energy_sum(masses, p_nodes1, coef1, t_nodes1, xi)
    log_sum = calculate_log_sum(masses, p_nodes2, Phi1, Phi2, mu, T, coef2, t_nodes2, xi)
    
    return -(chi + U + energy_sum + log_sum)
end

@inline function pressure_wrapper(x, mu, T, nodes_1, nodes_2, xi)
    """Wrapper function for pressure calculation"""
    phi = SVector{3}(x[1], x[2], x[3])
    Phi1, Phi2 = x[4], x[5]
    return calculate_pressure(phi, Phi1, Phi2, mu, T, nodes_1, nodes_2, xi)
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

@inline function calculate_thermo(x, mu, T, nodes_1, nodes_2, xi)
    """Calculate thermodynamic quantities"""
    rho = calculate_rho(x, mu, T, nodes_1, nodes_2, xi)

    f_T = T -> pressure_wrapper(x, mu, T, nodes_1, nodes_2, xi)
    entropy = ForwardDiff.derivative(f_T, T)

    pressure = pressure_wrapper(x, mu, T, nodes_1, nodes_2, xi)
    energy = -pressure + sum(mu .* rho) + T * entropy

    return pressure, rho, entropy, energy
end

function calculate_t_rho(x, T, rho, nodes_1, nodes_2, xi, fvec=Vector{eltype(x)}(undef, 8))
    """Calculate equations for fixed T and rho"""
    x_phi = SVector{5}(x[1:5])
    x_mu = SVector{3}(x[6:8])
    fvec[1:5] .= calculate_core(x_phi, x_mu, T, nodes_1, nodes_2, xi)
    fvec[6] = x_mu[1] - x_mu[2]  # μ_u - μ_d
    fvec[7] = x_mu[2] - x_mu[3]  # μ_d - μ_s
    fvec[8] = sum(calculate_rho(x_phi, x_mu, T, nodes_1, nodes_2, xi)) / (3.0*rho0) - rho
    return fvec
end

function pressure_solve_core(x, mu, T, nodes_1, nodes_2, xi)
    """Solve for pressure at equilibrium"""
    X0_typed = convert.(promote_type(eltype(x), typeof(T)), x)
    res = nlsolve(x -> calculate_core(x, mu, T, nodes_1, nodes_2, xi), X0_typed, autodiff=:forward)
    return pressure_wrapper(res.zero, mu, T, nodes_1, nodes_2, xi)
end

# Temperature derivatives of pressure
function dP_dT(x, mu, T, nodes_1, nodes_2, xi)
    """First derivative of pressure with respect to temperature"""
    f = T -> pressure_solve_core(x, mu, T, nodes_1, nodes_2, xi)
    return ForwardDiff.derivative(f, T)
end

function dP_dT2(x, mu, T, nodes_1, nodes_2, xi)
    """Second derivative of pressure with respect to temperature"""
    f = T -> dP_dT(x, mu, T, nodes_1, nodes_2, xi)
    return ForwardDiff.derivative(f, T)
end

function dP_dT3(x, mu, T, nodes_1, nodes_2, xi)
    """Third derivative of pressure with respect to temperature"""
    f = T -> dP_dT2(x, mu, T, nodes_1, nodes_2, xi)
    return ForwardDiff.derivative(f, T)
end

function dP_dT4(x, mu, T, nodes_1, nodes_2, xi)
    """Fourth derivative of pressure with respect to temperature"""
    f = T -> dP_dT3(x, mu, T, nodes_1, nodes_2, xi)
    return ForwardDiff.derivative(f, T)
end

function Trho_optimized(T_start, T_end, T_step, rho_start, rho_end, rho_step; 
                       nodes_1, nodes_2, xi=0.0, save_results=true, result_file="trho_results.dat")
    """
    Optimized T-rho calculation for anisotropic PNJL model
    
    Args:
        T_start, T_end, T_step: Temperature range and step
        rho_start, rho_end, rho_step: Density range and step
        nodes_1, nodes_2: Integration nodes
        xi: Anisotropy parameter
        save_results: Whether to save results
        result_file: Output file name
        
    Returns:
        Results dictionary or nothing if saving to file
    """
    # Create temperature and density sequences
    T_values = T_start:T_step:T_end
    rho_values = rho_start:rho_step:rho_end
    
    # Result storage
    results = Dict{Float64, Dict{Float64, Vector{Float64}}}()
    
    # Solution caching for better convergence
    prev_solutions = Dict{Float64, Vector{Float64}}()
    curr_solutions = Dict{Float64, Vector{Float64}}()
    
    # Initial solution
    x_initial = [-1.8, -1.8, -2.1, 0.8, 0.8, 320/hc, 320/hc, 320/hc]
    
    for (t_idx, T) in enumerate(T_values)
        empty!(curr_solutions)
        T_results = Dict{Float64, Vector{Float64}}()
        
        println("Processing T = $T ($(t_idx)/$(length(T_values)))")
        
        for rho in rho_values
            # Choose initial value based on previous solutions
            if haskey(prev_solutions, rho)
                x = copy(prev_solutions[rho])
            elseif !isempty(curr_solutions) && rho_step < 0
                closest_rho = findmax(filter(r -> r > rho, keys(curr_solutions)))[1]
                x = copy(curr_solutions[closest_rho])
            elseif !isempty(curr_solutions) && rho_step > 0
                closest_rho = findmin(filter(r -> r < rho, keys(curr_solutions)))[1]
                x = copy(curr_solutions[closest_rho])
            else
                x = copy(x_initial)
            end
            
            # Solve equations
            res = nlsolve(x -> calculate_t_rho(x, T, rho, nodes_1, nodes_2, xi), x)
            
            if res.f_converged
                curr_solutions[rho] = copy(res.zero)
                if save_results
                    T_results[rho] = copy(res.zero)
                end
            else
                @warn "Root finding did not converge for T=$T and rho=$rho"
                # Try with default initial value
                x = copy(x_initial)
                res = nlsolve(x -> calculate_t_rho(x, T, rho, nodes_1, nodes_2, xi), x)
                
                if res.f_converged
                    curr_solutions[rho] = copy(res.zero)
                    if save_results
                        T_results[rho] = copy(res.zero)
                    end
                else
                    @warn "Second attempt failed for T=$T and rho=$rho"
                end
            end
        end
        
        if save_results
            results[T] = T_results
        end
        
        # Update previous solutions
        empty!(prev_solutions)
        for (k, v) in curr_solutions
            prev_solutions[k] = v
        end
    end
    
    if save_results
        # Save results to file
        open(result_file, "w") do io
            println(io, "# T rho phi_u phi_d phi_s Phi1 Phi2 mu_u mu_d mu_s")
            for T in sort(collect(keys(results)))
                for rho in sort(collect(keys(results[T])))
                    sol = results[T][rho]
                    println(io, "$T $rho $(sol[1]) $(sol[2]) $(sol[3]) $(sol[4]) $(sol[5]) $(sol[6]) $(sol[7]) $(sol[8])")
                end
            end
        end
        return nothing
    else
        return results
    end
end
