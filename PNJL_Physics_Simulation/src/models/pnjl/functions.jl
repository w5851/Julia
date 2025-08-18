"""
Functions and calculations for the PNJL (Polyakov-Nambu-Jona-Lasinio) model.

This module contains computational functions for the PNJL model,
including chiral symmetry breaking, Polyakov loop dynamics, and 
thermodynamic calculations.
"""
module PNJLFunctions

using FastGaussQuadrature
using SpecialFunctions: log, exp
using ForwardDiff
using NLsolve
using BenchmarkTools
using StaticArrays
using FiniteDifferences
using ..MathUtils: safe_log
using ..Integration: gauleg  # Import Integration module
using ..IntegrationInterface: GaussLegendreIntegration, MomentumGrid, 
                             omega_thermal_integral, vacuum_energy_integral

# 直接定义需要的常数和函数
const π = 3.141592653589793
const hc = 197.33

# 从 PNJLConstants 导入
const rho0 = 0.16
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
const m0 = [m0_q_f, m0_q_f, m0_s_f]
const a0 = 3.51
const a1 = -2.47
const a2 = 15.2
const b3 = -1.75
const b4 = 7.5

export get_nodes, calculate_chiral, calculate_U, calculate_mass_vec,
       calculate_energy, calculate_log_term, calculate_pressure,
       calculate_core, calculate_rho, calculate_thermo, calculate_t_rho,
       Trho_optimized, pressure_solve_core

"""
    get_nodes(n_p::Int)

Generate integration nodes and weights for momentum integration in PNJL model.

# Arguments
- `n_p::Int`: Number of integration points

# Returns
- Array containing [p_nodes1, p_nodes2, coefficient1, coefficient2]
"""
function get_nodes(n_p::Int)
    p_nodes, p_weights = gauleg(0.0, Lambda_f, n_p)
    p_nodes2, p_weights2 = gauleg(0.0, 20.0, n_p)
    coefficient1 = p_weights .* p_nodes.^2 ./ π^2
    coefficient2 = p_weights2 .* p_nodes2.^2 ./ π^2
    return [p_nodes, p_nodes2, coefficient1, coefficient2]
end

"""
    calculate_chiral(phi)

Calculate chiral condensate contribution to the thermodynamic potential.
"""
@inline function calculate_chiral(phi)
    term1 = 2 * G_f * sum(phi .^ 2) - 4 * K_f * prod(phi)
    return term1
end

"""
    calculate_U(T, Phi1, Phi2)

Calculate Polyakov loop potential.
"""
@inline function calculate_U(T, Phi1, Phi2)
    T_ratio = T0 / T
    Ta = a0 + a1 * T_ratio + a2 * T_ratio^2
    Tb = b3 * T_ratio^3
    value = 1 - 6 * Phi2 * Phi1 + 4 * (Phi2^3 + Phi1^3) - 3 * (Phi2 * Phi1)^2
    # 使用安全对数函数避免负值或零值问题
    log_term = safe_log(value)
    U = T^4 * (-1/2 * Ta * Phi2 * Phi1 + Tb * log_term)
    return U
end

"""
    calculate_mass_vec(phi)

Calculate effective quark masses as a static vector compatible with automatic differentiation.
"""
@inline function calculate_mass_vec(phi)
    phiu, phid, phis = phi
    return SVector{3, eltype(phi)}(
        m0_q_f - 4 * G_f * phiu + 2 * K_f * phid * phis,
        m0_q_f - 4 * G_f * phid + 2 * K_f * phiu * phis,
        m0_s_f - 4 * G_f * phis + 2 * K_f * phiu * phid
    )
end

"""
    calculate_energy(mass_i, p)

Calculate quark energy dispersion relation.
"""
@inline function calculate_energy(mass_i, p)
    p2 = p^2
    mass_i2 = mass_i^2
    return sqrt(p2 + mass_i2)
end

"""
    calculate_log_term(E_i, mu_i, T, Phi1, Phi2)

Calculate logarithmic term in thermodynamic potential with Polyakov loop.
"""
@inline function calculate_log_term(E_i, mu_i, T, Phi1, Phi2)
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

"""
    calculate_energy_sum(masses, p_nodes, coefficient)

Calculate energy sum contribution to thermodynamic potential.

**DEPRECATED**: This function is deprecated. Use `vacuum_energy_integral` from IntegrationInterface instead.
"""
@inline function calculate_energy_sum(masses, p_nodes, coefficient)
    # 创建临时网格以兼容新接口
    domain = (0.0, maximum(p_nodes))
    grid = MomentumGrid(p_nodes, coefficient, domain, domain[2])
    
    # 使用新的积分接口
    result = vacuum_energy_integral(masses, grid, GaussLegendreIntegration())
    
    # 保持与原实现相同的归一化
    return result * (-Nc) / (1/(3.0 * π^2))  # 撤销新接口中的几何因子
end

"""
    calculate_log_sum(masses, p_nodes, Phi1, Phi2, mu, T, coefficient)

Calculate logarithmic sum contribution to thermodynamic potential.

**DEPRECATED**: This function is deprecated. Use `omega_thermal_integral` from IntegrationInterface instead.
"""
@inline function calculate_log_sum(masses, p_nodes, Phi1, Phi2, mu, T, coefficient)
    # 创建临时网格以兼容新接口
    domain = (0.0, maximum(p_nodes))
    grid = MomentumGrid(p_nodes, coefficient, domain, domain[2])
    
    # 使用新的积分接口
    result = omega_thermal_integral(masses, mu, T, Phi1, Phi2, grid, GaussLegendreIntegration())
    
    # 保持与原实现相同的归一化（新接口已包含-T和几何因子）
    return result * (3.0 * π^2)  # 撤销新接口中的几何因子，保持原始行为
end

"""
    calculate_pressure(phi, Phi1, Phi2, mu, T, nodes)

Calculate thermodynamic pressure in PNJL model.
"""
function calculate_pressure(phi, Phi1, Phi2, mu, T, nodes)
    p_nodes1 = @view nodes[1][:]
    p_nodes2 = @view nodes[2][:]
    coef1 = @view nodes[3][:]
    coef2 = @view nodes[4][:]

    chi = calculate_chiral(phi)
    U = calculate_U(T, Phi1, Phi2)
    
    masses = calculate_mass_vec(phi)
    energy_sum = calculate_energy_sum(masses, p_nodes1, coef1)
    log_sum = calculate_log_sum(masses, p_nodes2, Phi1, Phi2, mu, T, coef2)
    
    return -(chi + U + energy_sum + log_sum)
end

"""
    pressure_wrapper(x, mu, T, nodes)

Wrapper function for pressure calculation with vector input.
"""
@inline function pressure_wrapper(x, mu, T, nodes)
    phi = SVector{3}(x[1], x[2], x[3])
    Phi1, Phi2 = x[4], x[5]
    return calculate_pressure(phi, Phi1, Phi2, mu, T, nodes)
end

"""
    calculate_core(x, mu, T, nodes)

Calculate gradient of pressure for field equation solving.
"""
function calculate_core(x, mu, T, nodes)
    f = x -> pressure_wrapper(x, mu, T, nodes)
    return ForwardDiff.gradient(f, x)
end

"""
    calculate_rho(x, mu, T, nodes)

Calculate baryon number density from pressure derivative.
"""
@inline function calculate_rho(x, mu, T, nodes)
    f_mu = mu -> pressure_wrapper(x, mu, T, nodes)
    rho = ForwardDiff.gradient(f_mu, mu)
    return rho
end

"""
    calculate_thermo(x, mu, T, nodes)

Calculate thermodynamic quantities: pressure, density, entropy, energy.
"""
@inline function calculate_thermo(x, mu, T, nodes)
    rho = calculate_rho(x, mu, T, nodes)

    f_T = T -> pressure_wrapper(x, mu, T, nodes)
    entropy = ForwardDiff.derivative(f_T, T)

    pressure = pressure_wrapper(x, mu, T, nodes)
    energy = -pressure + sum(mu .* rho) + T * entropy

    return pressure, rho, entropy, energy
end

"""
    calculate_t_rho(x, T, rho, nodes, fvec)

Calculate residual function for T-ρ phase diagram construction.
"""
function calculate_t_rho(x, T, rho, nodes, fvec=Vector{eltype(x)}(undef, 8))
    x_phi = SVector{5}(x[1:5])
    x_mu = SVector{3}(x[6:8])
    fvec[1:5] .= calculate_core(x_phi, x_mu, T, nodes)
    fvec[6] = x_mu[1] - x_mu[2]  # μ_u - μ_d
    fvec[7] = x_mu[2] - x_mu[3]  # μ_d - μ_s
    fvec[8] = sum(calculate_rho(x_phi, x_mu, T, nodes)) / (3.0*rho0) - rho
    return fvec
end

"""
    Trho_optimized(T_start, T_end, T_step, rho_start, rho_end, rho_step; save_results, result_file)

Optimized calculation of T-ρ phase diagram with smart initial value selection.
"""
function Trho_optimized(T_start, T_end, T_step=1/hc, 
                         rho_start=3.00, rho_end=0.10, rho_step=-0.01; 
                         save_results=true, result_file="trho_results.jld2")
    nodes = get_nodes(128)
    
    T_values = T_start:T_step:T_end
    rho_values = rho_start:rho_step:rho_end
    
    results = Dict{Float64, Dict{Float64, Vector{Float64}}}()
    prev_solutions = Dict{Float64, Vector{Float64}}()
    curr_solutions = Dict{Float64, Vector{Float64}}()
    
    x_initial = [-1.8, -1.8, -2.1, 0.8, 0.8, 320/hc, 320/hc, 320/hc]
    
    for (t_idx, T) in enumerate(T_values)
        empty!(curr_solutions)
        T_results = Dict{Float64, Vector{Float64}}()
        
        println("Processing T = $T ($(t_idx)/$(length(T_values)))")
        
        for rho in rho_values
            # Smart initial value selection
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
            res = nlsolve(x -> calculate_t_rho(x, T, rho, nodes), x)
            
            if res.f_converged
                curr_solutions[rho] = copy(res.zero)
                if save_results
                    T_results[rho] = copy(res.zero)
                end
            else
                @warn "Root finding did not converge for T=$T and rho=$rho"
                # Retry with default initial value
                x = copy(x_initial)
                res = nlsolve(x -> calculate_t_rho(x, T, rho, nodes), x)
                
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
        
        # Slide solutions: current becomes previous
        empty!(prev_solutions)
        for (k, v) in curr_solutions
            prev_solutions[k] = v
        end
    end
    
    return save_results ? nothing : results
end

"""
    pressure_solve_core(x, mu, T, nodes)

Solve for equilibrium pressure by finding field equations' roots.
"""
function pressure_solve_core(x, mu, T, nodes)
    X0_typed = convert.(promote_type(eltype(x), typeof(T)), x)
    res = nlsolve(x -> calculate_core(x, mu, T, nodes), X0_typed, autodiff=:forward)
    return pressure_wrapper(res.zero, mu, T, nodes)
end

# Temperature derivative functions for thermodynamic analysis
function dP_dT(x, mu, T, nodes)
    f = T -> pressure_solve_core(x, mu, T, nodes)
    return ForwardDiff.derivative(f, T)
end

function dP_dT2(x, mu, T, nodes)
    f = T -> dP_dT(x, mu, T, nodes)
    return ForwardDiff.derivative(f, T)
end

function dP_dT3(x, mu, T, nodes)
    f = T -> dP_dT2(x, mu, T, nodes)
    return ForwardDiff.derivative(f, T)
end

function dP_dT4(x, mu, T, nodes)
    f = T -> dP_dT3(x, mu, T, nodes)
    return ForwardDiff.derivative(f, T)
end

"""
    dP_dT4_direct(x, mu, T, nodes, fdm)

Direct calculation of 4th order temperature derivative using finite differences.
"""
function dP_dT4_direct(x, mu, T, nodes, fdm)
    f = T -> pressure_solve_core(x, mu, T, nodes)
    return fdm(f, T)
end

end  # module PNJLFunctions
