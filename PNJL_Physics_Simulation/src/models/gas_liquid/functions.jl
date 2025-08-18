"""
Functions and calculations for the Gas-Liquid phase transition model.

This module contains all computational functions for the Gas-Liquid model,
including thermodynamic calculations, field equations, and phase analysis.
"""
module GasLiquidFunctions

# 使用简化的依赖
using FastGaussQuadrature
using SpecialFunctions: exp, log
using NLsolve
using FiniteDifferences
using BenchmarkTools

# 直接定义需要的常数和函数
const π = 3.141592653589793
const hc = 197.33

# 简化的积分函数
function gauleg(a, b, n)
    t_nodes, t_weights = gausslegendre(n)
    nodes = @. (b - a)/2 * t_nodes + (a + b)/2
    weights = @. (b - a)/2 * t_weights
    return nodes, weights
end

# 热力学函数
@inline function fermion_distribution(E, μ, T)
    return 1 / (exp((E - μ) / T) + 1)
end

@inline function fermion_anti_distribution(E, μ, T)
    return 1 / (exp((E + μ) / T) + 1)
end

# 从 GasLiquidConstants 导入
const m = 939.0 / hc      # Nucleon mass
const mσ = 550.0 / hc     # Sigma meson mass
const mω = 783.0 / hc     # Omega meson mass  
const mρ = 775.0 / hc     # Rho meson mass
const mδ = 980.0 / hc     # Delta meson mass

export get_nodes, calculate_mass, calculate_energy, calculate_ρ, calculate_ρ_s,
       calculate_σ_term, calculate_δ_term, calculate_ρ_term, calculate_ω_term,
       calculate_chemical_constraint, calculate_asymmetry_constraint,
       solve_fun_constraints, calculate_pressure, calculate_pressure_derivatives,
       calculate_thermodynamic_fluctuations, calculate_derivatives_batch

"""
    get_nodes(n_p)

Generate integration nodes and weights for momentum integration.
"""
function get_nodes(n_p)
    p_nodes, p_weights = gauleg(0.0, 20.0, n_p)
    coefficient = @. p_nodes^2 * p_weights / π^2
    return [p_nodes, coefficient]
end

"""
    calculate_mass(gσ, gδ)

Calculate effective masses for protons and neutrons.
"""
@inline function calculate_mass(gσ, gδ)
    m_p = m - gσ - gδ
    m_n = m - gσ + gδ
    return m_p, m_n
end

"""
    calculate_energy(gσ, gδ, p_nodes)

Calculate energy dispersion relations for protons and neutrons.
"""
@inline function calculate_energy(gσ, gδ, p_nodes)
    m_p, m_n = calculate_mass(gσ, gδ)
    p2 = @. p_nodes^2
    E_p = @. sqrt(p2 + m_p^2)
    E_n = @. sqrt(p2 + m_n^2)
    return E_p, E_n
end

"""
    calculate_ρ(E, μ, T, coef)

Calculate baryon number density.
"""
@inline function calculate_ρ(E, μ, T, coef)
    return mapreduce(i -> (fermion_distribution(E[i], μ, T) - 
                          fermion_anti_distribution(E[i], μ, T)) * coef[i], 
                    +, eachindex(E))
end

"""
    calculate_ρ_s(E, μ, T, coef, m)

Calculate scalar density.
"""
@inline function calculate_ρ_s(E, μ, T, coef, m)
    return mapreduce(i -> (fermion_distribution(E[i], μ, T) + 
                          fermion_anti_distribution(E[i], μ, T)) * coef[i] * m / E[i], 
                    +, eachindex(E))
end

"""
    calculate_σ_term(gσ, ρ_ps, ρ_ns, couplings)

Calculate sigma field equation residual.
"""
@inline function calculate_σ_term(gσ, ρ_ps, ρ_ns, couplings)
    fσ, _, _, _, b, c = couplings
    return -gσ + fσ * (ρ_ps + ρ_ns - b * m * gσ^2 - c * gσ^3)
end

"""
    calculate_δ_term(gδ, ρ_ps, ρ_ns, couplings)

Calculate delta field equation residual.
"""
@inline function calculate_δ_term(gδ, ρ_ps, ρ_ns, couplings)
    fδ = couplings[4]
    return -gδ + fδ * (ρ_ps - ρ_ns)
end

"""
    calculate_ρ_term(gρ, ρ_p, ρ_n, couplings)

Calculate rho field equation residual.
"""
@inline function calculate_ρ_term(gρ, ρ_p, ρ_n, couplings)
    fρ = couplings[3]
    return -gρ + fρ * (ρ_p - ρ_n)
end

"""
    calculate_ω_term(gω, ρ_p, ρ_n, couplings)

Calculate omega field equation residual.
"""
@inline function calculate_ω_term(gω, ρ_p, ρ_n, couplings)
    fω = couplings[2]
    return -gω + fω * (ρ_p + ρ_n)
end

"""
    calculate_chemical_constraint(μ_B, μ_n, ρ_p, ρ_n, couplings)

Calculate baryon chemical potential constraint for isospin asymmetric matter.
"""
@inline function calculate_chemical_constraint(μ_B, μ_n, ρ_p, ρ_n, couplings)
    gω = calculate_ω_term(0.0, ρ_p, ρ_n, couplings)
    gρ = calculate_ρ_term(0.0, ρ_p, ρ_n, couplings)
    return μ_B - μ_n - gω + gρ
end

"""
    calculate_asymmetry_constraint(ρ_n, ρ_p, target_asymmetry=0.198)

Calculate isospin asymmetry constraint.
"""
@inline function calculate_asymmetry_constraint(ρ_n, ρ_p, target_asymmetry=0.198)
    return target_asymmetry - (ρ_n - ρ_p)/(ρ_n + ρ_p)
end

"""
    calculate_fun_constraint(x, nodes, couplings, params)

Calculate residual equations for field and chemical potential constraints.
"""
function calculate_fun_constraint(x, nodes, couplings, params)
    gσ, gδ, μ_p, μ_n = x
    p_nodes, coefficient = nodes
    T = params[1]
    μ_B = params[2]
    
    m_p, m_n = calculate_mass(gσ, gδ)
    E_p, E_n = calculate_energy(gσ, gδ, p_nodes)
    ρ_p = calculate_ρ(E_p, μ_p, T, coefficient)
    ρ_n = calculate_ρ(E_n, μ_n, T, coefficient)
    ρ_ps = calculate_ρ_s(E_p, μ_p, T, coefficient, m_p)
    ρ_ns = calculate_ρ_s(E_n, μ_n, T, coefficient, m_n)

    σ_term = calculate_σ_term(gσ, ρ_ps, ρ_ns, couplings)
    δ_term = calculate_δ_term(gδ, ρ_ps, ρ_ns, couplings)
    chem_constraint = calculate_chemical_constraint(μ_B, μ_n, ρ_p, ρ_n, couplings)
    asymmetry_constraint = calculate_asymmetry_constraint(ρ_n, ρ_p)

    return [σ_term, δ_term, chem_constraint, asymmetry_constraint]
end

"""
    solve_fun_constraints(x0, nodes, couplings, params)

Solve field and chemical potential constraint equations.
"""
function solve_fun_constraints(x0, nodes, couplings, params)
    result = nlsolve(x -> calculate_fun_constraint(x, nodes, couplings, params), x0)
    return result.zero
end

"""
    calculate_init_term(E, μ, T, nodes)

Calculate pressure integral term.
"""
@inline function calculate_init_term(E, μ, T, nodes)
    p, coef = nodes
    return mapreduce(i -> (fermion_distribution(E[i], μ, T) + 
                          fermion_anti_distribution(E[i], μ, T)) * coef[i] * p[i]^2 / E[i], 
                    +, eachindex(E)) / 3.0
end

"""
    calculate_pressure(gσ, gδ, gω, gρ, μ_p, μ_n, T, nodes, couplings)

Calculate thermodynamic pressure.
"""
@inline function calculate_pressure(gσ, gδ, gω, gρ, μ_p, μ_n, T, nodes, couplings)
    fσ, fω, fρ, fδ, b, c = couplings
    p_nodes, _ = nodes
    
    E_p, E_n = calculate_energy(gσ, gδ, p_nodes)
    p_p = calculate_init_term(E_p, μ_p, T, nodes)
    p_n = calculate_init_term(E_n, μ_n, T, nodes)
    
    pressure = -(1.0/3.0) * b * m * gσ^3 - 
               (1.0/4.0) * c * gσ^4 - 
               (1.0/(2.0*fσ)) * gσ^2 + 
               (1.0/(2.0*fω)) * gω^2 + 
               p_p + p_n + 
               (1.0/(2.0*fρ)) * gρ^2 - 
               (1.0/(2.0*fδ)) * gδ^2
               
    return pressure
end

"""
    calculate_pressure_wrapper(x, nodes, couplings, params)

Calculate pressure with automatic field determination.
"""
@inline function calculate_pressure_wrapper(x, nodes, couplings, params)
    p_nodes, coef = nodes
    T, _ = params
    gσ, gδ, μ_p, μ_n = x

    E_p, E_n = calculate_energy(gσ, gδ, p_nodes)
    ρ_p = calculate_ρ(E_p, μ_p, T, coef)
    ρ_n = calculate_ρ(E_n, μ_n, T, coef)
    
    gω = calculate_ω_term(0.0, ρ_p, ρ_n, couplings)
    gρ = calculate_ρ_term(0.0, ρ_p, ρ_n, couplings)

    return calculate_pressure(gσ, gδ, gω, gρ, μ_p, μ_n, T, nodes, couplings)
end

"""
    calculate_pressure_solved(μ_B, T, x0, nodes, couplings)

Calculate pressure using solved field values.
"""
function calculate_pressure_solved(μ_B, T, x0, nodes, couplings)
    params = [T, μ_B]
    x = solve_fun_constraints(x0, nodes, couplings, params)
    return calculate_pressure_wrapper(x, nodes, couplings, params)
end

"""
    calculate_pressure_derivatives_efficient(μ_B, T, x0, nodes, couplings; h=1e-5)

Efficiently calculate pressure derivatives up to 4th order with respect to μ_B.
"""
function calculate_pressure_derivatives_efficient(μ_B, T, x0, nodes, couplings; h=1e-5)
    pressure_func = μ -> calculate_pressure_solved(μ, T, x0, nodes, couplings)
    
    fdm1 = central_fdm(5, 1)
    fdm2 = central_fdm(5, 2)
    fdm3 = central_fdm(7, 3)
    fdm4 = central_fdm(7, 4)
    
    pressure = pressure_func(μ_B)
    dpre_dmu1 = fdm1(pressure_func, μ_B)
    dpre_dmu2 = fdm2(pressure_func, μ_B)
    dpre_dmu3 = fdm3(pressure_func, μ_B)
    dpre_dmu4 = fdm4(pressure_func, μ_B)
    
    return pressure, dpre_dmu1, dpre_dmu2, dpre_dmu3, dpre_dmu4
end

"""
    calculate_thermodynamic_fluctuations(μ_B, T, x0, nodes, couplings)

Calculate thermodynamic fluctuation quantities and cumulants.
"""
function calculate_thermodynamic_fluctuations(μ_B, T, x0, nodes, couplings)
    pressure, dpre_dmu1, dpre_dmu2, dpre_dmu3, dpre_dmu4 = 
        calculate_pressure_derivatives_efficient(μ_B, T, x0, nodes, couplings)
    
    kappa1 = dpre_dmu1
    kappa2 = dpre_dmu2
    kappa3 = dpre_dmu3
    kappa4 = dpre_dmu4
    
    fluctuation_ratios = [
        kappa2 / kappa1,
        kappa3 / kappa2,
        kappa4 / kappa2
    ]
    
    return kappa1, kappa2, kappa3, kappa4, fluctuation_ratios
end

"""
    calculate_derivatives_batch(μ_B_array, T, x0, nodes, couplings; save_results=false, output_file)

Batch calculation of pressure derivatives for multiple chemical potential points.
"""
function calculate_derivatives_batch(μ_B_array, T, x0, nodes, couplings; 
                                   save_results=false, 
                                   output_file=joinpath(@__DIR__, "..", "..", "..", "output", "derivatives_output.dat"))
    n_points = length(μ_B_array)
    
    # Pre-allocate result arrays
    pressure_array = zeros(n_points)
    dpre_dmu1_array = zeros(n_points)
    dpre_dmu2_array = zeros(n_points)
    dpre_dmu3_array = zeros(n_points)
    dpre_dmu4_array = zeros(n_points)
    kappa1_array = zeros(n_points)
    kappa2_array = zeros(n_points)
    kappa3_array = zeros(n_points)
    kappa4_array = zeros(n_points)
    fluctuation_ratios_array = zeros(n_points, 3)
    
    println("Starting batch calculation for $(n_points) chemical potential points...")
    
    for (i, μ_B) in enumerate(μ_B_array)
        try
            pressure, dpre_dmu1, dpre_dmu2, dpre_dmu3, dpre_dmu4 = 
                calculate_pressure_derivatives_efficient(μ_B, T, x0, nodes, couplings)
            
            kappa1, kappa2, kappa3, kappa4, fluctuation_ratios = 
                calculate_thermodynamic_fluctuations(μ_B, T, x0, nodes, couplings)
            
            # Store results
            pressure_array[i] = pressure
            dpre_dmu1_array[i] = dpre_dmu1
            dpre_dmu2_array[i] = dpre_dmu2
            dpre_dmu3_array[i] = dpre_dmu3
            dpre_dmu4_array[i] = dpre_dmu4
            kappa1_array[i] = kappa1
            kappa2_array[i] = kappa2
            kappa3_array[i] = kappa3
            kappa4_array[i] = kappa4
            fluctuation_ratios_array[i, :] = fluctuation_ratios
            
            if i % 10 == 0 || i == n_points
                println("Completed: $(i)/$(n_points) ($(round(i/n_points*100, digits=1))%)")
            end
            
        catch e
            println("Warning: Calculation failed at μ_B = $(μ_B*hc) MeV: $e")
            pressure_array[i] = NaN
            dpre_dmu1_array[i] = NaN
            dpre_dmu2_array[i] = NaN
            dpre_dmu3_array[i] = NaN
            dpre_dmu4_array[i] = NaN
            kappa1_array[i] = NaN
            kappa2_array[i] = NaN
            kappa3_array[i] = NaN
            kappa4_array[i] = NaN
            fluctuation_ratios_array[i, :] .= NaN
        end
    end
    
    results = (
        μ_B = μ_B_array,
        T = T,
        pressure = pressure_array,
        dpre_dmu1 = dpre_dmu1_array,
        dpre_dmu2 = dpre_dmu2_array,
        dpre_dmu3 = dpre_dmu3_array,
        dpre_dmu4 = dpre_dmu4_array,
        kappa1 = kappa1_array,
        kappa2 = kappa2_array,
        kappa3 = kappa3_array,
        kappa4 = kappa4_array,
        fluctuation_ratios = fluctuation_ratios_array
    )
    
    if save_results
        save_derivatives_results(results, output_file)
        println("Results saved to file: $output_file")
    end
    
    println("Batch calculation completed!")
    return results
end

"""
    save_derivatives_results(results, filename)

Save derivative calculation results to file.
"""
function save_derivatives_results(results, filename)
    println("Saving results to file: $filename")
    
    output_dir = dirname(filename)
    if !isdir(output_dir)
        println("Creating output directory: $output_dir")
        mkpath(output_dir)
    end
    
    try
        open(filename, "w") do io
            write(io, "# Pressure derivative calculation results\n")
            write(io, "# Columns: μ_B(MeV) T(MeV) P ∂P/∂μ ∂²P/∂μ² ∂³P/∂μ³ ∂⁴P/∂μ⁴ κ₁ κ₂ κ₃ κ₄ κ₂/κ₁ κ₃/κ₂ κ₄/κ₂\n")
            
            for i in eachindex(results.μ_B)
                write(io, "$(results.μ_B[i]*hc) $(results.T*hc) $(results.pressure[i]) ")
                write(io, "$(results.dpre_dmu1[i]) $(results.dpre_dmu2[i]) $(results.dpre_dmu3[i]) $(results.dpre_dmu4[i]) ")
                write(io, "$(results.kappa1[i]) $(results.kappa2[i]) $(results.kappa3[i]) $(results.kappa4[i]) ")
                write(io, "$(results.fluctuation_ratios[i,1]) $(results.fluctuation_ratios[i,2]) $(results.fluctuation_ratios[i,3])\n")
            end
        end
        println("File saved successfully!")
    catch e
        println("Error saving file: $e")
        rethrow(e)
    end
end

end  # module GasLiquidFunctions
