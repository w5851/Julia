include("Constants_PNJL.jl")
include("init3.jl")
using .Constants_PNJL:hc, π, rho0, a0, a1, a2, b3, b4, T0, Nc, Lambda_f, G_f, K_f, m0,m0_q_f,m0_s_f
using .init3:gauleg
using SpecialFunctions: log, exp

# 安全数学函数
@inline function safe_log(x; min_val=1e-16, handle_negative=:clamp)
    """
    安全对数函数，处理负值和接近零的值以确保数值稳定性
    
    参数：
    - x: 输入值
    - min_val: 最小值限制 (默认: 1e-16)
    - handle_negative: 处理负值的方式 (:clamp, :error, :nan)
    
    返回值：
    - 安全的对数值
    """
    if x <= 0
        if handle_negative == :error
            error("safe_log: Input must be positive, got $x")
        elseif handle_negative == :nan
            return NaN
        else  # :clamp
            return log(min_val)
        end
    elseif x < min_val
        return log(min_val)
    else
        return log(x)
    end
end
using ForwardDiff
using NLsolve

using BenchmarkTools
using StaticArrays
using FiniteDifferences

function get_nodes(p_num::Int, t_num::Int)
    # 动量节点和权重
    nodes1, weights1 = gauleg(0.0, Lambda_f, p_num)
    nodes2, weights2 = gauleg(0.0, 20.0, p_num)
    # 角度节点和权重（cosθ ∈ [0,1]）
    t_nodes, t_weights = gauleg(0.0, 1.0, t_num)

    # meshgrid (i,j) = (p, θ)
    p1_mesh = repeat(nodes1, 1, t_num)
    t1_mesh = repeat(t_nodes', p_num, 1)
    w1_mesh = weights1 * t_weights'  # 外积
    coefficient1 = w1_mesh .* p1_mesh.^2 ./ π^2  # 转换为球坐标后相同的系数p^2*4π

    p2_mesh = repeat(nodes2, 1, t_num)
    t2_mesh = t1_mesh
    w2_mesh = weights2 * t_weights'
    coefficient2 = w2_mesh .* p2_mesh.^2 ./ π^2  # 转换为球坐标后相同的系数p^2*4π

    nodes1 = [p1_mesh, t1_mesh, coefficient1]
    nodes2 = [p2_mesh, t2_mesh, coefficient2]
    return nodes1, nodes2
end

@inline function calculate_chiral(phi)
    """计算手征相关量"""
    term1 = 2 * G_f * sum(phi .^ 2) - 4 * K_f * prod(phi)
    return term1
end

@inline function calculate_U(T, Phi1, Phi2)
    """计算极化Polyakov-loop势能"""
    T_ratio = T0 / T
    Ta = a0 + a1 * T_ratio + a2 * T_ratio^2
    Tb = b3 * T_ratio^3
    value = 1 - 6 * Phi2 * Phi1 + 4 * (Phi2^3 + Phi1^3) - 3 * (Phi2 * Phi1)^2
    # 使用安全对数函数避免负值或零值问题
    log_term = safe_log(value)
    U = T^4 * (-1/2 * Ta * Phi2 * Phi1 + Tb * log_term)  # 对数有效势
    return U
end

@inline function calculate_mass_vec(phi)
    """计算三种夸克的有效质量静态向量，兼容自动微分"""
    phiu, phid, phis = phi
    return SVector{3, eltype(phi)}(
        m0_q_f - 4 * G_f * phiu + 2 * K_f * phid * phis,
        m0_q_f - 4 * G_f * phid + 2 * K_f * phiu * phis,
        m0_s_f - 4 * G_f * phis + 2 * K_f * phiu * phid
    )
end

@inline function calculate_energy(mass_i, p, xi, t)
    p2 = p^2
    mass_i2 = mass_i^2
    term_xi = xi * (p*t)^2
    return  sqrt(p2 + mass_i2+ term_xi)
end

@inline function calculate_log_term(E_i, mu_i, T, Phi1, Phi2)
    invT = 1.0 / T  # 预计算倒数
    x_i = (E_i - mu_i) * invT
    x_i_anti = (E_i + mu_i) * invT
    
    # 一次性计算所有指数项
    exp1 =  exp(-x_i)
    exp2 = exp1 * exp1
    exp3 = exp1 * exp2
    exp1_anti =  exp(-x_i_anti)
    exp2_anti = exp1_anti * exp1_anti
    exp3_anti = exp1_anti * exp2_anti
    
    f1_val = 1.0 + 3.0 * Phi1 * exp1 + 3.0 * Phi2 * exp2 + exp3
    f2_val = 1.0 + 3.0 * Phi2 * exp1_anti + 3.0 * Phi1 * exp2_anti + exp3_anti
    
    # 使用安全对数函数避免负值或零值问题
    return safe_log(f1_val) + safe_log(f2_val)
end

@inline function calculate_energy_sum(masses, p_nodes, coefficient, t_nodes, xi)
    """逐元素计算能量和"""
    total = 0.0
    # 完全展开嵌套循环，每次操作单个元素
    @inbounds for i in eachindex(masses)
        mass_i = masses[i]
        
        @inbounds @simd for j in eachindex(p_nodes)
            p = p_nodes[j]
            t = t_nodes[j]
            coefficient_j = coefficient[j]
            E = calculate_energy(mass_i, p, xi, t)
            total += E * coefficient_j
        end
    end
    return total * (-Nc)
end

@inline function calculate_log_sum(masses, p_nodes, Phi1, Phi2, mu, T, coefficient, t_nodes, xi)
    """逐元素计算对数项"""
    total = 0.0
    
    # 完全展开嵌套循环，每次操作单个元素
    @inbounds for i in eachindex(masses)
        mass_i = masses[i]
        mu_i = mu[i]  # 缓存数组元素
        @inbounds @simd for j in eachindex(p_nodes)
            p = p_nodes[j]
            t = t_nodes[j]
            coefficient_j = coefficient[j]
            E_i = calculate_energy(mass_i, p, xi, t)
            # 直接计算单个元素
            log_term = calculate_log_term(E_i, mu_i, T, Phi1, Phi2)
            total += log_term * coefficient_j
        end
    end
    return total * (-T)
end

function calculate_pressure(phi,Phi1,Phi2,mu,T,nodes_1,nodes_2, xi=0.0)
    """计算压力=-omega,T和mu传入前需归一化"""
    # 在函数开始时解包 nodes，并将数组部分转换为视图
    p_nodes1 = @view nodes_1[1][:]  # 假设 nodes[1] 是数组
    t_nodes1 = @view nodes_1[2][:]  # 假设 nodes[1] 是数组
    coef1 = @view nodes_1[3][:]  # 假设 nodes[1] 是数组
    p_nodes2 = @view nodes_2[1][:]  # 假设 nodes[2] 是数组
    t_nodes2 = @view nodes_2[2][:]  # 假设 nodes[2] 是数组
    coef2 = @view nodes_2[3][:]  # 假设 nodes[2] 是数组
    

    chi = calculate_chiral(phi)
    U = calculate_U(T,Phi1,Phi2)
    
    masses = calculate_mass_vec(phi)
    # 计算能量部分
    energy_sum = calculate_energy_sum(masses, p_nodes1, coef1, t_nodes1, xi)
    # 计算 log 部分
    log_sum = calculate_log_sum(masses, p_nodes2, Phi1, Phi2, mu, T, coef2, t_nodes2, xi)
    
    return -(chi+U+energy_sum+log_sum)
end

@inline function pressure_wrapper(x, mu, T, nodes_1, nodes_2, xi)
    phi = SVector{3}(x[1], x[2], x[3])
    Phi1, Phi2 = x[4], x[5]
    return calculate_pressure(phi, Phi1, Phi2, mu, T, nodes_1, nodes_2, xi)
end

function calculate_core(x, mu, T, nodes_1, nodes_2, xi)
    # 创建闭包函数，它捕获T、mu和nodes值
    f = x -> pressure_wrapper(x, mu, T, nodes_1, nodes_2, xi)
    
    return ForwardDiff.gradient(f, x)
end

@inline function calculate_rho(x,mu,T,nodes_1, nodes_2, xi)
    f_mu = mu -> pressure_wrapper(x, mu, T, nodes_1, nodes_2, xi)
    rho = ForwardDiff.gradient(f_mu, mu)
    return rho
end

@inline function calculate_thermo(x , mu,T,nodes_1, nodes_2, xi)
    rho = sum(calculate_rho(x, mu, T, nodes_1, nodes_2, xi)) / (3.0*rho0)

    f_T = T -> pressure_wrapper(x, mu, T, nodes_1, nodes_2, xi)
    entropy = ForwardDiff.derivative(f_T, T)

    pressure = pressure_wrapper(x, mu, T, nodes_1, nodes_2, xi)
    energy = -pressure + sum(mu .* rho) + T * entropy  # 使用热力学关系计算能量

    return pressure,rho, entropy,energy
end

function calculate_t_rho(x,T,rho,nodes_1, nodes_2, xi, fvec=Vector{eltype(x)}(undef, 8))
    x_phi = SVector{5}(x[1:5])
    x_mu = SVector{3}(x[6:8])
    fvec[1:5] .= calculate_core(x_phi, x_mu, T, nodes_1, nodes_2, xi)
    fvec[6] = x_mu[1] - x_mu[2]  # μ_u - μ_d
    fvec[7] = x_mu[2] - x_mu[3]  # μ_d - μ_s
    fvec[8] = sum(calculate_rho(x_phi, x_mu, T, nodes_1, nodes_2, xi)) / (3.0*rho0) - rho
    return fvec
end



function Trho(T_start, T_end)
    # 获取节点（p_num=128, t_num=16）
    nodes_1, nodes_2 = get_nodes(128, 16)

    # 输出目录和文件
    outdir = joinpath(@__DIR__, "..", "output")
    mkpath(outdir)
    outfile = joinpath(outdir, "trho_aniso.csv")

    # 初始x值
    x_initial = [-1.8, -1.8, -2.1, 0.8, 0.8, 320 / hc, 320 / hc, 320 / hc]
    # 保存每个T下rho=3.00的解
    x_rho_3 = copy(x_initial)

    # 打开文件并写入表头，然后逐行追加扫描得到的解
    open(outfile, "w") do io
    println(io, "T,rho,phi_u,phi_d,phi_s,Phi1,Phi2,mu_u,mu_d,mu_s,pressure,entropy,energy,converged")

        # 主循环：按 T 扫描
        for T in T_start:1/hc:T_end
            # 使用上一个T的rho=3.00的解作为初始值
            x = copy(x_rho_3)

            # 首先单独计算rho=3.00的情况
            rho = 3.00
            converged = false
            try
                res = nlsolve(x -> calculate_t_rho(x, T, rho, nodes_1, nodes_2, 0.0), x)
                converged = res.f_converged
                if converged
                    copyto!(x, res.zero)
                    # 保存当前T下rho=3.00的解，供下一个T循环使用
                    copyto!(x_rho_3, x)
                else
                    @warn "Root finding did not converge for T=$T and rho=$rho"
                end
            catch err
                @warn "Exception in root finding for T=$T and rho=$rho: $err"
                converged = false
            end

            # 计算热力学量（若收敛）并写入 csv 行
            if converged
                x_phi = SVector{5}(x[1:5])
                x_mu = SVector{3}(x[6:8])
                pressure, _, entropy, energy = calculate_thermo(x_phi, x_mu, T, nodes_1, nodes_2, 0.0)
            else
                pressure = NaN
                entropy = NaN
                energy = NaN
            end
            println(io, join([T*hc, rho, x..., pressure, entropy, energy, converged], ","))
            flush(io)

            # 然后计算剩余的rho值
            for rho in 2.99:-0.01:0.10
                try
                    res = nlsolve(x -> calculate_t_rho(x, T, rho, nodes_1, nodes_2, 0.0), x)
                    converged = res.f_converged
                    if converged
                        copyto!(x, res.zero)
                    else
                        @warn "Root finding did not converge for T=$T and rho=$rho"
                    end
                catch err
                    @warn "Exception in root finding for T=$T and rho=$rho: $err"
                    converged = false
                end

                if converged
                    x_phi = SVector{5}(x[1:5])
                    x_mu = SVector{3}(x[6:8])
                    pressure, _, entropy, energy = calculate_thermo(x_phi, x_mu, T, nodes_1, nodes_2, 0.0)
                else
                    pressure = NaN
                    entropy = NaN
                    energy = NaN
                end
                println(io, join([T*hc, rho, x..., pressure, entropy, energy, converged], ","))
                flush(io)
            end
        end
    end

    return nothing
end


function pressure_solve_core(x, mu, T, nodes)
    X0_typed = convert.(promote_type(eltype(x), typeof(T)), x)
    res = nlsolve(x -> calculate_core(x,mu,T,nodes), X0_typed, autodiff=:forward)
    return pressure_wrapper(res.zero,mu,T,nodes)
end


function Tmu(;T_start, T_end, T_step, mu_start, mu_end, mu_step)
    # 节点
    nodes_1, nodes_2 = get_nodes(256, 16)

    # 输出文件
    outdir = joinpath(@__DIR__, "..", "output")
    mkpath(outdir)
    outfile = joinpath(outdir, "tmu_aniso.csv")

    # 初始 x（只有 5 个变量：phi_u, phi_d, phi_s, Phi1, Phi2）
    x_initial = [-1.8, -1.8, -2.1, 0.8, 0.8]
    x_prev = copy(x_initial)

    open(outfile, "w") do io
        println(io, "T,mu,phi_u,phi_d,phi_s,Phi1,Phi2,pressure,rho,entropy,energy,converged")

        for T in T_start:T_step:T_end
            # 对每个 T，在 mu 方向扫描
            for mu in mu_start:mu_step:mu_end
                x = copy(x_prev)
                converged = false
                pressure = NaN
                entropy = NaN
                energy = NaN
                rho =NaN
                try
                    # mu_vec 为三夸克相同的化学势
                    mu_vec = SVector{3}(mu, mu, mu)
                    # 使用 calculate_core 求解 5 个未知量
                    res = nlsolve(x -> calculate_core(x, mu_vec, T, nodes_1, nodes_2, 0.0), x; autodiff = :forward)
                    converged = res.f_converged
                    if converged
                        copyto!(x, res.zero)
                        x_prev .= x  # 用当前解作为下一点的初始值
                        # 计算热力学量（pressure, entropy, energy）
                        pressure, rho, entropy, energy = calculate_thermo(x, mu_vec, T, nodes_1, nodes_2, 0.0)
                    else
                        @warn "Root finding did not converge for T=$T and mu=$mu"
                    end
                catch err
                    @warn "Exception in root finding for T=$T and mu=$mu: $err"
                    converged = false
                end

                # 写入 csv：将 T, mu 转换回物理单位（乘 hc）以与 Trho 保持一致
                println(io, join([T*hc, mu*hc, x..., pressure, rho,entropy, energy, converged], ","))
                flush(io)
            end
        end
    end

    return nothing
end

#pressure_solve_core(x, mu, T, nodes)
#@show p = pressure_solve_core(x, mu, T, nodes)
#@code_warntype pressure_solve_core(x, mu, T, nodes)
#res = @benchmark pressure_solve_core(x, mu, T, nodes) samples=100 seconds=10
#display(res)
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
#res = @benchmark dP_dT(x,mu,T,nodes) samples=100 seconds=10
#res = @benchmark dP_dT2(x,mu,T,nodes) samples=100 seconds=10
#res = @benchmark dP_dT3(x,mu,T,nodes) samples=100 seconds=10
#@show dP_dT4(x, mu, T, nodes)
#res = @benchmark dP_dT4(x,mu,T,nodes) samples=100 seconds=10
#display(res)

#fdm = central_fdm(5, 4)  # 五点中心差分，步长自动选择
function dP_dT4_direct(x, mu, T, nodes,fdm)
    f = T -> pressure_solve_core(x, mu, T, nodes)
    return fdm(f, T)
end
#@show dP_dT4_direct(x, mu, T, nodes,fdm)
#res = @benchmark dP_dT4_direct(x, mu, T, nodes,fdm) samples=100 seconds=10
#display(res)

Tmu(T_start=130/hc, T_end=131/hc, T_step=1/hc, mu_start=400/hc, mu_end=0.0, mu_step=-1/hc)
#Trho(T_start=100/hc, T_end=101/hc)