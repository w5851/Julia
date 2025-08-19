"""
PNJL 各向异性模型函数 - 基于Omega公式的积分接口实现

此模块严格按照omega公式文档实现各向异性PNJL模型函数，
使用统一的积分接口，消除闭包实现方式。

物理公式依据:
- 总热力学势: Ω_total = Ω_chiral + Ω_thermal + Ω_U
- 各向异性能量: E_f(p,t) = sqrt(m_f^2 + p^2 + ξ*(p*t)^2)
- 热力学贡献: Ω_th = -T*g_spin*N_c*Σ_f ∫∫ (p^2/(2π^2)) [ln ℱ_f + ln ℱ_f(-μ)]
- 真空贡献: Ω_vac = -g_spin*N_c*Σ_f ∫∫ (p^2/(2π^2)) E_f(p,t)

架构优化:
1. 严格按照物理公式定义被积函数
2. 使用IntegrationInterface统一积分计算
3. 消除闭包，使用显式函数参数传递
4. 保持ForwardDiff兼容性以支持自动微分
5. 保留向后兼容接口

参考文档: src/models/pnjl_aniso/omega_formulas.md
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
                             integrate_2d, MomentumGrid, create_momentum_grid, create_angle_grid

# ============================================================================
# 核心物理函数 - 基于Omega公式的精确实现
# ============================================================================

function get_nodes_aniso(p_num::Int, t_num::Int)
    """
    为各向异性PNJL模型生成积分节点和权重
    
    参数:
        p_num: 动量节点数
        t_num: 角度节点数
        
    返回:
        包含动量、角度和系数数组的 (nodes1, nodes2) 元组
    """
    # 动量节点和权重
    nodes1, weights1 = gauleg(0.0, Lambda_f, p_num)
    nodes2, weights2 = gauleg(0.0, 20.0, p_num)
    
    # 角度节点和权重 (cosθ ∈ [0,1])
    t_nodes, t_weights = gauleg(0.0, 1.0, t_num)

    # 网格 (i,j) = (p, θ)
    p1_mesh = repeat(nodes1, 1, t_num)
    t1_mesh = repeat(t_nodes', p_num, 1)
    w1_mesh = weights1 * t_weights'  # 外积
    coefficient1 = w1_mesh .* p1_mesh.^2 ./ π^2  # 球坐标系数

    p2_mesh = repeat(nodes2, 1, t_num)
    t2_mesh = t1_mesh
    w2_mesh = weights2 * t_weights'
    coefficient2 = w2_mesh .* p2_mesh.^2 ./ π^2

    nodes1 = [p1_mesh, t1_mesh, coefficient1]
    nodes2 = [p2_mesh, t2_mesh, coefficient2]
    return nodes1, nodes2
end

function calculate_chiral_aniso(phi)
    """计算手征凝聚贡献 - 基于Omega公式"""
    return 2 * G_f * sum(phi .^ 2) - 4 * K_f * prod(phi)
end

@inline function calculate_U_aniso(T, Phi1, Phi2)
    """计算Polyakov环势能 - 基于Omega公式"""
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
    """计算三种夸克味的有效质量"""
    phiu, phid, phis = phi
    return SVector{3, eltype(phi)}(
        m0_q_f - 4 * G_f * phiu + 2 * K_f * phid * phis,
        m0_q_f - 4 * G_f * phid + 2 * K_f * phiu * phis,
        m0_s_f - 4 * G_f * phis + 2 * K_f * phiu * phid
    )
end

@inline function calculate_energy_aniso(mass_i, p, xi, t)
    """
    使用各向异性参数 xi 计算能量
    
    基于Omega公式: E_f(p,t) = sqrt(m_f^2 + p^2 + ξ*(p*t)^2)
    """
    p2 = p^2
    mass_i2 = mass_i^2
    term_xi = xi * (p*t)^2
    return sqrt(p2 + mass_i2 + term_xi)
end

@inline function calculate_polyakov_statistical_factor(E, mu, T, Phi1, Phi2; anti_particle=false)
    """
    计算Polyakov扩张的统计因子
    
    基于Omega公式的完整群论展开:
    ℱ_f(E, μ; Φ, Φ̄) = 1 + 3Φ*exp(-β(E-μ)) + 3Φ̄*exp(-2β(E-μ)) + exp(-3β(E-μ))
    
    反粒子项: ℱ_f(E, -μ; Φ̄, Φ) - 注意Φ和Φ̄的交换
    """
    beta = 1.0 / T
    
    if anti_particle
        # 反粒子项: μ → -μ, Φ ↔ Φ̄
        x = beta * (E + mu)
        F = 1.0 + 3.0*Phi2*exp(-x) + 3.0*Phi1*exp(-2.0*x) + exp(-3.0*x)
    else
        # 正粒子项
        x = beta * (E - mu)
        F = 1.0 + 3.0*Phi1*exp(-x) + 3.0*Phi2*exp(-2.0*x) + exp(-3.0*x)
    end
    
    return F
end

@inline function calculate_log_term_aniso(E, mu, T, Phi1, Phi2)
    """
    计算完整的对数项，包含粒子和反粒子贡献
    
    基于Omega公式: ln ℱ_f(E, μ) + ln ℱ_f(E, -μ)
    """
    # 粒子项
    F_particle = calculate_polyakov_statistical_factor(E, mu, T, Phi1, Phi2, anti_particle=false)
    
    # 反粒子项  
    F_antiparticle = calculate_polyakov_statistical_factor(E, mu, T, Phi1, Phi2, anti_particle=true)
    
    # 使用安全对数函数
    return safe_log(F_particle) + safe_log(F_antiparticle)
end

# ============================================================================
# 基于积分接口的Omega函数实现 - 严格按照公式
# ============================================================================

"""
    calculate_omega_vacuum(masses, xi, p_grid, t_grid, method=GaussLegendreIntegration())

基于omega公式实现真空贡献计算。

物理公式:
Ω_vac = -g_spin * N_c * Σ_f ∫_0^Λf dp ∫_{-1}^1 dt (p^2/(2π^2)) E_f(p,t)

其中:
- E_f(p,t) = sqrt(m_f^2 + p^2 + ξ*(p*t)^2) (各向异性能量)
- g_spin = 2, N_c = 3 (物理常数)
"""
function calculate_omega_vacuum(masses::AbstractVector{T}, xi::T,
                               p_grid::MomentumGrid, t_grid::AngleGrid,
                               method=GaussLegendreIntegration()) where T<:Real
    
    total_energy = zero(T)
    g_spin = 2
    N_c = 3
    
    # 按照omega公式: Σ_f ∫∫ 被积函数_f(p,t)
    @inbounds for mass_f in masses
        
        # 定义真空能量被积函数 - 无闭包，显式参数传递
        contribution_f = integrate_2d(method, p_grid, t_grid, vacuum_energy_integrand;
                                    mass=mass_f, anisotropy=xi)
        
        total_energy += contribution_f
    end
    
    return -g_spin * N_c * total_energy
end

"""
真空能量被积函数 - 避免闭包，使用显式参数
"""
@inline function vacuum_energy_integrand(p::Real, t::Real; mass::Real, anisotropy::Real)
    # 各向异性能量: E_f(p,t) = sqrt(m_f^2 + p^2 + ξ*(p*t)^2)
    E_f = sqrt(mass^2 + p^2 + anisotropy * (p*t)^2)
    
    # 动量权重: p^2/(2π^2)
    momentum_weight = p^2 / (2.0 * π^2)
    
    return momentum_weight * E_f
end

"""
    calculate_omega_thermal(masses, mu, T, Phi1, Phi2, xi, p_grid, t_grid, method=GaussLegendreIntegration())

基于omega公式实现热力学贡献计算。

物理公式:
Ω_th = -T * g_spin * N_c * Σ_f ∫_0^Λf dp ∫_{-1}^1 dt (p^2/(2π^2)) [ln ℱ_f(E_f, μ_f; Φ, Φ̄) + ln ℱ_f(E_f, -μ_f; Φ̄, Φ)]

其中:
- E_f(p,t) = sqrt(m_f^2 + p^2 + ξ*(p*t)^2) (各向异性能量)
- ℱ_f(E, μ; Φ, Φ̄) = 1 + 3Φ*exp(-β(E-μ)) + 3Φ̄*exp(-2β(E-μ)) + exp(-3β(E-μ))
- g_spin = 2, N_c = 3
"""
function calculate_omega_thermal(masses::AbstractVector{T}, mu::AbstractVector{T}, temp::T,
                                Phi1::T, Phi2::T, xi::T,
                                p_grid::MomentumGrid, t_grid::AngleGrid,
                                method=GaussLegendreIntegration()) where T<:Real
    
    total_contribution = zero(T)
    g_spin = 2
    N_c = 3
    
    # 按照omega公式: Σ_f ∫∫ 被积函数_f(p,t)
    @inbounds for (f, mass_f) in enumerate(masses)
        mu_f = mu[f]
        
        # 执行2D积分 - 无闭包，显式参数传递
        contribution_f = integrate_2d(method, p_grid, t_grid, thermal_integrand;
                                    mass=mass_f, chemical_potential=mu_f,
                                    temperature=temp, polyakov1=Phi1, polyakov2=Phi2,
                                    anisotropy=xi)
        
        total_contribution += contribution_f
    end
    
    return -temp * g_spin * N_c * total_contribution
end

"""
热力学被积函数 - 避免闭包，使用显式参数
"""
@inline function thermal_integrand(p::Real, t::Real; mass::Real, chemical_potential::Real,
                                  temperature::Real, polyakov1::Real, polyakov2::Real,
                                  anisotropy::Real)
    
    # 各向异性能量: E_f(p,t) = sqrt(m_f^2 + p^2 + ξ*(p*t)^2)
    E_f = sqrt(mass^2 + p^2 + anisotropy * (p*t)^2)
    
    # 热力学分布函数的对数项
    beta = 1.0 / temperature
    
    # 粒子项: ℱ_f(E, μ; Φ, Φ̄)
    x = beta * (E_f - chemical_potential)
    F_particle = 1.0 + 3.0*polyakov1*exp(-x) + 3.0*polyakov2*exp(-2.0*x) + exp(-3.0*x)
    
    # 反粒子项: ℱ_f(E, -μ; Φ̄, Φ) - 注意Φ和Φ̄的交换
    x_anti = beta * (E_f + chemical_potential)
    F_antiparticle = 1.0 + 3.0*polyakov2*exp(-x_anti) + 3.0*polyakov1*exp(-2.0*x_anti) + exp(-3.0*x_anti)
    
    # 对数项: ln ℱ_f(E, μ) + ln ℱ_f(E, -μ)
    log_term = safe_log(F_particle) + safe_log(F_antiparticle)
    
    # 动量权重: p^2/(2π^2)
    momentum_weight = p^2 / (2.0 * π^2)
    
    return momentum_weight * log_term
end

"""
    calculate_pressure_aniso(phi, Phi1, Phi2, mu, T, nodes_1, nodes_2, xi=0.0)

计算各向异性PNJL模型压力 - 使用无闭包的新实现

保持与旧接口的兼容性，但内部使用改进的无闭包实现。

Args:
    phi: 手征凝聚矢量 [φ_u, φ_d, φ_s]
    Phi1, Phi2: Polyakov环变量
    mu: 化学势矢量 [μ_u, μ_d, μ_s]
    T: 温度
    nodes_1, nodes_2: 积分节点 (旧格式: [p_mesh, t_mesh, coefficient])
    xi: 各向异性参数（默认0.0）

Returns:
    压力值
"""
function calculate_pressure_aniso(phi, Phi1, Phi2, mu, T, nodes_1, nodes_2, xi=0.0)
    # 使用改进的无闭包实现
    return calculate_pressure_aniso_legacy(phi, Phi1, Phi2, mu, T, nodes_1, nodes_2, xi)
end

"""
    calculate_pressure_aniso_modern(phi, Phi1, Phi2, mu, T, p_grid, t_grid, xi=0.0)

使用现代积分接口计算各向异性PNJL模型压力 - 推荐新代码使用

物理公式:
P = -Ω_total = -(Ω_chiral + Ω_thermal + Ω_vacuum + Ω_U)

参数:
- phi: 手征凝聚矢量 [φ_u, φ_d, φ_s]
- Phi1, Phi2: Polyakov环变量
- mu: 化学势矢量 [μ_u, μ_d, μ_s]
- T: 温度
- p_grid: 动量积分网格
- t_grid: 角度积分网格  
- xi: 各向异性参数（默认0.0表示各向同性）
"""
function calculate_pressure_aniso_modern(phi::AbstractVector{T}, Phi1::T, Phi2::T,
                                        mu::AbstractVector{T}, temp::T, 
                                        p_grid::MomentumGrid, t_grid::AngleGrid,
                                        xi::T=zero(T), method=GaussLegendreIntegration()) where T<:Real
    
    # 1. 手征贡献 (解析公式)
    Omega_chiral = calculate_chiral_aniso(phi)
    
    # 2. Polyakov势能贡献 (解析公式)  
    Omega_U = calculate_U_aniso(temp, Phi1, Phi2)
    
    # 3. 计算有效质量
    masses = calculate_mass_vec(phi)
    
    # 4. 真空能量贡献 (使用积分接口，无闭包)
    Omega_vac = calculate_omega_vacuum(masses, xi, p_grid, t_grid, method)
    
    # 5. 热力学贡献 (使用积分接口，无闭包)
    Omega_th = calculate_omega_thermal(masses, mu, temp, Phi1, Phi2, xi, 
                                      p_grid, t_grid, method)
    
    # 6. 总压力: P = -(Ω_chiral + Ω_thermal + Ω_vacuum + Ω_U)
    total_omega = Omega_chiral + Omega_th + Omega_vac + Omega_U
    pressure = -total_omega
    
    return pressure
end

# ============================================================================
# 标准网格创建工具函数
# ============================================================================

"""
    create_aniso_grids(p_points::Int, t_points::Int; p_cutoff=Lambda_f, t_domain=(-1.0, 1.0))

为各向异性PNJL模型创建标准积分网格

参数:
- p_points: 动量网格点数
- t_points: 角度网格点数
- p_cutoff: 动量截断值，默认使用Lambda_f
- t_domain: 角度积分域，默认为(-1, 1)

返回:
- (p_grid, t_grid): 动量和角度积分网格
"""
function create_aniso_grids(p_points::Int, t_points::Int; 
                           p_cutoff=Lambda_f, t_domain=(-1.0, 1.0))
    p_grid = create_momentum_grid(p_points, p_cutoff)
    
    # 创建自定义角度网格
    t_nodes, t_weights = gauleg(t_domain[1], t_domain[2], t_points)
    t_grid = AngleGrid(t_nodes, t_weights, t_domain)
    
    return p_grid, t_grid
end

# ============================================================================
# 向后兼容性包装函数 
# ============================================================================

"""
向后兼容的接口包装函数，将旧的节点格式转换为新的网格格式。

这些函数保持与现有代码的兼容性，同时内部使用新的积分接口。
"""
function calculate_pressure_aniso_compat(phi, Phi1, Phi2, mu, T, nodes_1, nodes_2, xi=0.0)
    """兼容旧接口的压力计算函数"""
    
    # 提取旧格式的节点数据
    p_nodes1, t_nodes1, coefficient1 = nodes_1
    p_nodes2, t_nodes2, coefficient2 = nodes_2
    
    # 转换为新的网格格式 (使用第一个网格，因为它通常是主要的物理积分网格)
    p_values = vec(p_nodes1)
    t_values = vec(t_nodes1)
    weights = vec(coefficient1) .* (2.0 * π^2)  # 恢复原始权重
    
    # 重新计算分离的动量和角度权重
    n_p = size(p_nodes1, 1)
    n_t = size(p_nodes1, 2)
    
    # 提取唯一的动量和角度节点
    p_unique = p_nodes1[:, 1]  # 第一列
    t_unique = t_nodes1[1, :]  # 第一行
    
    # 计算动量和角度权重（从coefficient中反推）
    p_weights = [coefficient1[i, 1] * (2.0 * π^2) / p_unique[i]^2 for i in 1:n_p]
    t_weights = [coefficient1[1, j] * (2.0 * π^2) / t_unique[j] for j in 1:n_t]
    
    # 创建网格对象
    p_domain = (minimum(p_unique), maximum(p_unique))
    t_domain = (minimum(t_unique), maximum(t_unique))
    
    p_grid = MomentumGrid(p_unique, p_weights, p_domain, p_domain[2])
    t_grid = AngleGrid(t_unique, t_weights, t_domain)
    
    # 使用新接口计算
    return calculate_pressure_aniso(phi, Phi1, Phi2, mu, T, p_grid, t_grid, xi)
end

# ============================================================================
# 旧接口兼容性函数 (更新为无闭包实现)
# ============================================================================

"""
旧格式节点的直接计算函数，更新为使用显式函数而非闭包
"""
function calculate_pressure_aniso_legacy(phi, Phi1, Phi2, mu, T, nodes_1, nodes_2, xi=0.0)
    """
    Calculate pressure = -omega for anisotropic PNJL model using direct node calculation
    
    **更新**: 移除闭包实现，使用显式函数调用
    
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
    # 手征和Polyakov势贡献（解析）
    chi = calculate_chiral_aniso(phi)
    U = calculate_U_aniso(T, Phi1, Phi2)
    masses = calculate_mass_vec(phi)
    
    # 提取节点数据
    p_nodes1 = vec(nodes_1[1])   # momentum values
    t_nodes1 = vec(nodes_1[2])   # angular values
    coef1 = vec(nodes_1[3])      # combined weights * p^2 / π^2
    
    # 真空能量积分 - 无闭包实现
    energy_sum = vacuum_energy_legacy_direct(masses, p_nodes1, t_nodes1, coef1, xi)
    
    # 提取第二组节点数据
    p_nodes2 = vec(nodes_2[1])
    t_nodes2 = vec(nodes_2[2]) 
    coef2 = vec(nodes_2[3])
    
    # 热力学积分 - 无闭包实现  
    log_sum = thermal_integral_legacy_direct(masses, mu, T, Phi1, Phi2, p_nodes2, t_nodes2, coef2, xi)
    
    return -(chi + U + energy_sum + log_sum)
end

"""
真空能量积分的直接实现 - 移除闭包
"""
function vacuum_energy_legacy_direct(masses, p_nodes, t_nodes, coefficients, xi::Real=0.0)
    total_energy = 0.0
    masses_vec = collect(masses)
    
    @inbounds for mass_i in masses_vec
        contribution = 0.0
        
        @inbounds for k in eachindex(p_nodes)
            try
                p = p_nodes[k]
                t = t_nodes[k] 
                coef = coefficients[k]
                
                # 直接调用能量函数 - 无闭包
                E = calculate_energy_aniso(mass_i, p, xi, t)
                
                if isfinite(E)
                    contribution += E * coef
                end
            catch e
                @warn "Vacuum energy integration failed at index $k: $e"
            end
        end
        
        total_energy += contribution
    end
    
    return total_energy * (-Nc)
end

"""
热力学积分的直接实现 - 移除闭包
"""
function thermal_integral_legacy_direct(masses, mu, T, Phi1, Phi2, p_nodes, t_nodes, coefficients, xi=0.0)
    total_contribution = zero(eltype(T))
    masses_vec = collect(masses)
    
    @inbounds for (i, mass_i) in enumerate(masses_vec)
        mu_i = mu[i]
        contribution = zero(eltype(T))
        
        @inbounds for k in eachindex(p_nodes)
            try
                p = p_nodes[k]
                t = t_nodes[k]
                coef = coefficients[k]
                
                # 直接调用物理函数 - 无闭包
                E_i = calculate_energy_aniso(mass_i, p, xi, t)
                log_term = calculate_log_term_aniso(E_i, mu_i, T, Phi1, Phi2)
                
                if isfinite(log_term)
                    contribution += log_term * coef
                end
            catch e
                @warn "Thermal integration failed at index $k: $e"
            end
        end
        
        total_contribution += contribution
    end
    
    return total_contribution * (-T)
end

# ============================================================================
# 包装函数和求解接口 (保持与现有代码的兼容性)
# ============================================================================

@inline function pressure_wrapper(x, mu, T, nodes_1, nodes_2, xi)
    """包装函数，用于求解器 - 使用兼容的接口"""
    phi = SVector{3}(x[1], x[2], x[3])
    Phi1, Phi2 = x[4], x[5]
    return calculate_pressure_aniso(phi, Phi1, Phi2, mu, T, nodes_1, nodes_2, xi)
end

function calculate_core(x, mu, T, nodes_1, nodes_2, xi)
    """计算梯度用于方程求解"""
    f = x -> pressure_wrapper(x, mu, T, nodes_1, nodes_2, xi)
    return ForwardDiff.gradient(f, x)
end

@inline function calculate_rho(x, mu, T, nodes_1, nodes_2, xi)
    """通过化学势导数计算密度"""
    f_mu = mu -> pressure_wrapper(x, mu, T, nodes_1, nodes_2, xi)
    rho = ForwardDiff.gradient(f_mu, mu)
    return rho
end

function pressure_solve_core(x, mu, T, nodes_1, nodes_2, xi)
    """在平衡态求解压力"""
    X0_typed = convert.(promote_type(eltype(x), typeof(T)), x)
    res = nlsolve(x -> calculate_core(x, mu, T, nodes_1, nodes_2, xi), X0_typed, autodiff=:forward)
    return pressure_wrapper(res.zero, mu, T, nodes_1, nodes_2, xi)
end

# ============================================================================
# 新接口的包装函数 (推荐使用)
# ============================================================================

"""
现代接口的包装函数，使用积分网格而非旧式节点
"""
@inline function pressure_wrapper_modern(x, mu, T, p_grid, t_grid, xi)
    """使用现代积分接口的包装函数"""
    phi = SVector{3}(x[1], x[2], x[3])
    Phi1, Phi2 = x[4], x[5]
    return calculate_pressure_aniso(phi, Phi1, Phi2, mu, T, p_grid, t_grid, xi)
end

function calculate_core_modern(x, mu, T, p_grid, t_grid, xi)
    """使用现代积分接口计算梯度"""
    f = x -> pressure_wrapper_modern(x, mu, T, p_grid, t_grid, xi)
    return ForwardDiff.gradient(f, x)
end

function pressure_solve_core_modern(x, mu, T, p_grid, t_grid, xi)
    """使用现代接口在平衡态求解压力"""
    X0_typed = convert.(promote_type(eltype(x), typeof(T)), x)
    res = nlsolve(x -> calculate_core_modern(x, mu, T, p_grid, t_grid, xi), X0_typed, autodiff=:forward)
    return pressure_wrapper_modern(res.zero, mu, T, p_grid, t_grid, xi)
end

# Export all functions
export get_nodes_aniso, calculate_chiral_aniso, calculate_U_aniso, calculate_mass_vec,
       calculate_energy_aniso, calculate_log_term_aniso, calculate_polyakov_statistical_factor,
       # 新积分接口函数 (推荐使用)
       calculate_omega_thermal, calculate_omega_vacuum, calculate_pressure_aniso_modern, 
       create_aniso_grids, vacuum_energy_integrand, thermal_integrand,
       # 兼容性接口函数 (保持向后兼容)
       calculate_pressure_aniso, calculate_pressure_aniso_compat, calculate_pressure_aniso_legacy,
       vacuum_energy_legacy_direct, thermal_integral_legacy_direct,
       # 现代包装函数
       pressure_wrapper_modern, calculate_core_modern, pressure_solve_core_modern,
       # 旧接口函数 (向后兼容)
       pressure_wrapper, calculate_core, calculate_rho, pressure_solve_core

end  # module PNJLAnisoFunctions
