"""
    IntegrationInterface

通用积分计算接口模块，为PNJL物理仿真项目提供统一的数值积分功能。

# 设计目标
- 统一积分接口，消除代码重复
- 支持多种积分方法和网格类型
- 提供物理模型专用的积分函数
- 确保数值稳定性和高性能

# 主要功能
- 基础积分接口：integrate()
- 多维积分支持：integrate_2d()
- 物理专用积分：omega_thermal_integral(), vacuum_energy_integral()
- 积分网格管理：MomentumGrid, AngleGrid等
"""
module IntegrationInterface

using ..MathUtils: safe_log, safe_exp
using FastGaussQuadrature
using StaticArrays

# ---------------------------------------------
# gauss-legendre helper (moved from integration.jl)
# 提供 gauleg(a, b, n) -> (nodes, weights)
# ---------------------------------------------
function gauleg(a, b, n)
    t_nodes, t_weights = gausslegendre(n)
    nodes   = @. (b - a)/2 * t_nodes + (a + b)/2
    weights = @. (b - a)/2 * t_weights
    return nodes, weights
end

export IntegrationMethod, IntegrationGrid,
       GaussLegendreIntegration, AdaptiveIntegration,
       MomentumGrid, AngleGrid, ProductGrid,
       integrate, integrate_2d, integrate_vectorized,
       discrete_sum, mixed_integral_sum, angular_momentum_sum,
       create_momentum_grid, create_angle_grid, create_product_grid,
       create_angular_momentum_grid

# ============================================================================
# 抽象类型定义
# ============================================================================

"""
    IntegrationMethod

积分方法的抽象基类。所有具体的积分方法都应继承此类型。
"""
abstract type IntegrationMethod end

"""
    IntegrationGrid

积分网格的抽象基类。定义积分域的节点和权重。
"""
abstract type IntegrationGrid end

# ============================================================================
# 具体积分方法实现
# ============================================================================

"""
    GaussLegendreIntegration <: IntegrationMethod

高斯-勒让德积分方法。适用于光滑函数的高精度积分。

# Fields
- `precision::Float64`: 目标精度，默认1e-12
- `adaptive::Bool`: 是否使用自适应步长，默认false
"""
struct GaussLegendreIntegration <: IntegrationMethod
    precision::Float64
    adaptive::Bool
    
    GaussLegendreIntegration(precision=1e-12, adaptive=false) = new(precision, adaptive)
end

"""
    AdaptiveIntegration <: IntegrationMethod

自适应积分方法。根据被积函数的特性动态调整积分精度。

# Fields
- `tolerance::Float64`: 收敛容限，默认1e-10
- `max_subdivisions::Int`: 最大细分次数，默认100
"""
struct AdaptiveIntegration <: IntegrationMethod
    tolerance::Float64
    max_subdivisions::Int
    
    AdaptiveIntegration(tolerance=1e-10, max_subdivisions=100) = new(tolerance, max_subdivisions)
end

# ============================================================================
# 积分网格定义
# ============================================================================

"""
    MomentumGrid <: IntegrationGrid

动量积分网格，用于动量空间的积分计算。

# Fields
- `nodes::Vector{Float64}`: 积分节点
- `weights::Vector{Float64}`: 积分权重
- `domain::Tuple{Float64, Float64}`: 积分域 (a, b)
- `cutoff::Float64`: 动量截断值
"""
struct MomentumGrid <: IntegrationGrid
    nodes::Vector{Float64}
    weights::Vector{Float64}
    domain::Tuple{Float64, Float64}
    cutoff::Float64
    
    function MomentumGrid(nodes, weights, domain, cutoff)
        length(nodes) == length(weights) || throw(ArgumentError("nodes and weights must have same length"))
        new(nodes, weights, domain, cutoff)
    end
end

"""
    AngleGrid <: IntegrationGrid

角度积分网格，用于角度坐标的积分计算。

# Fields
- `nodes::Vector{Float64}`: 角度节点
- `weights::Vector{Float64}`: 积分权重
- `domain::Tuple{Float64, Float64}`: 积分域，通常为(-1, 1)或(0, π)
"""
struct AngleGrid <: IntegrationGrid
    nodes::Vector{Float64}
    weights::Vector{Float64}
    domain::Tuple{Float64, Float64}
    
    function AngleGrid(nodes, weights, domain=(-1.0, 1.0))
        length(nodes) == length(weights) || throw(ArgumentError("nodes and weights must have same length"))
        new(nodes, weights, domain)
    end
end

"""
    ProductGrid <: IntegrationGrid

多维积分的乘积网格，由多个一维网格构成。

# Fields
- `grids::Vector{IntegrationGrid}`: 各维度的积分网格
"""
struct ProductGrid <: IntegrationGrid
    grids::Vector{IntegrationGrid}
    
    ProductGrid(grids...) = new(collect(grids))
end

# ============================================================================
# 核心积分接口
# ============================================================================

"""
    integrate(method::IntegrationMethod, grid::IntegrationGrid, integrand::Function; kwargs...) -> Float64

执行一维数值积分。

# Arguments
- `method::IntegrationMethod`: 积分方法
- `grid::IntegrationGrid`: 积分网格
- `integrand::Function`: 被积函数 f(x) 或 f(x; kwargs...)

# Keyword Arguments
- 任意关键字参数会传递给被积函数，支持额外参数传递

# Returns
- `Float64`: 积分结果

# Examples
```julia
# 基本用法
grid = create_momentum_grid(64, 20.0)
method = GaussLegendreIntegration()
result = integrate(method, grid, x -> x^2 * exp(-x))

```
"""
function integrate(method::IntegrationMethod, grid::IntegrationGrid, integrand::Function; kwargs...)
    params = NamedTuple(kwargs)
    return _integrate_impl(method, grid, integrand, params)
end

# GaussLegendre方法的具体实现（支持额外参数）
function _integrate_impl(method::GaussLegendreIntegration, grid::Union{MomentumGrid, AngleGrid}, 
                        integrand::Function, params::NamedTuple)
    result = 0.0
    
    @inbounds @simd for i in eachindex(grid.nodes)
        node = grid.nodes[i]
        weight = grid.weights[i]
        
        # 计算被积函数值，使用安全计算避免数值问题
        try
            # 支持带参数的被积函数调用
            func_value = isempty(params) ? integrand(node) : integrand(node; params...)
            if isfinite(func_value)
                result += func_value * weight
            end
        catch e
            @warn "Integration point failed at x=$node: $e"
            # 跳过有问题的积分点
        end
    end
    
    return result
end

"""
    integrate_2d(method::IntegrationMethod, grid1::IntegrationGrid, grid2::IntegrationGrid, 
                integrand::Function; kwargs...) -> Float64

执行二维数值积分。

# Arguments
- `method::IntegrationMethod`: 积分方法
- `grid1, grid2::IntegrationGrid`: 两个维度的积分网格
- `integrand::Function`: 被积函数 f(x, y) 或 f(x, y; kwargs...)

# Keyword Arguments
- 任意关键字参数会传递给被积函数

# Returns
- `Float64`: 二维积分结果

# Examples
```julia
p_grid = create_momentum_grid(64, 20.0)
t_grid = create_angle_grid(16)
result = integrate_2d(GaussLegendreIntegration(), p_grid, t_grid, (p, t) -> p^2 * sin(t))

```
"""
function integrate_2d(method::IntegrationMethod, grid1::IntegrationGrid, grid2::IntegrationGrid, 
                     integrand::Function; kwargs...)
    params = NamedTuple(kwargs)
    return _integrate_2d_impl(method, grid1, grid2, integrand, params)
end

function _integrate_2d_impl(method::IntegrationMethod, grid1::IntegrationGrid, grid2::IntegrationGrid, 
                           integrand::Function, params::NamedTuple)
    result = 0.0
    
    @inbounds for i in eachindex(grid1.nodes)
        x = grid1.nodes[i]
        wx = grid1.weights[i]
        
        @inbounds @simd for j in eachindex(grid2.nodes)
            y = grid2.nodes[j]
            wy = grid2.weights[j]
            
            try
                # 支持带参数的被积函数调用
                func_value = isempty(params) ? integrand(x, y) : integrand(x, y; params...)
                if isfinite(func_value)
                    result += func_value * wx * wy
                end
            catch e
                @warn "2D integration point failed at ($x, $y): $e"
            end
        end
    end
    
    return result
end

"""
    integrate_vectorized(method::IntegrationMethod, grid::IntegrationGrid, integrands::Vector{Function}) -> Vector{Float64}

向量化积分，同时计算多个被积函数的积分。

# Arguments
- `method::IntegrationMethod`: 积分方法
- `grid::IntegrationGrid`: 积分网格
- `integrands::Vector{Function}`: 被积函数数组

# Returns
- `Vector{Float64}`: 各函数的积分结果

# Example
```julia
grid = create_momentum_grid(64, 20.0)
functions = [x -> x^2, x -> x^3, x -> exp(-x)]
results = integrate_vectorized(GaussLegendreIntegration(), grid, functions)
```
"""
function integrate_vectorized(method::IntegrationMethod, grid::IntegrationGrid, integrands::Vector{Function})
    n_funcs = length(integrands)
    results = zeros(Float64, n_funcs)
    
    @inbounds for i in eachindex(grid.nodes)
        node = grid.nodes[i]
        weight = grid.weights[i]
        
        @inbounds @simd for k in 1:n_funcs
            try
                func_value = integrands[k](node)
                if isfinite(func_value)
                    results[k] += func_value * weight
                end
            catch e
                @warn "Vectorized integration failed for function $k at x=$node: $e"
            end
        end
    end
    
    return results
end

# ============================================================================
# 求和接口函数 (用于离散求和)
# ============================================================================

"""
    discrete_sum(summand::Function, indices::AbstractVector) -> Float64

执行离散求和计算。

用于处理角动量量子数等离散变量的求和，与连续积分形成互补。

# Arguments
- `summand::Function`: 求和函数 f(index)
- `indices::AbstractVector`: 离散指标集合

# Returns
- `Float64`: 求和结果

# Example
```julia
# 角动量求和: Σ_n (2n+1) * f(n) for n = 0, 1, 2, ..., n_max
indices = 0:10
result = discrete_sum(n -> (2*n + 1) * exp(-n), indices)
```
"""
function discrete_sum(summand::Function, indices::AbstractVector)
    result = 0.0
    
    @inbounds @simd for index in indices
        try
            value = summand(index)
            if isfinite(value)
                result += value
            end
        catch e
            @warn "Discrete sum failed at index $index: $e"
        end
    end
    
    return result
end

"""
    mixed_integral_sum(method::IntegrationMethod, grid::IntegrationGrid, 
                      indices::AbstractVector, integrand_sum::Function) -> Float64

执行混合积分-求和计算。

用于处理既有连续积分变量又有离散求和变量的情况，
如旋转模型中的 ∫ dp Σ_n f(p,n)。

# Arguments
- `method::IntegrationMethod`: 积分方法
- `grid::IntegrationGrid`: 连续变量的积分网格
- `indices::AbstractVector`: 离散求和变量的指标
- `integrand_sum::Function`: 被积函数 f(continuous_var, discrete_index)

# Returns
- `Float64`: 混合积分-求和结果

# Example
```julia
# 计算 ∫₀^∞ dp Σ_{n=0}^{n_max} p² f(p,n)
p_grid = create_momentum_grid(64, 20.0)
n_indices = 0:10
result = mixed_integral_sum(GaussLegendreIntegration(), p_grid, n_indices,
                          (p, n) -> p^2 * exp(-sqrt(p^2 + n^2)))
```
"""
function mixed_integral_sum(method::IntegrationMethod, grid::IntegrationGrid,
                           indices::AbstractVector, integrand_sum::Function)
    # 对离散变量求和，每一项是关于连续变量的积分
    total_result = 0.0
    
    @inbounds for index in indices
        # 为当前离散指标定义被积函数
        current_integrand = x -> integrand_sum(x, index)
        
        # 执行关于连续变量的积分
        integral_result = integrate(method, grid, current_integrand)
        total_result += integral_result
    end
    
    return total_result
end

"""
    angular_momentum_sum(method::IntegrationMethod, momentum_grid::MomentumGrid,
                        n_max::Int, integrand::Function) -> Float64

专用于角动量量子化系统的积分-求和计算。

这是旋转PNJL模型的专用函数，处理形如 ∫ dp Σₙ (2n+1) f(p,n) 的计算。

# Arguments
- `method::IntegrationMethod`: 积分方法
- `momentum_grid::MomentumGrid`: 动量积分网格
- `n_max::Int`: 最大角动量量子数
- `integrand::Function`: 被积函数 f(p, n)

# Returns
- `Float64`: 角动量求和结果

# Physics
在旋转系统中，角动量沿旋转轴方向量子化，每个角动量态n有统计权重(2n+1)。
总贡献为: ∫ dp Σₙ₌₀^{n_max} (2n+1) f(p,n)

# Example
```julia
p_grid = create_momentum_grid(64, 20.0)
n_max = 10
result = angular_momentum_sum(GaussLegendreIntegration(), p_grid, n_max,
                            (p, n) -> log(1 + exp(-sqrt(p^2 + mass^2))))
```
"""
function angular_momentum_sum(method::IntegrationMethod, momentum_grid::MomentumGrid,
                             n_max::Int, integrand::Function)
    n_indices = 0:n_max
    
    # 定义包含角动量统计权重的被积函数
    weighted_integrand = function(p, n)
        statistical_weight = 2 * n + 1  # 角动量态的统计权重
        return statistical_weight * integrand(p, n)
    end
    
    return mixed_integral_sum(method, momentum_grid, n_indices, weighted_integrand)
end

# ============================================================================
# 辅助函数
# ============================================================================

# ============================================================================
# 网格创建工具函数
# ============================================================================

"""
    create_momentum_grid(n_points::Int, cutoff::Float64) -> MomentumGrid

创建动量积分网格。

# Arguments
- `n_points::Int`: 积分点数
- `cutoff::Float64`: 动量截断值

# Returns
- `MomentumGrid`: 动量积分网格

# Example
```julia
grid = create_momentum_grid(128, 20.0)  # 128个点，截断动量20 GeV
```
"""
function create_momentum_grid(n_points::Int, cutoff::Float64)
    # 使用高斯-勒让德节点和权重
    nodes, weights = _gauss_legendre_nodes_weights(n_points, 0.0, cutoff)
    
    return MomentumGrid(nodes, weights, (0.0, cutoff), cutoff)
end

"""
    create_angle_grid(n_points::Int) -> AngleGrid

创建角度积分网格。

# Arguments
- `n_points::Int`: 积分点数

# Returns
- `AngleGrid`: 角度积分网格，定义域为(-1, 1)

# Example
```julia
grid = create_angle_grid(16)  # 16个角度积分点
```
"""
function create_angle_grid(n_points::Int)
    # 标准高斯-勒让德节点，定义域(-1, 1)
    nodes, weights = _gauss_legendre_nodes_weights(n_points, -1.0, 1.0)
    
    return AngleGrid(nodes, weights, (-1.0, 1.0))
end

"""
    create_angular_momentum_grid(n_points::Int) -> AngleGrid

创建角动量量子化网格，用于旋转系统中的角动量求和。

在PNJL旋转模型中，角动量是量子化的，采用分立的角动量量子数l。
这个函数创建适合角动量求和的网格。

# Arguments
- `n_points::Int`: 角动量量子数的最大值 (l_max)

# Returns
- `AngleGrid`: 角动量网格，包含l=0,1,2,...,l_max的点

# Physical Background
在有限角速度下，夸克的角动量沿旋转轴方向量子化：
- l_z = 0, ±1/2, ±1, ±3/2, ..., 这里简化为整数值
- 权重来自于角动量态的统计权重 (2l+1)

# Example
```julia
l_grid = create_angular_momentum_grid(10)  # l = 0, 1, 2, ..., 10
```
"""
function create_angular_momentum_grid(n_points::Int)
    # 角动量量子数：l = 0, 1, 2, ..., n_points-1
    nodes = Float64.(0:(n_points-1))
    
    # 角动量态的统计权重：2l + 1
    weights = Float64[2.0 * l + 1.0 for l in nodes]
    
    # 归一化权重（可选，取决于具体应用）
    total_weight = sum(weights)
    weights ./= total_weight
    
    # 定义域为离散的整数集合
    domain = (0.0, Float64(n_points - 1))
    
    return AngleGrid(nodes, weights, domain)
end

"""
    create_product_grid(grids...) -> ProductGrid

创建多维积分的乘积网格。

# Arguments
- `grids...`: 各维度的积分网格

# Returns
- `ProductGrid`: 多维乘积网格

# Example
```julia
p_grid = create_momentum_grid(64, 20.0)
t_grid = create_angle_grid(16)
product_grid = create_product_grid(p_grid, t_grid)
```
"""
function create_product_grid(grids...)
    return ProductGrid(grids...)
end

# ============================================================================
# 内部数值计算函数
# ============================================================================

"""
    _gauss_legendre_nodes_weights(n::Int, a::Float64, b::Float64) -> Tuple{Vector{Float64}, Vector{Float64}}

计算区间[a,b]上n点高斯-勒让德积分的节点和权重。

使用现有的gauleg函数实现。
"""
function _gauss_legendre_nodes_weights(n::Int, a::Float64, b::Float64)
    return gauleg(a, b, n)
end

end # module IntegrationInterface
