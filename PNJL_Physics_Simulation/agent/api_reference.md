# PNJL 物理仿真项目 API 参考文档

> 注：通用流程与维护规范已合并到 `agent/README.md`，请先阅读 README 获取首读顺序与维护规则。

## 核心积分接口 (IntegrationInterface)

### 基础积分函数

#### `integrate(method, grid, integrand; kwargs...)`

**功能**: 执行一维数值积分，支持参数传递

**参数**:

**返回值**:

**使用示例**:
```julia
# 基本用法
grid = create_momentum_grid(64, 20.0)
method = GaussLegendreIntegration()
result = integrate(method, grid, x -> x^2 * exp(-x))

# 带参数用法（推荐）
result = integrate(method, grid, (x; mass, temp) -> x^2 * exp(-sqrt(x^2 + mass^2)/temp);
                  mass=1.5, temp=0.2)
```

**自动微分兼容性**: ✅ 完全支持 ForwardDiff

---

#### `integrate_2d(method, grid1, grid2, integrand; kwargs...)`

**功能**: 执行二维数值积分，支持参数传递

**参数**:
- `method::IntegrationMethod`: 积分方法
- `grid1, grid2::IntegrationGrid`: 两个维度的积分网格
- `integrand::Function`: 被积函数 `f(x, y)` 或 `f(x, y; kwargs...)`
- `kwargs...`: 传递给被积函数的关键字参数

**返回值**:
- `Float64`: 二维积分结果

**使用示例**:
```julia
p_grid = create_momentum_grid(64, 20.0)
t_grid = create_angle_grid(16)

# 各向异性能量积分
result = integrate_2d(method, p_grid, t_grid, 
                     (p, t; mass, xi) -> p^2 * sqrt(p^2 + mass^2 + xi*(p*t)^2);
                     mass=1.0, xi=0.5)
```

**物理应用**: 各向异性模型的核心计算接口

---

### 求和接口函数

#### `discrete_sum(summand, indices)`

**功能**: 执行离散求和计算

**参数**:
- `summand::Function`: 求和函数 `f(index)`
- `indices::AbstractVector`: 离散指标集合

**返回值**:
- `Float64`: 求和结果

**使用示例**:
```julia
# 角动量求和: Σₙ (2n+1) * f(n)
indices = 0:10
result = discrete_sum(n -> (2*n + 1) * exp(-n^2), indices)
```

---

#### `mixed_integral_sum(method, grid, indices, integrand_sum)`

**功能**: 执行混合积分-求和计算

**参数**:
- `method::IntegrationMethod`: 积分方法
- `grid::IntegrationGrid`: 连续变量的积分网格
- `indices::AbstractVector`: 离散求和变量的指标
- `integrand_sum::Function`: 被积函数 `f(continuous_var, discrete_index)`

**返回值**:
- `Float64`: 混合积分-求和结果

**物理应用**: 旋转系统中的 `∫ dp Σₙ f(p,n)` 计算

---

#### `angular_momentum_sum(method, momentum_grid, n_max, integrand)`

**功能**: 专用于角动量量子化系统的积分-求和计算

**参数**:
- `method::IntegrationMethod`: 积分方法
- `momentum_grid::MomentumGrid`: 动量积分网格
- `n_max::Int`: 最大角动量量子数
- `integrand::Function`: 被积函数 `f(p, n)`

**返回值**:
- `Float64`: 角动量求和结果（自动包含统计权重）

**物理背景**: 
在旋转系统中，计算 `∫ dp Σₙ₌₀^{n_max} (2n+1) f(p,n)`，
其中 `(2n+1)` 是角动量态的统计权重。

---

### 网格创建函数

#### `create_momentum_grid(n_points, cutoff)`

**功能**: 创建动量积分网格

**参数**:
- `n_points::Int`: 积分点数
- `cutoff::Float64`: 动量截断值

**返回值**:
- `MomentumGrid`: 动量积分网格

---

#### `create_angle_grid(n_points)`

**功能**: 创建角度积分网格

**参数**:
- `n_points::Int`: 积分点数

**返回值**:
- `AngleGrid`: 角度积分网格，定义域为(-1, 1)

---

## 模型专用积分工具

### PNJL模型 (PNJLIntegrationUtils)

#### `omega_thermal_integral(masses, mu, temperature, Phi1, Phi2, grid)`

**功能**: 计算PNJL模型的热力学Omega积分贡献

**物理公式**: `Ω_thermal = -T ∑ᵢ ∫ dp p² log[ℱ(E)] * Polyakov_factor`

### PNJL各向异性模型 (PNJLAnisoIntegrationUtils)

#### `omega_thermal_integral_aniso(masses, mu, T, Phi1, Phi2, p_grid, t_grid, xi)`

**功能**: 计算各向异性PNJL模型的热力学积分

**物理公式**: 包含各向异性能量 `E(p,t) = √(p² + m² + ξ(pt)²)`

---

## 当前接口函数参考（保持兼容性）

### 核心计算函数

#### `calculate_pressure(phi, Phi1, Phi2, mu, T, nodes_1, nodes_2, xi=0.0)`

**功能**: 计算系统压力（负的大正则势）

**参数**:
- `phi::Union{Vector,SVector}`: 手征凝聚态参数 [φᵤ, φᵈ, φₛ]
- `Phi1::Real`: Polyakov Loop 参数 1  
- `Phi2::Real`: Polyakov Loop 参数 2
- `mu::Union{Vector,SVector}`: 化学势 [μᵤ, μᵈ, μₛ] (已归一化)
- `T::Real`: 温度 (已归一化)
- `nodes_1::Tuple`: 第一个积分网格 (动量范围 [0, Λf])
- `nodes_2::Tuple`: 第二个积分网格 (动量范围 [0, 20])  
- `xi::Real`: 各向异性参数 (默认 0.0)

**返回值**:
- `Float64`: 系统压力值

**物理背景**: 计算 PNJL 模型的热力学压力，包括手征项、Polyakov Loop 势能项、费米子真空贡献和有限温度/密度贡献。

**数学表达式**:
```
P = -(χ + U + E_vacuum + Ω_thermal)
```

**使用示例**:
```julia
phi = [-0.1, -0.1, -1.7]
Phi1, Phi2 = 0.5, 0.5
mu = [320/hc, 320/hc, 320/hc]
T = 150/hc
nodes_1, nodes_2 = get_nodes(128, 16)
pressure = calculate_pressure(phi, Phi1, Phi2, mu, T, nodes_1, nodes_2, 0.4)
```

**注意事项**:
- 所有输入参数应已按 `hc` 归一化
- `xi` 参数控制动量各向异性强度
- 函数内部使用数值积分，计算开销较大

---

#### `calculate_core(x, mu, T, nodes_1, nodes_2, xi)`

**功能**: 计算压力对所有场变量的梯度（用于寻找临界点）

**参数**:
- `x::Vector{5}`: 场变量 [φᵤ, φᵈ, φₛ, Φ₁, Φ₂]
- `mu::Vector{3}`: 化学势向量
- `T::Real`: 温度
- `nodes_1, nodes_2::Tuple`: 积分网格
- `xi::Real`: 各向异性参数

**返回值**:
- `Vector{5}`: 梯度向量 [∂P/∂φᵤ, ∂P/∂φᵈ, ∂P/∂φₛ, ∂P/∂Φ₁, ∂P/∂Φ₂]

**物理背景**: 在热力学平衡态，压力对所有场变量的偏导数应为零。这些梯度用于 gap equation 的求解。

**自动微分**: 使用 ForwardDiff 自动计算梯度，确保数值精度。

**使用示例**:
```julia
x = [-0.1, -0.1, -1.7, 0.5, 0.5]
gradient = calculate_core(x, mu, T, nodes_1, nodes_2, xi)
# 平衡态: gradient ≈ [0, 0, 0, 0, 0]
```

---

#### `calculate_thermo(x, mu, T, nodes_1, nodes_2, xi)`

**功能**: 计算完整的热力学量

**参数**: 与 `calculate_core` 相同

**返回值**:
- `pressure::Float64`: 压力
- `rho::Vector{3}`: 数密度 [nᵤ, nᵈ, nₛ]  
- `entropy::Float64`: 熵密度
- `energy::Float64`: 能量密度

**热力学关系**:
```
ρᵢ = ∂P/∂μᵢ
s = ∂P/∂T  
ε = -P + Σᵢ μᵢρᵢ + Ts
```

**使用示例**:
```julia
P, rho, s, epsilon = calculate_thermo(x, mu, T, nodes_1, nodes_2, xi)
println("Pressure: $P, Baryon density: $(sum(rho)/(3*rho0))")
```

---

#### `calculate_t_rho(x, T, rho, nodes_1, nodes_2, xi, fvec)`

**功能**: 计算固定温度和重子数密度下的约束方程组

**参数**:
- `x::Vector{8}`: 扩展场变量 [φᵤ, φᵈ, φₛ, Φ₁, Φ₂, μᵤ, μᵈ, μₛ]
- `T::Real`: 固定温度
- `rho::Real`: 目标重子数密度 (相对于核物质密度 ρ₀)
- `nodes_1, nodes_2::Tuple`: 积分网格
- `xi::Real`: 各向异性参数
- `fvec::Vector{8}`: 输出向量（原地修改）

**约束条件**:
1. Gap equations: ∇P = 0 (前5个方程)
2. Isospin symmetry: μᵤ = μᵈ (第6个方程)  
3. Strangeness neutrality: μᵈ = μₛ (第7个方程)
4. Baryon density constraint: Σᵢnᵢ/(3ρ₀) = rho (第8个方程)

**使用示例**:
```julia
x = [-1.8, -1.8, -2.1, 0.8, 0.8, 320/hc, 320/hc, 320/hc]
fvec = zeros(8)
calculate_t_rho(x, 150/hc, 1.0, nodes_1, nodes_2, 0.0, fvec)
# 解: fvec ≈ [0, 0, 0, 0, 0, 0, 0, 0]
```

---

### 辅助函数

#### `get_nodes(p_num::Int, t_num::Int)`

**功能**: 生成 Gauss-Legendre 积分节点和权重

**参数**:
- `p_num::Int`: 动量方向节点数
- `t_num::Int`: 角度方向节点数

**返回值**:
- `nodes_1::Tuple`: 低能积分网格 (动量 ∈ [0, Λf])
- `nodes_2::Tuple`: 高能积分网格 (动量 ∈ [0, 20])

**网格结构**: 每个 nodes 包含 [p_mesh, t_mesh, coefficient_mesh]

---

#### `calculate_chiral(phi)`

**功能**: 计算手征相互作用项

**数学表达式**:
```
χ = 2Gf Σᵢφᵢ² - 4Kf φᵤφᵈφₛ
```

**使用**:
```julia
phi = [-0.1, -0.1, -1.7]
chiral_term = calculate_chiral(phi)
```

---

#### `calculate_U(T, Phi1, Phi2)`

**功能**: 计算 Polyakov Loop 势能

**⚠️ 已知问题**: 当 `value = 1 - 6Φ₂Φ₁ + 4(Φ₂³ + Φ₁³) - 3(Φ₂Φ₁)²` 为负数时，`log(value)` 会导致数值错误。

**需要修复**: 实现安全对数函数处理负值情况。

**数学表达式**:
```
U = T⁴[-½Tₐ Φ₂Φ₁ + Tᵦ ln(1 - 6Φ₂Φ₁ + 4(Φ₂³ + Φ₁³) - 3(Φ₂Φ₁)²)]
```

---

#### `calculate_mass_vec(phi)`

**功能**: 计算夸克的组成质量

**返回**: `SVector{3}` 包含 [mᵤ, mᵈ, mₛ]

**数学表达式**:
```
mᵤ = m₀ᵤ - 4Gfφᵤ + 2Kfφᵈφₛ
mᵈ = m₀ᵈ - 4Gfφᵈ + 2Kfφᵤφₛ  
mₛ = m₀ₛ - 4Gfφₛ + 2Kfφᵤφᵈ
```

---

### 高级应用函数

#### `Trho_optimized(T_start, T_end, T_step, rho_start, rho_end, rho_step)`

**功能**: 执行温度-密度扫描计算相图

**优化特性**:
- 使用滑动窗口算法提高收敛性
- 智能初值选择减少迭代次数
- 支持结果保存和加载

**参数**:
- `T_start, T_end, T_step`: 温度扫描范围
- `rho_start, rho_end, rho_step`: 密度扫描范围
- `save_results::Bool`: 是否保存结果
- `result_file::String`: 结果文件名

**算法策略**:
1. 滑动窗口：使用相邻点的解作为初值
2. 多层回退：收敛失败时尝试不同初值
3. 进度跟踪：显示计算进度

---

### 微分计算函数

#### `dP_dT(x, mu, T, nodes)`, `dP_dT2(x, mu, T, nodes)`, ...

**功能**: 计算压力对温度的高阶导数

**方法**: 
- 嵌套自动微分 (ForwardDiff)
- 有限差分 (FiniteDifferences.jl)

**性能比较**: 有限差分在高阶导数计算中通常更快

**使用示例**:
```julia
# 自动微分方法
d4P_dT4_AD = dP_dT4(x, mu, T, nodes)

# 有限差分方法  
fdm = central_fdm(5, 4)
d4P_dT4_FD = dP_dT4_direct(x, mu, T, nodes, fdm)
```

---

## 待实现的接口函数

### 1. 安全数学函数

```julia
"""
安全对数函数，处理负值和零值输入

# 参数
- `x::Real`: 输入值
- `min_val::Real`: 最小正值阈值 (默认 1e-16)
- `handle_negative::Symbol`: 负值处理策略 (:clamp, :complex, :error)

# 返回值  
- `Real`: 对数值或处理后的替代值

# 使用示例
```julia
safe_value = safe_log(-0.1)  # 返回 log(1e-16) 而不是报错
```
"""
function safe_log(x; min_val=1e-16, handle_negative=:clamp)
    # 待实现
end
```

### 2. 抽象模型接口

```julia
abstract type PhysicsModel end

"""
通用计算接口

适用于所有物理模型的统一计算接口

# 参数
- `model::PhysicsModel`: 物理模型实例
- `quantity::Symbol`: 要计算的物理量 (:pressure, :entropy, :energy_density)
- `state`: 系统状态参数
- `config`: 计算配置

# 返回值
根据 quantity 类型返回对应的物理量数值

# 扩展性
新模型通过实现此接口即可集成到框架中
"""
function calculate(model::PhysicsModel, quantity::Symbol, state, config)
    error("Must implement calculate for $(typeof(model))")
end
```

### 3. 配置管理接口

```julia
"""
创建计算配置

# 参数
- `grid_points::Tuple{Int,Int}`: 网格点数 (动量, 角度)
- `momentum_cutoff::Float64`: 动量截断
- `numerical_tolerance::Float64`: 数值容差
- `solver_method::Symbol`: 求解器类型

# 返回值
- `ComputationConfig`: 标准化配置对象
"""
function create_config(; grid_points=(128,16), momentum_cutoff=20.0, 
                       numerical_tolerance=1e-12, solver_method=:nlsolve)
    # 待实现
end
```

---

## 使用模式和最佳实践

### 1. 典型计算流程

```julia
# 1. 初始化
nodes_1, nodes_2 = get_nodes(128, 16)
x = [-0.1, -0.1, -1.7, 0.5, 0.5]
mu = [320/hc, 320/hc, 320/hc]
T = 150/hc

# 2. 寻找平衡态
f = x -> calculate_core(x, mu, T, nodes_1, nodes_2, 0.0)
result = nlsolve(f, x)

if result.f_converged
    equilibrium_state = result.zero
    
    # 3. 计算热力学量
    P, rho, s, ε = calculate_thermo(equilibrium_state, mu, T, nodes_1, nodes_2, 0.0)
    
    println("平衡态结果:")
    println("压力: $P")
    println("重子数密度: $(sum(rho)/(3*rho0))")
end
```

### 2. 相图计算模式

```julia
# 温度-密度扫描
Trho_optimized(100/hc, 200/hc, 1/hc,     # 温度范围
               3.0, 0.1, -0.01,           # 密度范围
               save_results=true,         # 保存结果
               result_file="phase_diagram.jld2")
```

### 3. 性能优化建议

- 预分配数组避免内存分配开销
- 使用 `@inbounds` 和 `@simd` 优化循环
- 合理选择网格点数平衡精度和速度
- 缓存重复计算的结果

---

**文档版本**: v1.0  
**最后更新**: 2025年8月18日  
**维护者**: 项目开发团队

**注意**: 此文档反映当前代码状态。随着重构进行，接口将逐步标准化和改进。
