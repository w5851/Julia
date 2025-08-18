# PNJL 物理仿真项目架构设计

## 总体架构原则

### 设计理念
- **模块化**: 功能相关的代码组织在一起
- **抽象化**: 使用接口隐藏实现细节
- **可扩展**: 新功能可以通过实现接口轻松添加
- **可测试**: 每个模块可以独立测试
- **高性能**: 保持 Julia 的计算性能优势

### 依赖层次结构
```
应用层 (Application Layer)
    ↓
物理模型层 (Physics Models Layer)
    ↓
核心计算层 (Core Computation Layer)
    ↓
基础工具层 (Foundation Layer)
```

## 模块架构设计

### 1. Foundation Layer (基础工具层)

#### 1.1 数学工具模块 (`src/core/math_utils.jl`)
```julia
module MathUtils

# 安全数学函数
function safe_log(x; min_val=1e-16, handle_negative=:clamp)
function safe_exp(x; max_val=700.0)
function safe_sqrt(x; min_val=0.0)
function safe_division(numerator, denominator; epsilon=1e-16)

# 数值积分工具
function adaptive_quadrature(f, a, b; tolerance=1e-10)
function gauss_legendre_nodes(n::Int, a, b)

# 有限差分工具
function finite_difference(f, x; method=:central, order=1, h=1e-6)

end
```

#### 1.2 积分计算接口 (`src/core/integration.jl`)
```julia
module Integration

# 积分方法抽象
abstract type IntegrationMethod end
abstract type IntegrationGrid end

# 具体积分方法
struct GaussLegendre <: IntegrationMethod
    n_points::Int
    tolerance::Float64
end

struct AdaptiveQuadrature <: IntegrationMethod
    tolerance::Float64
    max_subdivisions::Int
end

# 积分网格定义
struct UniformGrid <: IntegrationGrid
    nodes::Vector{Float64}
    weights::Vector{Float64}
    domain::Tuple{Float64, Float64}
end

struct ProductGrid <: IntegrationGrid
    grids::Vector{IntegrationGrid}
end

# 核心积分接口
function integrate(method::IntegrationMethod, grid::IntegrationGrid, 
                  integrand::Function)
    
function integrate_multidimensional(method::IntegrationMethod, 
                                   grids::Vector{<:IntegrationGrid},
                                   integrand::Function)

# PNJL专用积分函数
function omega_thermal_integral(masses::Vector, chemical_potentials::Vector,
                               T::Float64, Phi1::Float64, Phi2::Float64,
                               momentum_grid::IntegrationGrid,
                               method::IntegrationMethod=GaussLegendre())

function vacuum_energy_integral(masses::Vector, momentum_grid::IntegrationGrid,
                               method::IntegrationMethod=GaussLegendre())

# 网格生成工具
function create_momentum_grid(n_points::Int, cutoff::Float64)
function create_angle_grid(n_points::Int)
function create_product_grid(grids...)

end
```

#### 1.4 常数管理模块 (`src/core/constants.jl`)
```julia
module PhysicsConstants

abstract type ConstantSet end

struct PNJLConstants <: ConstantSet
    hc::Float64
    rho0::Float64
    # ... 其他常数
end

struct GasLiquidConstants <: ConstantSet
    # 气液相变相关常数
end

# 常数访问接口
function get_constant(set::ConstantSet, name::Symbol)

end
```

### 2. Core Computation Layer (核心计算层)

#### 2.1 热力学计算核心 (`src/core/thermodynamics.jl`)
```julia
module ThermodynamicsCore

abstract type ThermodynamicState end
abstract type ThermodynamicSystem end

# 状态表示
struct SystemState{T<:Real} <: ThermodynamicState
    temperature::T
    chemical_potential::Vector{T}
    order_parameters::Vector{T}
end

# 通用热力学计算接口
function calculate_pressure(system::ThermodynamicSystem, state::ThermodynamicState)
function calculate_entropy(system::ThermodynamicSystem, state::ThermodynamicState)
function calculate_energy_density(system::ThermodynamicSystem, state::ThermodynamicState)

# 导数计算接口
function calculate_derivatives(system::ThermodynamicSystem, state::ThermodynamicState, 
                             variable::Symbol, order::Int=1)

end
```

#### 2.2 数值求解核心 (`src/core/solvers.jl`)
```julia
module NumericalSolvers

abstract type SolutionMethod end

struct NonlinearSolver <: SolutionMethod
    tolerance::Float64
    max_iterations::Int
    method::Symbol  # :newton, :nlsolve, etc.
end

# 通用求解接口
function solve_system(method::SolutionMethod, equations, initial_guess, parameters...)
function solve_with_continuation(method::SolutionMethod, equations, param_range, initial_state)

end
```

### 3. Physics Models Layer (物理模型层)

#### 3.1 PNJL 模型 (`src/models/pnjl/PNJLModel.jl`)
```julia
module PNJLModel

using ..ThermodynamicsCore
using ..PhysicsConstants

struct PNJL <: ThermodynamicSystem
    constants::PNJLConstants
    grid_config::GridConfiguration
    anisotropy_parameter::Float64
end

# PNJL 特定的状态
struct PNJLState <: ThermodynamicState
    temperature::Float64
    chemical_potentials::SVector{3,Float64}  # u, d, s
    condensates::SVector{3,Float64}         # φᵤ, φᵈ, φₛ  
    polyakov_loops::SVector{2,Float64}      # Φ₁, Φ₂
end

# PNJL 模型实现 (使用新积分接口)
function calculate_pressure(system::PNJL, state::PNJLState)
    # 解包状态参数
    T, mu, phi, Phi1, Phi2 = state.temperature, state.chemical_potentials, 
                             state.condensates, state.polyakov_loops...
    
    # 计算各个贡献项
    chi = calculate_chiral_condensate(phi)
    U = calculate_polyakov_potential(T, Phi1, Phi2)
    
    # 使用新积分接口计算费米子贡献
    masses = calculate_effective_masses(phi)
    
    # 真空能量积分
    vacuum_energy = vacuum_energy_integral(masses, system.grid_config.momentum_grid)
    
    # 热贡献积分 (替代原来的 calculate_log_sum)
    thermal_energy = omega_thermal_integral(masses, mu, T, Phi1, Phi2, 
                                          system.grid_config.momentum_grid)
    
    return -(chi + U + vacuum_energy + thermal_energy)
end

function calculate_chiral_condensate(system::PNJL, state::PNJLState)
function calculate_polyakov_potential(system::PNJL, state::PNJLState)
function calculate_fermionic_contribution(system::PNJL, state::PNJLState)

end
```

#### 3.2 气液相变模型 (`src/models/gas_liquid/GasLiquidModel.jl`)
```julia
module GasLiquidModel

using ..ThermodynamicsCore

struct GasLiquid <: ThermodynamicSystem
    constants::GasLiquidConstants
    interaction_parameters::Vector{Float64}
end

# 气液相变特定实现
function calculate_pressure(system::GasLiquid, state::ThermodynamicState)
function find_phase_boundary(system::GasLiquid, temperature_range)

end
```

### 4. Application Layer (应用层)

#### 4.1 高级计算接口 (`src/physics/HighLevelAPI.jl`)
```julia
module HighLevelAPI

# 温度-密度扫描
function temperature_density_scan(model::ThermodynamicSystem, 
                                T_range, ρ_range;
                                solver_config=default_solver_config())

# 相图计算
function phase_diagram(model::ThermodynamicSystem, 
                      T_range, μ_range;
                      resolution=(100, 100))

# 热力学量计算
function thermodynamic_quantities(model::ThermodynamicSystem, 
                                state::ThermodynamicState)

end
```

## 接口设计规范

### 1. 类型系统设计

#### 抽象类型层次
```julia
# 基础抽象类型
abstract type PhysicalQuantity end
abstract type SystemParameter end
abstract type ComputationMethod end

# 具体物理量
struct Pressure <: PhysicalQuantity end
struct Entropy <: PhysicalQuantity end
struct EnergyDensity <: PhysicalQuantity end

# 计算方法
struct AutomaticDifferentiation <: ComputationMethod end
struct FiniteDifference <: ComputationMethod end
```

#### 泛型接口设计
```julia
# 通用计算接口
function compute(quantity::PhysicalQuantity, 
               system::ThermodynamicSystem, 
               state::ThermodynamicState,
               method::ComputationMethod=AutomaticDifferentiation())

# 专门化实现
function compute(::Pressure, system::PNJL, state::PNJLState, ::AutomaticDifferentiation)
    # PNJL 压力的自动微分实现
end
```

### 2. 配置管理

#### 计算配置 (`src/core/config.jl`)
```julia
module Configuration

@kwdef struct GridConfig
    momentum_points::Int = 128
    angle_points::Int = 16
    momentum_cutoff::Float64 = 20.0
end

@kwdef struct SolverConfig
    tolerance::Float64 = 1e-12
    max_iterations::Int = 1000
    method::Symbol = :nlsolve
end

@kwdef struct NumericalConfig
    finite_difference_step::Float64 = 1e-6
    safe_log_threshold::Float64 = 1e-16
    max_exponential::Float64 = 700.0
end

end
```

### 3. 错误处理策略

#### 分层错误处理
```julia
# 基础异常类型
abstract type PNJLException <: Exception end

struct ConvergenceError <: PNJLException
    message::String
    iterations::Int
    residual::Float64
end

struct NumericalInstabilityError <: PNJLException
    function_name::String
    problematic_values::Vector{Float64}
end

# 错误处理接口
function handle_convergence_failure(error::ConvergenceError, fallback_strategy)
function validate_physical_constraints(state::ThermodynamicState)
```

## 数据流设计

### 计算流水线
```
输入参数 → 参数验证 → 状态初始化 → 方程求解 → 结果计算 → 输出格式化
    ↓           ↓            ↓           ↓          ↓           ↓
错误检查 → 约束检查 → 数值稳定性 → 收敛检查 → 物理检查 → 结果验证
```

### 状态管理
```julia
# 计算上下文
struct ComputationContext
    model::ThermodynamicSystem
    config::Configuration
    cache::Dict{String, Any}
    history::Vector{ThermodynamicState}
end

# 状态转换
function evolve_state(context::ComputationContext, 
                     current_state::ThermodynamicState,
                     new_parameters...)
```

## 性能优化策略

### 1. 编译时优化
- 使用类型稳定的函数
- 适当的 `@inline` 标记
- 静态数组用于小向量计算

### 2. 运行时优化
- 预分配内存
- 缓存重复计算结果
- 向量化操作

### 3. 算法优化
- 使用高效的数值算法
- 适应性求解策略
- 并行计算支持（未来）

## 扩展性设计

### 新模型集成
1. 定义模型结构体，继承 `ThermodynamicSystem`
2. 实现必要的接口函数
3. 添加模型特定的计算函数
4. 编写测试和文档

### 新计算方法集成
1. 定义方法结构体，继承 `ComputationMethod`
2. 实现 `compute` 函数的特化版本
3. 集成到现有计算流程中

## 测试架构

### 测试层次
1. **单元测试**: 每个函数独立测试
2. **集成测试**: 模块间交互测试  
3. **系统测试**: 完整计算流程测试
4. **基准测试**: 性能回归测试

### 测试组织
```
test/
├── unit/              # 单元测试
│   ├── test_math_utils.jl
│   ├── test_pnjl_model.jl
│   └── ...
├── integration/       # 集成测试
│   ├── test_model_solver.jl
│   └── ...
├── benchmarks/        # 性能测试
└── fixtures/          # 测试数据
```

---

**注意**: 这个架构设计是一个渐进式重构的目标。现有代码应该逐步迁移到这个架构，而不是一次性重写。
