# PNJL 物理仿真项目需求文档

积分设计思路：在每个models文件夹的functions开头用.md格式写出Omega函数的核心公式，然后据此定义需要的被积函数，再使用公共的积分接口函数和求和接口函数进行计算
接口函数思路：由于Omega函数的多参数性质，接口函数可能需要允许其他参数的传入(被积函数的参数)，同时由于需要用到Omega的导数性质，接口函数要保证自动微分的兼容性

> 注：请先阅读 `agent/README.md` 以获得文档索引与推荐阅读顺序。
## #### 2. Omega函数积分接口重构
  - 创建通用积分接口：输入(nodes, weights, functions) → 输出(numerical results)
  - 将 `calculate_log_sum`, `calculate_energy_sum` 重构为标准积分调用
  - 实现模块化的积分计算核心
  - [x] 设计通用积分接口 `IntegrationInterface` 模块
  - [x] 实现基础积分方法和网格系统
  - [x] 实现物理专用积分函数
  - [x] 创建完整的单元测试套件（36个测试全部通过）
  - [x] 重构 PNJL 模型的积分调用（包装器方式完成）
  - [x] **已完成** 重构 PNJL_aniso 模型的积分调用（2025年8月18日完成）
  - [ ] 重构 Rotation 模型的积分调用
  - [ ] 验证重构后的数值一致性
  - ✅ 阶段一：积分接口基础架构 (100%完成)
  - 🟡 阶段二：各模型适配 (66%完成，2/3模型)  
  - ❌ 阶段三：数值一致性验证 (待开始)
  - `julia --project=. test/test_integration_interface.jl` - 36/36测试通过
  - `julia --project=. test/test_pnjl_aniso.jl` - 31/31测试通过

## 当前需求状态

### 🔴 高优先级需求（紧急）

#### 1. PNJL各向异性模型公共接口设计与实现
- **状态**: 🔄 实现中 (优先级最高)
- **描述**: 为PNJL各向异性模型实现4个核心公共接口，支持自动微分和方程求解
- **文件**: `src/core/pnjl_aniso_public_interface.jl` (新建)
- **具体接口需求**:
  - [ ] **Interface 1**: 自动微分接口 - 计算Omega函数对所有变量的偏导数
    - 输入: (φᵤ, φᵈ, φˢ, Φ₁, Φ₂), 温度T, 化学势μ, 各向异性参数ξ
    - 输出: (∂Ω/∂φᵤ, ∂Ω/∂φᵈ, ∂Ω/∂φˢ, ∂Ω/∂Φ₁, ∂Ω/∂Φ₂)
    - 要求: ∂Ω/∂Φ₁ = ∂Ω/∂Φ₂ = 0 (热力学平衡条件)
  - [ ] **Interface 2**: NLsolve包装接口 - 统一方程组求解器
    - 功能: 求解平衡态方程组 ∇Ω = 0
    - 支持: 自定义初值、收敛性控制、多种求解算法
  - [ ] **Interface 3**: 解评价接口 - 从平衡态计算物理量
    - 输入: 平衡解 (φ*, Φ₁*, Φ₂*)
    - 输出: 压力、能量密度、熵密度、重子密度等
  - [ ] **Interface 4**: 相图扫描接口 - 批量计算并导出结果
    - 功能: T-μ或T-ρ空间的系统扫描
    - 输出: 文件导出 (.dat/.csv格式)
- **技术要求**:
  - 完全基于已有的现代化functions.jl实现
  - 保持ForwardDiff兼容性
  - 统一错误处理和边界条件检查
  - 高性能设计，支持大规模计算
- **预期完成时间**: 2天
- **责任人**: 开发者
- **优先级**: 🔴 最高 (新功能开发的关键基础)

#### 2. 数值稳定性修复
- **状态**: ✅ 已完成
- **描述**: 修复 `calculate_U` 函数中的对数负值问题
- **文件**: `Function_PNJL_aniso.jl`
- **具体任务**:
  - [x] 识别 `log(value)` 可能出现负值的情况
  - [x] 实现安全对数函数 `safe_log`
  - [x] 替换所有不安全的对数调用
  - [x] 测试边界条件处理
- **完成时间**: 2025年8月18日
- **责任人**: 开发者
- **影响范围**: 核心计算稳定性

#### 2. Omega函数积分接口重构
- **状态**: � 阶段二已完成 (整体完成度90%)
- **描述**: 重构Omega函数实现，采用独立积分接口设计，消除闭包实现方式
- **文件**: `src/core/integration_interface.jl`, 各模型函数文件
- **设计方案**: 
  - 创建通用积分接口：输入(nodes, weights, functions) → 输出(numerical results)
  - 将被积函数从闭包重构为纯函数调用
  - 实现模块化的积分计算核心
- **具体任务**:
  - [x] 设计通用积分接口 `IntegrationInterface` 模块
  - [x] 实现基础积分方法和网格系统
  - [x] 实现物理专用积分函数
  - [x] 创建完整的单元测试套件（36个测试全部通过）
  - [x] 重构 PNJL 模型的积分调用（包装器方式完成）
  - [x] **已完成** 重构 PNJL_aniso 模型的积分调用（2025年8月19日完成）
  - [ ] 重构 Rotation 模型的积分调用
  - [x] 验证重构后的数值一致性（相对误差<0.01%）
- **阶段完成情况**:
  - ✅ 阶段一：积分接口基础架构 (100%完成)
  - ✅ 阶段二：PNJL_aniso模型重构 (100%完成，消除闭包)  
  - ❌ 阶段三：Rotation模型重构 (待开始)
- **测试状态**: 
  - `julia --project=. test/test_integration_interface.jl` - 36/36测试通过
  - `julia --project=. test/test_pnjl_aniso.jl` - 31/31测试通过
- **性能改进**: PNJL_aniso模型计算速度提升约2.3倍
- **完成报告**: `docs/pnjl_aniso_refactor_completion_report.md`
- **预期完成时间**: 0.5天 (仅剩Rotation模型)
- **优先级**: 高 (架构改进的基础)
- **最后更新**: 2025年8月19日

#### 3. 有限差分步长优化
- **状态**: 待开始
- **描述**: 优化 `central_fdm` 方法的步长控制
- **文件**: `Function_PNJL_aniso.jl`
- **具体任务**:
  - [ ] 分析当前步长选择的问题
  - [ ] 实现自适应步长或手动步长控制
  - [ ] 验证数值精度改进
  - [ ] 性能基准测试
- **预期完成时间**: 2天
- **依赖**: Omega函数积分接口重构完成

### 🟡 中优先级需求（重要）

#### 4. 代码架构重构
- **状态**: ✅ 常量分离已完成，设计阶段进行中
- **描述**: 将现有代码重构为高内聚低耦合的模块化架构
- **具体任务**:
  - [x] **常量分离重构**: 将所有与模型绑定的常量分离到各模型的constants.jl中
  - [ ] 设计抽象接口和类型系统
  - [ ] 拆分 `Function_PNJL_aniso.jl` 为多个功能模块
  - [ ] 实现 Core、Models、Physics 模块分离
  - [ ] 定义统一的计算接口
- **预期完成时间**: 1周
- **影响范围**: 整个项目结构
- **依赖**: Omega函数积分接口重构完成
- **最新进展**: 2025年8月18日完成常量分离，消除重复定义，建立模块化常数架构

#### 5. 接口函数文档化
- **状态**: 待开始
- **描述**: 为所有公开接口编写详细文档
- **具体任务**:
  - [ ] 识别所有公开接口函数
  - [ ] 编写标准化的函数文档
  - [ ] 创建使用示例
  - [ ] 添加类型签名和约束说明
- **预期完成时间**: 3天
- **依赖**: 架构重构进行中

### 🟢 低优先级需求（改进）

#### 5. 性能优化
- **状态**: 分析阶段
- **描述**: 优化关键计算路径的性能
- **具体任务**:
  - [ ] 使用 `@benchmark` 识别性能瓶颈
  - [ ] 优化内存分配和类型稳定性
  - [ ] 考虑并行计算可能性
  - [ ] 缓存计算结果
- **预期完成时间**: 1周
- **影响范围**: 计算效率

#### 6. 测试套件完善
- **状态**: 部分完成
- **描述**: 建立完整的测试覆盖
- **具体任务**:
  - [ ] 单元测试：每个函数独立测试
  - [ ] 集成测试：模块间交互测试
  - [ ] 边界测试：极限条件测试
  - [ ] 回归测试：确保修改不破坏现有功能
- **预期完成时间**: 1周

## 具体功能需求

### 核心计算模块需求

#### Omega函数积分接口重构
**设计目标**: 创建通用、可复用的积分计算接口，解决当前代码中积分逻辑分散和重复的问题。

**现状分析**:
- `calculate_log_sum` 函数在3个模型中重复实现
- 积分参数传递复杂：nodes, weights, coefficient 等多参数耦合
- 物理计算与数值积分逻辑混合，难以测试和维护

**重构方案**:
```julia
# 通用积分接口设计
module IntegralInterface

abstract type IntegrationMethod end
struct GaussLegendre <: IntegrationMethod end

# 核心积分接口
function integrate(method::IntegrationMethod, 
                  nodes::Vector, weights::Vector, 
                  integrand::Function) 
    return numerical_result
end

# 多维积分支持
function integrate_2d(method::IntegrationMethod,
                     p_nodes, p_weights, t_nodes, t_weights,
                     integrand::Function)
    return numerical_result
end

# PNJL模型特化
function omega_thermal_contribution(masses::Vector, mu::Vector, T, Phi1, Phi2,
                                  nodes, method::IntegrationMethod)
    integrand = (p, i) -> begin
        mass_i, mu_i = masses[i], mu[i]
        E = calculate_energy(mass_i, p)
        return calculate_log_term(E, mu_i, T, Phi1, Phi2)
    end
    
    return sum(i -> integrate(method, nodes[1], nodes[2], 
                             p -> integrand(p, i)), eachindex(masses)) * (-T)
end

end
```

**实现步骤**:
1. 创建 `src/core/integration_interface.jl` 模块
2. 重构现有 `calculate_log_sum`, `calculate_energy_sum` 为接口调用  
3. 统一节点权重管理：`get_nodes()` → `IntegrationGrid`
4. 实现向后兼容的包装函数
5. 添加完整的单元测试覆盖

#### 安全数学函数库
```julia
# 需要实现的安全函数
function safe_log(x; min_val=1e-16, handle_negative=:clamp)
function safe_exp(x; max_val=700.0)
function safe_sqrt(x; min_val=0.0)
```

#### 热力学计算接口
```julia
abstract type ThermodynamicSystem end
function calculate_pressure(system::ThermodynamicSystem, state...)
function calculate_derivatives(system::ThermodynamicSystem, state...)
```

### 物理模型需求

#### PNJL 模型抽象
```julia
struct PNJLModel <: PhysicsModel
    constants::PNJLConstants
    parameters::PNJLParameters
end

function solve_equations(model::PNJLModel, conditions...)
```

#### 相变计算模块
```julia
function phase_transition_analysis(model, T_range, rho_range)
function critical_point_finder(model)
```

## 技术债务清单

### 代码质量问题
1. **类型不稳定性**: `Function_PNJL_aniso.jl` 中存在类型推断问题
2. **硬编码常数**: 魔法数字散布在代码中，需要集中管理
3. **函数过长**: 部分函数超过100行，需要拆分
4. **错误处理不足**: 缺少异常情况的优雅处理

### 性能问题
1. **内存分配**: 循环中创建临时数组
2. **重复计算**: 相同表达式在循环中重复计算
3. **自动微分开销**: ForwardDiff 在某些场景下效率低

### 可维护性问题
1. **文档缺失**: 大部分函数缺少使用文档
2. **测试覆盖不足**: 核心功能缺少测试
3. **依赖混乱**: 模块间循环依赖

## 长期规划需求

### 第一阶段（当前）：稳定性和重构
- 修复数值稳定性问题
- 建立模块化架构
- 完善文档和测试

### 第二阶段：功能扩展
- 支持更多物理模型
- 实现图形化结果展示
- 添加并行计算支持

### 第三阶段：性能和易用性
- 全面性能优化
- 提供用户友好的 API
- 集成到 Julia 包生态系统

## 质量标准

### 代码质量指标
- **测试覆盖率**: > 80%
- **文档覆盖率**: 100%（所有公开接口）
- **类型稳定性**: 无类型推断警告
- **性能回归**: < 5%

### 文档标准
- 每个公开函数必须有完整的 docstring
- 包含参数说明、返回值、使用示例
- 物理背景和数学公式说明
- 边界条件和错误处理说明

---

## 更新日志
- **2025-08-18**: 初始需求文档创建
- **2025-08-18**: 识别对数函数负值问题，添加到高优先级需求
- **2025-08-18**: ✅ 完成常量分离重构
  - 移除重复的物理常数定义（π, hc, Nc等）
  - 建立通用物理常数模块(`core/constants.jl`)
  - 各模型常数独立管理：`models/{model}/constants.jl`
  - 消除 `UnifiedConstants` 模块的冗余
  - 清理过时测试文件和临时文件
  - 验证包成功加载和常数访问

---

**注意**: 本文档是活文档，每次代码更新后都需要相应更新需求状态。完成的需求标记为 [x]，进行中的需求保持 [ ] 并更新进度。
