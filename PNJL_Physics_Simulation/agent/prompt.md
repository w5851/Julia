# PNJL 物理仿真项目开发提示词

## 项目概述

本项目是一个基于 Julia 的物理仿真系统，专注于 PNJL (Polyakov Nambu Jona-Lasinio) 模型和旋转物理现象的计算。项目包含气液相变、手征对称性破缺、Polyakov Loop 势能计算等复杂物理过程的数值模拟。
项目文件夹名:PNJL_Physics_Simulation
项目注释使用中文

## 开发工作流程

### 1. 代码更新前的必读文档
每次开始编码前，必须按顺序阅读以下文档：
- `agent/requirements.md` - 当前需求和待实现功能
- `agent/architecture.md` - 项目架构和设计原则
- `agent/api_reference.md` - 所有接口函数文档
- `USAGE_GUIDE.md` - 项目使用指南

### 2. 代码更新后的文档维护
每次代码更新后，必须：
- 更新 `agent/requirements.md` 中对应的需求状态
- 如有新增接口函数，更新 `agent/api_reference.md`
- 如有架构变更，更新 `agent/architecture.md`
- 记录变更到 `agent/changelog.md`

## 代码设计原则

### 1. 高内聚低耦合架构
- **模块化设计**：每个物理过程独立为一个模块
- **接口抽象**：定义清晰的函数接口，隐藏实现细节
- **依赖注入**：通过参数传递依赖，避免硬编码
- **单一职责**：每个函数只负责一个具体功能

### 2. 抽象接口设计
```julia
# 示例：抽象计算接口
abstract type PhysicsModel end
abstract type ThermodynamicQuantity end

# 通用计算接口
function calculate(model::PhysicsModel, quantity::ThermodynamicQuantity, params...)
    error("Must implement calculate for $(typeof(model)) and $(typeof(quantity))")
end
```

### 3. 功能模块划分
- **Core 模块**：基础数值计算（积分、常数、热力学）
- **Models 模块**：物理模型实现（PNJL、气液相变、旋转）
- **Physics 模块**：高层物理计算接口
- **Utils 模块**：工具函数和辅助计算

### 4. 错误处理和数值稳定性
- 所有数学函数必须处理边界条件
- 使用安全的对数函数避免负值输入
- 实现数值稳定的算法（如安全的指数计算）
- 提供详细的错误信息和调试信息

## 接口函数编写规范

### 1. 函数签名规范
```julia
"""
函数功能简述

# 参数
- `param1::Type`: 参数1描述
- `param2::Type`: 参数2描述（可选，默认值）

# 返回值
- `ReturnType`: 返回值描述

# 示例
```julia
result = function_name(arg1, arg2)
```

# 注意事项
- 特殊情况处理
- 数值稳定性考虑
- 性能注意事项
"""
function function_name(param1::Type, param2::Type=default_value)
    # 实现
end
```

### 2. 类型安全和性能
- 使用类型注解提高性能
- 使用 `@inline` 优化热点函数
- 使用 `StaticArrays` 处理小向量
- 避免不必要的内存分配

### 3. 测试覆盖
每个接口函数必须有对应的测试：
- 单元测试：验证函数正确性
- 边界测试：测试极限情况
- 性能测试：确保计算效率

## 模块解耦指导

### 1. 数据流设计
```
输入参数 → 预处理 → 核心计算 → 后处理 → 输出结果
```

### 2. 依赖管理
- 核心计算模块不依赖具体物理模型
- 物理模型通过接口与核心模块交互
- 配置和常数集中管理

### 3. 可扩展性考虑
- 新物理模型通过实现标准接口即可集成
- 计算方法可插拔替换
- 支持不同精度和性能需求的实现

## 代码质量要求

### 1. 性能优化
- 使用 `@benchmark` 测量关键函数性能
- 避免类型不稳定的代码
- 合理使用编译时计算和缓存

### 2. 代码可读性
- 变量名称清晰表达物理意义
- 适当的注释解释物理背景
- 保持函数长度适中（建议 <50 行）

### 3. 版本控制
- 小步骤提交，每次提交有明确目的
- 提交信息包含变更类型和影响范围
- 重要功能更新需要标记版本

## 常见问题和解决方案

### 1. 数值稳定性问题
- 对数函数负值输入 → 使用 `safe_log` 函数
- 指数溢出 → 使用稳定的指数计算
- 除零错误 → 添加适当的检查和处理

### 2. 性能优化问题
- 自动微分性能 → 考虑使用有限差分或解析导数
- 内存分配过多 → 使用预分配和 in-place 操作
- 类型不稳定 → 添加类型注解和类型断言

### 3. 物理正确性问题
- 热力学一致性检查
- 单位制统一使用
- 物理约束条件验证

## 开发检查清单

### 代码开发前 ✓
- [ ] 阅读当前需求文档
- [ ] 了解相关接口定义
- [ ] 确认设计符合架构原则

### 代码开发中 ✓
- [ ] 遵循编码规范
- [ ] 实现错误处理
- [ ] 添加适当注释

### 代码开发后 ✓
- [ ] 编写/更新测试
- [ ] 更新文档
- [ ] 性能测试
- [ ] 更新需求状态

## 常用测试命令

### 积分接口测试
```bash
# 运行积分接口完整测试套件（36个测试）
julia --project=. test/test_integration_interface.jl

# 运行所有测试
julia --project=. test/run_all_tests.jl

# 运行核心功能测试
julia --project=. test/test_core_issues.jl
```

### 模块测试
```bash
# 测试PNJL模型
julia --project=. test/test_pnjl.jl

# 测试PNJL各向异性模型  
julia --project=. test/test_pnjl_aniso.jl

# 测试旋转模型
julia --project=. test/test_rotation.jl
```

## 正确使用 PNJLPhysicsSimulation 包的方法

### 1. 环境准备和激活
```bash
# 进入项目目录
cd "d:\Desktop\Julia\PNJL_Physics_Simulation"

# 激活项目环境并安装依赖包
julia --project=. -e 'using Pkg; Pkg.activate("."); Pkg.instantiate()'
```

### 2. 基本使用方法
```julia
# 方法1：作为开发包使用（推荐）
julia --project=.

# 在Julia REPL中
using PNJLPhysicsSimulation

# 检查包是否正确加载
println("Package loaded successfully!")

# 查看可用的模块
println(names(PNJLPhysicsSimulation))
# 输出：[:FunctionRegistry, :GasLiquidConstants, :GasLiquidFunctions, 
#        :Integration, :IntegrationInterface, :MathUtils, :ModelConfiguration, 
#        :PNJLAnisoConstants, :PNJLAnisoFunctions, :PNJLConstants, :PNJLFunctions, 
#        :PhysicalConstants, :RotationConstants, :RotationFunctions, :Thermodynamics]
```

### 3. 访问物理常数
```julia
using PNJLPhysicsSimulation

# 访问物理常数
π_value = PhysicalConstants.π              # π = 3.141592653589793
hc_value = PhysicalConstants.hc            # ℏc = 197.33 MeV⋅fm
```

### 4. 使用具体模块功能
```julia
using PNJLPhysicsSimulation

# 使用Gas-Liquid模型
# 通过GasLiquidConstants访问常数
# 通过GasLiquidFunctions访问函数

# 使用PNJL模型  
# 通过PNJLConstants访问常数
# 通过PNJLFunctions访问函数

# 使用PNJL各向异性模型
# 通过PNJLAnisoConstants访问常数  
# 通过PNJLAnisoFunctions访问函数

# 使用旋转模型
# 通过RotationConstants访问常数
# 通过RotationFunctions访问函数

# 使用积分功能
# 通过Integration访问数值积分函数

# 使用热力学功能
# 通过Thermodynamics访问热力学计算函数
```

### 5. 命令行快速测试
```bash
# 测试包加载
julia --project=. -e 'using PNJLPhysicsSimulation; println("Package loaded successfully!")'

# 测试物理常数
julia --project=. -e 'using PNJLPhysicsSimulation; println("π = ", PhysicalConstants.π); println("ℏc = ", PhysicalConstants.hc, " MeV⋅fm")'

# 查看导出的模块
julia --project=. -e 'using PNJLPhysicsSimulation; println("Available exports:"); println(names(PNJLPhysicsSimulation))'
```

### 6. 重要注意事项
- **必须在项目目录下使用** `--project=.` 参数激活项目环境
- **包名为** `PNJLPhysicsSimulation`（注意大小写）
- **所有功能都通过子模块访问**，如 `PhysicalConstants.π`、`Integration`、`PNJLConstants` 等
- **项目已配置为开发包**，修改源码后无需重新安装，直接重新加载即可

### 7. 开发工作流程
```julia
# 1. 激活环境
julia --project=.

# 2. 加载包
using PNJLPhysicsSimulation

# 3. 进行开发和测试
# 修改源码后，使用以下命令重新加载：
using Pkg; Pkg.activate("."); 
```

这是经过实际测试验证的正确使用方法，确保每次都能成功加载和使用包。

---

**记住：好的代码不仅要正确，还要易于理解、维护和扩展。每一行代码都是对未来的自己和同事的负责。**
