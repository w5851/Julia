# 快速开发参考
> 注：建议先阅读 `agent/README.md` 获取索引与推荐阅读顺序。

## 🚀 包使用快速指南

### 激活项目环境
```bash
cd "d:\Desktop\Julia\PNJL_Physics_Simulation"
julia --project=. -e 'using Pkg; Pkg.activate("."); Pkg.instantiate()'
```

### 基本使用
```julia
julia --project=.
using PNJLPhysicsSimulation

# 访问物理常数
PhysicalConstants.π        # π值
PhysicalConstants.hc       # ℏc = 197.33 MeV⋅fm

# 可用模块
# [:FunctionRegistry, :GasLiquidConstants, :GasLiquidFunctions, 
#  :Integration, :IntegrationInterface, :MathUtils, :ModelConfiguration, 
#  :PNJLAnisoConstants, :PNJLAnisoFunctions, :PNJLConstants, :PNJLFunctions, 
#  :PhysicalConstants, :RotationConstants, :RotationFunctions, :Thermodynamics]
```

### 快速测试
```bash
# 测试包加载
julia --project=. -e 'using PNJLPhysicsSimulation; println("Package loaded successfully!")'
```

## 每次开发前的检查清单 ✓

1. **📖 读取需求** - 查看 `agent/requirements.md` 当前待处理任务
2. **🏗️ 检查架构** - 确认修改符合 `agent/architecture.md` 设计原则  
3. **📋 查看API** - 了解相关接口规范 `agent/api_reference.md`
4. **⚠️ 识别问题** - 检查已知技术债务和限制条件

## 每次开发后的更新清单 ✓

1. **📝 更新需求** - 在 `agent/requirements.md` 中标记完成状态
2. **📚 更新文档** - 如有新接口，更新 `agent/api_reference.md`
3. **📰 记录变更** - 在 `agent/changelog.md` 中记录修改
4. **🧪 运行测试** - 确保修改没有破坏现有功能

## 当前紧急问题 🚨

### 1. 数值稳定性 - `calculate_U` 函数
```julia
# 问题代码 (Function_PNJL_aniso.jl:58)
log_term = log(value)  # value 可能为负数

# 需要修复为:  
log_term = safe_log(value, min_val=1e-16)
```

### 2. 步长控制 - `central_fdm`
```julia
# 当前代码
fdm = central_fdm(5, 1)  # 步长自动选择

# 需要支持:
fdm = central_fdm(5, 1; step_size=1e-6)  # 手动步长
```

## 常用代码模板

### 安全数学函数模板
```julia
@inline function safe_log(x; min_val=1e-16, handle_negative=:clamp)
    if x > min_val
        return log(x)
    elseif x < 0 && handle_negative == :clamp
        return log(min_val)
    else
        error("Invalid input for safe_log: $x")
    end
end
```

### 函数文档模板
```julia
"""
函数功能简述

计算具体的物理量，处理特定的数值情况

# 参数
- `x::Type`: 参数描述，包括单位和范围
- `config::Type`: 配置参数（可选）

# 返回值
- `ReturnType`: 返回值描述和单位

# 物理背景
简述相关的物理概念和数学公式

# 使用示例
```julia
result = function_name(input_data, config)
println("Result: \$result")
```

# 注意事项
- 数值稳定性考虑
- 边界条件处理
- 性能特征

# 参见
相关函数或文档链接
"""
function function_name(x::Type, config::Type=default_config())
    # 参数验证
    @assert x > 0 "x must be positive"
    
    # 核心计算
    result = compute_something(x, config)
    
    # 结果验证
    validate_result(result)
    
    return result
end
```

## 性能优化检查点

### Julia 性能最佳实践
- ✅ 类型稳定性：所有函数返回类型可推断
- ✅ 避免全局变量：使用参数传递
- ✅ 内存预分配：循环外分配数组  
- ✅ 使用 `@inbounds` 和 `@simd`：优化热点循环
- ✅ 静态数组：小向量使用 `SVector`

### 数值计算最佳实践
- ✅ 避免下溢：使用安全数学函数
- ✅ 控制上溢：限制指数计算范围
- ✅ 精度控制：合适的容差设置
- ✅ 稳定算法：选择数值稳定的计算方法

## 测试驱动开发模式

### 1. 编写测试
```julia
@testset "safe_log function tests" begin
    # 正常情况
    @test safe_log(1.0) ≈ 0.0
    @test safe_log(exp(1)) ≈ 1.0
    
    # 边界情况  
    @test safe_log(1e-16) ≈ log(1e-16)
    @test safe_log(-1.0) ≈ log(1e-16)  # 默认clamp行为
    
    # 异常情况
    @test_throws ArgumentError safe_log(-1.0, handle_negative=:error)
end
```

### 2. 实现功能
```julia
function safe_log(x; min_val=1e-16, handle_negative=:clamp)
    # 实现代码...
end
```

### 3. 验证测试
```julia
julia> using Test
julia> include("test_safe_log.jl")
Test Summary: | Pass  Total
safe_log function tests |    5      5
```

## 调试和诊断工具

### 类型稳定性检查
```julia
julia> @code_warntype calculate_pressure(phi, Phi1, Phi2, mu, T, nodes_1, nodes_2, xi)
```

### 性能分析
```julia
julia> using BenchmarkTools
julia> @benchmark calculate_pressure($phi, $Phi1, $Phi2, $mu, $T, $nodes_1, $nodes_2, $xi)
```

### 内存分配检查
```julia
julia> @time calculate_pressure(phi, Phi1, Phi2, mu, T, nodes_1, nodes_2, xi)
```

## 快速修复指南

### 常见错误及解决方案

#### 1. `DomainError: log will only return a real result if its argument is positive`
```julia
# 错误代码
result = log(value)

# 修复方案
result = safe_log(value)
```

#### 2. `MethodError: no method matching`
```julia
# 可能的类型不匹配问题
# 检查函数签名和参数类型
julia> typeof(your_argument)
```

#### 3. `BoundsError: attempt to access`
```julia  
# 数组越界
# 添加边界检查或使用 @boundscheck
@boundscheck 1 <= i <= length(array)
```

## 项目结构导航

```
Rotation_PNJL/
├── agent/          # 🤖 AI代理提示词和文档
│   ├── prompt.md   # 主开发提示词
│   ├── requirements.md  # 需求跟踪
│   ├── architecture.md  # 架构设计
│   ├── api_reference.md # API文档
│   └── changelog.md     # 变更日志
├── src/            # 📁 源代码
│   ├── Function_PNJL_aniso.jl  # 主计算函数
│   ├── Constants_*.jl          # 物理常数
│   └── Function_*.jl           # 各模型函数
└── test/           # 🧪 测试文件
    └── *.jl        # 各种测试
```

---

**💡 记住**: 
- 小步骤，频繁提交
- 测试驱动开发  
- 文档与代码同步更新
- 遵循既定架构原则
