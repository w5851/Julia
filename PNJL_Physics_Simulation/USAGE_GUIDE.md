# PNJL Physics Simulation - 安装和使用指南

## 快速开始

### 1. 安装依赖包

首先启动 Julia 并安装必要的依赖包：

```julia
using Pkg

# 安装所需的包
Pkg.add("FastGaussQuadrature")
Pkg.add("ForwardDiff")
Pkg.add("NLsolve") 
Pkg.add("FiniteDifferences")
Pkg.add("BenchmarkTools")
Pkg.add("StaticArrays")
Pkg.add("SpecialFunctions")
Pkg.add("Test")
```

### 2. 运行测试

进入项目目录并运行测试：

```julia
# 在 Julia REPL 中
cd("d:/Desktop/Julia/PNJL_Physics_Simulation")
include("test/runtests.jl")
```

或者从命令行运行：

```bash
cd "d:\Desktop\Julia\PNJL_Physics_Simulation"
julia test/runtests.jl
```

### 3. 使用示例

#### 方法一：直接包含模块文件

```julia
# 进入项目目录
cd("d:/Desktop/Julia/PNJL_Physics_Simulation")

# 加载核心模块
include("src/core/constants.jl")
include("src/core/integration.jl")

# 加载气液模型
include("src/models/gas_liquid/constants.jl")
include("src/models/gas_liquid/functions.jl")

# 使用模块
using .PhysicalConstants
using .Integration
using .GasLiquidConstants
using .GasLiquidFunctions

# 运行计算
nodes = GasLiquidFunctions.get_nodes(256)
couplings = GasLiquidConstants.calculate_couplings(0.16, -16.0/197.33, 240.0/197.33, 0.75, 31.3/197.33)
```

#### 方法二：加载完整包（推荐）

```julia
# 进入项目目录
cd("d:/Desktop/Julia/PNJL_Physics_Simulation")

# 添加源目录到加载路径
push!(LOAD_PATH, "src")

# 加载主包
include("src/PNJLPhysicsSimulation.jl")
using .PNJLPhysicsSimulation
```

### 4. 示例计算

#### 气液相变模型示例

```julia
# 加载模块
include("src/models/gas_liquid/constants.jl")
include("src/models/gas_liquid/functions.jl")
using .GasLiquidConstants
using .GasLiquidFunctions

# 设置参数
nodes = get_nodes(256)
T = 50.0 / 197.33  # 50 MeV 温度
couplings = calculate_couplings(0.16, -16.0/197.33, 240.0/197.33, 0.75, 31.3/197.33)

# 计算压强导数
μ_B = 1001.0 / 197.33
x0 = [1.25, 0.01, μ_B/2, μ_B/2]
pressure, dpre_dmu1, dpre_dmu2, dpre_dmu3, dpre_dmu4 = 
    calculate_pressure_derivatives_efficient(μ_B, T, x0, nodes, couplings)

println("Pressure: $pressure")
println("First derivative: $dpre_dmu1")
```

#### PNJL 模型示例

```julia
# 加载模块
include("src/models/pnjl/constants.jl")
include("src/models/pnjl/functions.jl")
using .PNJLConstants
using .PNJLFunctions

# 设置参数
nodes = get_nodes(128)
T = 150.0 / 197.33
x = [-0.1, -0.1, -1.7, 0.5, 0.5]
mu = [320.0/197.33, 320.0/197.33, 320.0/197.33]

# 计算压强
pressure = pressure_solve_core(x, mu, T, nodes)
println("PNJL Pressure: $pressure")
```

## 故障排除

### 问题 1: "Package not found" 错误

**解决方案**: 不要使用 `using PNJLPhysicsSimulation`，而是使用直接包含的方法：

```julia
cd("d:/Desktop/Julia/PNJL_Physics_Simulation")
include("src/core/constants.jl")
# ... 其他模块
```

### 问题 2: FastGaussQuadrature 未安装

**解决方案**: 
```julia
using Pkg
Pkg.add("FastGaussQuadrature")
```

### 问题 3: 模块依赖错误

**解决方案**: 按顺序加载模块，先加载核心模块，再加载模型模块。

## 项目结构

```
PNJL_Physics_Simulation/
├── src/
│   ├── core/                    # 核心模块
│   │   ├── constants.jl         # 物理常数
│   │   ├── integration.jl       # 数值积分
│   │   └── thermodynamics.jl    # 热力学函数
│   ├── models/                  # 物理模型
│   │   ├── gas_liquid/          # 气液相变模型
│   │   ├── pnjl/                # PNJL 模型
│   │   └── pnjl_aniso/          # 各向异性 PNJL
│   └── PNJLPhysicsSimulation.jl # 主包文件
├── test/                        # 测试文件
├── scripts/                     # 示例脚本
└── docs/                        # 文档
```

## 更多示例

查看 `scripts/` 目录中的示例文件：
- `gas_liquid_example.jl` - 气液模型完整示例
- `pnjl_example.jl` - PNJL 模型示例

运行示例：
```julia
cd("d:/Desktop/Julia/PNJL_Physics_Simulation")
include("scripts/gas_liquid_example.jl")
```
