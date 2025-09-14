# Gas_Liquid 脚本文件夹

这个文件夹包含了与气液相变程序相关的ForwardDiff自动微分计算脚本。

## 📁 文件说明

### 1. `forwarddiff_temperature_scan.jl` - 主计算脚本
**功能**: 使用ForwardDiff方法进行PNJL模型的温度扫描计算

**主要特性**:
- 🔥 **高精度自动微分**: 使用ForwardDiff计算一到四阶压强导数
- 📊 **温度扫描**: 固定重子化学势μ_B，扫描温度范围20-200 MeV
- 🎯 **热力学涨落**: 计算κ₃/κ₁和κ₄/κ₂比值，用于相变研究
- 💾 **数据输出**: 自动保存带元数据的CSV格式结果文件

**计算内容**:
- 压强归一化值 P/T⁴
- 一阶到四阶累积量 κ₁, κ₂, κ₃, κ₄
- 热力学涨落比值 κ₃/κ₁ 和 κ₄/κ₂
- 化学势温度比 μ/T

### 2. `plot_temperature_scan.jl` - 数据可视化脚本
**功能**: 读取温度扫描结果并生成可视化图表

**主要特性**:
- 📈 **数据可视化**: 绘制κ₃/κ₁和κ₄/κ₂随温度的变化曲线
- 📋 **元数据解析**: 自动读取CSV文件中的运行参数元数据
- 🖼️ **图像保存**: 生成高分辨率PNG图像文件
- 📊 **多图模式**: 支持比值图和各个κ值的独立展示

## 🚀 快速开始

### 第一步：运行温度扫描计算
```bash
# 在项目根目录下运行
cd D:\Desktop\Julia\Rotation_PNJL
julia --project=. scripts/Gas_Liquid/forwarddiff_temperature_scan.jl
```

**计算参数**:
- 重子化学势: μ_B = 697 MeV (固定)
- 温度范围: 20 - 200 MeV
- 温度步长: 1 MeV
- 输出文件: `output/Gas_Liquid/forwarddiff_temperature_scan.csv`

### 第二步：生成可视化图表
```bash
# 绘制结果图表
julia --project=. scripts/Gas_Liquid/plot_temperature_scan.jl
```

**输出图像**:
- `output/Gas_Liquid/kappa_ratios_temperature_scan.png` - κ比值图
- `output/Gas_Liquid/individual_kappas_temperature_scan.png` - 各κ值图

## 📋 详细使用说明

### `forwarddiff_temperature_scan.jl` 使用方法

#### 基本运行
脚本使用预设参数运行，无需额外配置：

```julia
# 直接运行脚本
julia --project=. scripts/Gas_Liquid/forwarddiff_temperature_scan.jl
```

#### 参数配置
可以在脚本中修改以下参数：

```julia
# 在脚本末尾的主程序部分
μ_B_fixed = 697.0/hc      # 重子化学势 (MeV)
T_min = 20.0/hc           # 最小温度 (MeV) 
T_max = 200.0/hc          # 最大温度 (MeV)
T_step = 1.0/hc           # 温度步长 (MeV)
```

#### 模型参数
PNJL模型参数设置在 `forwarddiff_temperature_scan()` 函数中：

```julia
# 模型参数
gsigma = 1.25     # sigma场初值
gdelta = 0.01     # delta场初值
fs = 17.28476     # sigma耦合常数
fo = 11.66174     # omega耦合常数
fr = 0.89363      # rho耦合常数
fd = 0.0          # delta耦合常数
b = 0.00210       # 三次项系数
c = -0.00297      # 四次项系数
```

### `plot_temperature_scan.jl` 使用方法

#### 基本运行（自动模式）
```julia
julia --project=. scripts/Gas_Liquid/plot_temperature_scan.jl
```
自动查找 `output/Gas_Liquid/forwarddiff_temperature_scan.csv` 文件并生成图表。

#### 指定输入文件
```julia
julia --project=. scripts/Gas_Liquid/plot_temperature_scan.jl custom_data.csv
```

#### 函数调用方式
```julia
# 在Julia REPL中使用
include("scripts/Gas_Liquid/plot_temperature_scan.jl")

# 绘制比值图
p1 = plot_temperature_scan_results("output/Gas_Liquid/forwarddiff_temperature_scan.csv", 
                                   "my_ratios.png")

# 绘制各个κ值
p2 = plot_individual_kappas("output/Gas_Liquid/forwarddiff_temperature_scan.csv", 
                            "my_kappas.png")
```

## 📊 输出文件格式

### CSV数据文件
**文件位置**: `output/Gas_Liquid/forwarddiff_temperature_scan.csv`

**元数据头部** (以 # 开头的注释行):
```
# ForwardDiff Temperature Scan Results
# Generated on: 2024-01-15T10:30:00.000
# Model Parameters:
# gsigma = 1.25
# gdelta = 0.01
# fs = 17.28476
# fo = 11.66174
# fr = 0.89363
# fd = 0.0
# b = 0.00210
# c = -0.00297
# mu_B = 697.0 MeV
# T_range = 20.0 - 200.0 MeV
# T_step = 1.0 MeV
# nodes = 256
```

**数据列说明**:
| 列名 | 单位 | 说明 |
|------|------|------|
| `T_MeV` | MeV | 温度 |
| `P_T4` | 无量纲 | 归一化压强 P/T⁴ |
| `kappa1` | 无量纲 | 一阶累积量 κ₁ |
| `kappa2` | 无量纲 | 二阶累积量 κ₂ |
| `kappa3` | 无量纲 | 三阶累积量 κ₃ |
| `kappa4` | 无量纲 | 四阶累积量 κ₄ |
| `kappa3_over_kappa1` | 无量纲 | 热力学涨落比值 κ₃/κ₁ |
| `kappa4_over_kappa2` | 无量纲 | 热力学涨落比值 κ₄/κ₂ |
| `mu_over_T` | 无量纲 | 化学势温度比 μ_B/T |

### 图像文件
**输出位置**: `output/Gas_Liquid/`

1. **`kappa_ratios_temperature_scan.png`**
   - 主要结果图：κ₃/κ₁ 和 κ₄/κ₂ 随温度变化
   - 包含水平参考线 (比值 = 1)
   - 高分辨率 (300 DPI)

2. **`individual_kappas_temperature_scan.png`**
   - 各个κ值的独立展示
   - 对数坐标显示
   - 用于诊断和详细分析

## ⚙️ 技术细节

### ForwardDiff自动微分
脚本使用Julia的ForwardDiff.jl包进行自动微分：

1. **一阶导数**: `ForwardDiff.derivative()`
2. **高阶导数**: 递归应用自动微分
3. **类型安全**: 使用类型提升确保数值稳定性
4. **增强四阶导数**: 使用5点中心差分提高精度

### 数值方法特点
- **约束方程求解**: 使用NLsolve.jl的Newton方法
- **类型保护**: 动态类型提升避免类型不匹配
- **错误处理**: 计算失败时返回NaN而不是中断
- **内存优化**: 预分配结果数组减少内存分配

### 物理意义
- **κ₃/κ₁**: 反映重子数涨落的偏度，相变信号
- **κ₄/κ₂**: 反映重子数涨落的峰度，临界点标识
- **μ_B/T**: 化学势温度比，控制相图位置

## 🔧 故障排除

### 常见问题及解决方案

**1. 包依赖问题**
```julia
# 确保安装了必要的包
using Pkg
Pkg.add(["NLsolve", "ForwardDiff", "CSV", "DataFrames", "Plots"])
```

**2. 输出目录不存在**
```bash
# 手动创建输出目录
mkdir -p output/Gas_Liquid
```

**3. 计算收敛失败**
- 调整模型参数 (gsigma, gdelta)
- 减小温度步长
- 检查化学势范围的合理性

**4. 绘图失败 (SSH环境)**
```julia
# 确保使用GR后端且保存文件而不显示
using Plots
gr()  # 设置后端
# 脚本会自动保存PNG文件
```

**5. 内存不足**
- 减少温度扫描范围
- 增大温度步长
- 监控系统内存使用

### 性能优化建议

**1. 计算性能**
- 使用较少的积分节点 (如128而不是256)
- 增大温度步长减少计算点数
- 在多核系统上可考虑并行化

**2. 数据管理**
- 定期清理output目录中的旧文件
- 使用压缩格式保存大数据集
- 备份重要的计算结果

## 📚 相关文档

- **主项目文档**: `../../docs/`
- **核心函数文档**: `../../src/Function_Gas_Liquid.jl`
- **测试文件**: `../../test/test_gas_liquid.jl`
- **使用指南**: `../../USAGE_GUIDE.md`

## 🤝 贡献指南

如需修改或扩展这些脚本：

1. **备份原始文件**
2. **测试修改后的功能**  
3. **更新文档和注释**
4. **添加适当的错误处理**
5. **保持代码风格一致性**

---
**最后更新**: 2025年9月14日  
**版本**: v1.0  
**维护者**: Rotation_PNJL项目组
