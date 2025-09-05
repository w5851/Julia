# Gas_Liquid 脚本文件夹

这个文件夹包含了与气液相变程序相关的示例脚本和演示程序。

## 文件说明

### 1. function_gas_liquid_examples.jl
完整的示例脚本，展示了 `Function_Gas_Liquid.jl` 中各种函数的用法，包括：
- 单点计算示例
- 批量化学势扫描
- 固定μ_B的温度扫描
- 热力学涨落计算

### 2. demo_temperature_scan.jl
专门用于演示温度扫描功能的脚本，展示如何：
- 固定重子化学势μ_B
- 改变温度T
- 计算不同温度下的κ₃/κ₁和κ₄/κ₂比值

## 使用方法

在项目根目录下运行：

```bash
# 运行完整示例
julia --project=. scripts/Gas_Liquid/function_gas_liquid_examples.jl

# 运行温度扫描演示
julia --project=. scripts/Gas_Liquid/demo_temperature_scan.jl
```

## 输出文件

所有脚本的输出文件都保存在 `../../output/` 目录下，格式为CSV：

- `pressure_derivatives_example.csv` - 压强导数批量计算结果
- `fluctuation_ratios_vs_T_example.csv` - 涨落比值温度扫描结果
- `demo_fluctuation_ratios_vs_T.csv` - 温度扫描演示结果

## CSV文件格式

### 压强导数文件格式
包含以下列：
- `mu_B_MeV`: 重子化学势 (MeV)
- `T_MeV`: 温度 (MeV)
- `pressure`: 压强
- `dpre_dmu1` 到 `dpre_dmu4`: 一到四阶导数
- `kappa1` 到 `kappa4`: 累积量
- `kappa2_over_kappa1`, `kappa3_over_kappa2`, `kappa4_over_kappa2`: 涨落比值

### 涨落比值文件格式
包含以下列：
- `T_MeV`: 温度 (MeV)
- `kappa3_over_kappa1`: κ₃/κ₁比值
- `kappa4_over_kappa2`: κ₄/κ₂比值

## 注意事项

- 确保在项目根目录下运行脚本
- 脚本会自动创建输出目录
- CSV格式便于后续的数据分析和绘图
