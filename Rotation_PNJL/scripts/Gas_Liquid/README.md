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
- 使用迭代初解策略提高收敛性

### 3. test_iterative_strategy.jl
性能对比测试脚本，比较：
- 固定初解策略 vs 迭代初解策略
- 收敛性和计算效率
- 结果精度分析

## 新功能：迭代初解策略

### 概述
新版本的温度扫描函数 `calculate_fluctuation_ratios_vs_temperature` 现在支持迭代初解策略：
- 每个温度点使用上一个温度点的收敛解作为初始猜测值
- 可选的线性外推功能，进一步改善初始猜测
- 显著提高计算效率和收敛稳定性

### 主要优势
1. **更快的收敛速度**：利用解的连续性，减少迭代次数
2. **更好的稳定性**：避免因初始猜测不当导致的收敛失败
3. **物理合理性**：利用了温度变化下物理量平滑变化的特性

### 函数版本
- `calculate_fluctuation_ratios_vs_temperature`：标准版本，默认使用迭代初解
- `calculate_fluctuation_ratios_vs_temperature_advanced`：高级版本，提供更多控制选项

### 使用示例
```julia
# 标准用法（推荐）
temperature_array, kappa3_over_kappa1, kappa4_over_kappa2, results_matrix = 
    calculate_fluctuation_ratios_vs_temperature(μ_B, T_min, T_max, x0, nodes, couplings)

# 高级用法
temperature_array, kappa3_over_kappa1, kappa4_over_kappa2, results_matrix, solution_history = 
    calculate_fluctuation_ratios_vs_temperature_advanced(μ_B, T_min, T_max, x0, nodes, couplings,
                                                        use_iterative_guess=true,
                                                        extrapolation_weight=0.3,
                                                        return_solution_history=true)
```

## 使用方法

在项目根目录下运行：

```bash
# 运行完整示例
julia --project=. scripts/Gas_Liquid/function_gas_liquid_examples.jl

# 运行温度扫描演示
julia --project=. scripts/Gas_Liquid/demo_temperature_scan.jl
```

## 输出文件

所有脚本的输出文件都保存在 `../../output/Gas_Liquid/` 目录下，格式为CSV：

- `pressure_derivatives_example.csv` - 压强导数批量计算结果
- `fluctuation_ratios_vs_T_example.csv` - 涨落比值温度扫描结果
- `demo_fluctuation_ratios_vs_T.csv` - 温度扫描演示结果
- `demo_fluctuation_ratios_vs_T_advanced.csv` - 高级温度扫描演示结果

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
- 脚本会自动创建输出目录 `output/Gas_Liquid/`
- CSV格式便于后续的数据分析和绘图
- 所有Gas_Liquid相关的输出文件都集中保存在专用文件夹中，便于管理
