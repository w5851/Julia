# 温度差平方和计算函数使用说明

## 概述
基于 `Advanced_FindTforDiff.jl` 中的批量计算逻辑，新增了两个专门用于优化的函数，用于计算多组κ比值对应的温度差的平方和。

## 主要函数

### 1. `calculate_temperature_difference_sum_of_squares`

计算多组κ比值在给定优化参数下对应的温度差的平方和，返回标量值（MeV²单位）。

#### 函数签名
```julia
calculate_temperature_difference_sum_of_squares(
    kappa_pairs, μ_B, optimization_params, T_min, T_max;
    T_step_scan=1.0/hc, gsigma=1.25, gdelta=0.01, n_nodes=256,
    penalty_for_missing=1e6, verbose=false)
```

#### 参数说明
- `kappa_pairs`: κ比值对数组，格式 `[(κ₃/κ₁, κ₄/κ₂), ...]`
- `μ_B`: 重子化学势（无量纲单位）
- `optimization_params`: 优化参数元组 `(ρ₀, B_A, K, m_ratio, E_sym)`
- `T_min, T_max`: 温度搜索范围（无量纲单位）
- `penalty_for_missing`: 计算失败时的惩罚值（默认1e6）
- `verbose`: 是否显示详细信息（优化时建议设为false）

#### 返回值
- `sum_of_squares`: 所有组温度差平方和（MeV²单位）

### 2. `calculate_temperature_difference_sum_of_squares_with_weights`

计算加权温度差平方和，允许为不同的κ比值对分配不同权重。

#### 函数签名
```julia
calculate_temperature_difference_sum_of_squares_with_weights(
    kappa_pairs, weights, μ_B, optimization_params, T_min, T_max;
    T_step_scan=1.0/hc, gsigma=1.25, gdelta=0.01, n_nodes=256,
    penalty_for_missing=1e6, verbose=false)
```

#### 新增参数
- `weights`: 权重数组，与 `kappa_pairs` 长度相同

## 使用示例

```julia
# 包含模块
include("src/Gas_Liquid/Advanced_FindTforDiff.jl")

# 定义输入参数
kappa_pairs = [(1.2, 2.5), (1.5, 3.0), (1.8, 3.5)]
μ_B = 300.0 / hc  # 300 MeV
optimization_params = (0.15, 16.0, 240.0, 0.7, 32.0)  # (ρ₀, B_A, K, m_ratio, E_sym)
T_min, T_max = 80.0/hc, 200.0/hc  # 80-200 MeV

# 计算温度差平方和
sum_sq = calculate_temperature_difference_sum_of_squares(
    kappa_pairs, μ_B, optimization_params, T_min, T_max,
    verbose=false)

println("温度差平方和: $sum_sq MeV²")

# 使用权重版本
weights = [1.0, 2.0, 0.5]
weighted_sum_sq = calculate_temperature_difference_sum_of_squares_with_weights(
    kappa_pairs, weights, μ_B, optimization_params, T_min, T_max,
    verbose=false)

println("加权温度差平方和: $weighted_sum_sq MeV²")
```

## 运行演示

```julia
# 运行内置演示函数
demo_temperature_difference_sum_of_squares()
```

## 错误处理

函数具有完善的错误处理机制：
1. **计算失败**: 当某组κ比值的温度计算失败时，自动应用惩罚值
2. **无效温度**: 当找不到有效温度时，同样应用惩罚值
3. **优化友好**: 即使部分计算失败，函数仍能返回有效的目标值用于优化

## 优化应用

这些函数特别适用于：
1. **参数优化**: 作为优化算法的目标函数
2. **敏感性分析**: 评估参数变化对温度一致性的影响
3. **模型验证**: 检验理论模型的内在一致性

## 性能考虑

- 建议在优化过程中设置 `verbose=false` 以减少输出
- 可以调整 `penalty_for_missing` 值来控制优化行为
- `T_step_scan` 参数影响计算精度和速度的平衡