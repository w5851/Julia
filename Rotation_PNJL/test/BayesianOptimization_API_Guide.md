# BayesianOptimization.jl 官方API使用指南

## 📖 概述

BayesianOptimization.jl 是一个专业的贝叶斯优化库，基于高斯过程的代理模型来优化昂贵的黑箱函数。

## 🚀 核心组件

### 1. 主要API结构

```julia
# 基本工作流
opt = BOpt(func, model, acquisition, modeloptimizer, lowerbounds, upperbounds; kwargs...)
result = boptimize!(opt)
```

### 2. 主要导出类型和函数

**核心函数：**
- `BOpt` - 贝叶斯优化器构造函数
- `boptimize!` - 执行优化
- `optimize` - 快速优化接口（使用默认参数）

**采集函数：**
- `ExpectedImprovement()` - 期望改进
- `ProbabilityOfImprovement()` - 改进概率
- `UpperConfidenceBound()` - 上置信界
- `ThompsonSamplingSimple()` - 简单汤普森采样
- `MutualInformation()` - 互信息

**模型优化器：**
- `MAPGPOptimizer` - 最大后验估计优化器
- `NoModelOptimizer` - 不优化模型超参数

**其他工具：**
- `ScaledSobolIterator`, `ScaledLHSIterator` - 初始化采样
- `Min`, `Max` - 优化方向
- `Silent`, `Timings`, `Progress` - 详细程度
- `maxduration!`, `maxiterations!` - 动态调整

## 🔧 详细用法

### BOpt 构造函数

```julia
BOpt(func, model, acquisition, modeloptimizer, lowerbounds, upperbounds;
     sense = Max,                    # 优化方向 (Max/Min)
     maxiterations = 10^4,           # 最大迭代次数
     maxduration = Inf,              # 最大运行时间
     acquisitionoptions = NamedTuple(), # 采集函数优化选项
     repetitions = 1,                # 每个点的重复评估次数
     verbosity = Progress,           # 输出详细程度
     initializer_iterations = 5*length(lowerbounds), # 初始采样点数
     initializer = ScaledSobolIterator(...))  # 初始化采样器
```

### MAPGPOptimizer 详细配置

```julia
# 关键：kernbounds 的正确格式
MAPGPOptimizer(
    every = 20,                     # 每20步优化一次超参数
    noisebounds = [-4, 3],          # 对数噪声边界 [下界, 上界]
    kernbounds = [
        [-3*ones(d); -3],           # 下界：[核参数下界..., 对数信号方差下界]
        [4*ones(d); 3]              # 上界：[核参数上界..., 对数信号方差上界]
    ],
    maxeval = 100                   # 超参数优化的最大评估次数
)
```

**重要说明：**
- `kernbounds` 是 `[下界向量, 上界向量]` 格式
- 对于 `SEArd` 核：需要 `d+1` 个参数（d个长度尺度 + 1个信号方差）
- 边界必须满足 `下界[i] <= 上界[i]`

### 高斯过程模型设置

```julia
using GaussianProcesses

# 创建弹性高斯过程模型
model = ElasticGPE(
    d,                              # 输入维度
    mean = MeanConst(0.0),         # 均值函数
    kernel = SEArd(zeros(d), 0.0), # 核函数：自动相关确定的平方指数核
    logNoise = -2.0,               # 对数噪声
    capacity = 3000                # 容量
)

# 预加载数据的情况
X_init = # d×n 矩阵（特征维度×样本数）
y_init = # n维向量
gp = GP(X_init, y_init, MeanConst(0.0), SEArd(ones(d), 0.0))
```

## 📋 完整示例

### 示例1：基本用法

```julia
using BayesianOptimization, GaussianProcesses

# 目标函数
f(x) = -(x[1] - 2)^2 - (x[2] + 1)^2 + 5

# 模型
model = ElasticGPE(2, mean = MeanConst(0.0), 
                   kernel = SEArd([0.0, 0.0], 0.0), 
                   logNoise = -2.0)

# 模型优化器 - 正确的kernbounds格式
modeloptimizer = MAPGPOptimizer(
    every = 10,
    noisebounds = [-4, 3],
    kernbounds = [[-2, -2, -3], [3, 3, 2]], # [x1_scale, x2_scale, signal_var]
    maxeval = 100
)

# 优化器
opt = BOpt(f, model, ExpectedImprovement(), modeloptimizer,
           [-5.0, -5.0], [5.0, 5.0],
           sense = Max,
           maxiterations = 50,
           verbosity = Progress)

# 执行优化
result = boptimize!(opt)
println("最优解: ", result.observed_optimizer)
println("最优值: ", result.observed_optimum)
```

### 示例2：预加载数据的热启动

```julia
# 已有的数据点
X_init = [[-1.0, 1.0], [0.0, 0.0], [2.0, -1.0]]  # 初始点列表
y_init = [f(x) for x in X_init]                   # 对应的函数值

# 转换为GP需要的格式
X_matrix = hcat(X_init...)  # 2×3 矩阵
y_vector = Vector{Float64}(y_init)

# 创建预加载的GP
gp = GP(X_matrix, y_vector, MeanConst(0.0), SEArd([1.0, 1.0], 0.0))

# 优化器（设置 initializer_iterations = 0）
opt = BOpt(f, gp, UpperConfidenceBound(), 
           NoModelOptimizer(),  # 使用固定超参数
           [-5.0, -5.0], [5.0, 5.0],
           sense = Max,
           maxiterations = 30,
           initializer_iterations = 0,  # 不进行额外初始化
           verbosity = Progress)

result = boptimize!(opt)
```

### 示例3：一维优化

```julia
# 一维函数
f_1d(x) = -(x[1] - 3)^2 + 5

# 一维初始数据
X_1d = reshape([1.0, 2.0, 4.0], 1, 3)  # 1×3 矩阵
y_1d = [f_1d([x]) for x in [1.0, 2.0, 4.0]]

# 一维GP
gp_1d = GP(X_1d, y_1d, MeanConst(0.0), SEArd([1.0], 0.0))

# 一维优化器
opt_1d = BOpt(f_1d, gp_1d, ExpectedImprovement(),
              MAPGPOptimizer(every = 5, 
                           noisebounds = [-4, 3],
                           kernbounds = [[-2, -3], [3, 2]]), # [length_scale, signal_var]
              [0.0], [6.0],
              sense = Max,
              maxiterations = 20,
              initializer_iterations = 0)

result_1d = boptimize!(opt)
```

## ⚙️ 高级配置

### 采集函数选项

```julia
# 采集函数优化配置
acquisitionoptions = (
    method = :LD_LBFGS,    # NLopt优化方法
    restarts = 10,         # 随机重启次数
    maxeval = 2000,        # 最大评估次数
    maxtime = 0.1          # 最大时间限制
)

opt = BOpt(f, model, acquisition, modeloptimizer, bounds...,
           acquisitionoptions = acquisitionoptions)
```

### 不同采集函数的特点

```julia
# 期望改进（最常用）
ExpectedImprovement()

# 上置信界（适合探索）
UpperConfidenceBound(scaling = BrochuBetaScaling(0.1), βt = 1.0)

# 改进概率（保守）
ProbabilityOfImprovement()

# 汤普森采样（随机）
ThompsonSamplingSimple()

# 互信息（理论最优）
MutualInformation()
```

## 🔍 常见问题和解决方案

### 1. kernbounds 边界错误
**错误**: `ArgumentError("invalid NLopt arguments: bounds 1 fail -1 <= 1 <= -2")`

**解决**: 确保下界 <= 上界
```julia
# 错误
kernbounds = [[-1, 3], [-2, 2]]  # -2 < 3 违反了第二个参数的边界

# 正确
kernbounds = [[-2, -3], [3, 2]]  # 所有下界都小于对应上界
```

### 2. GP 数据格式
```julia
# GP需要的数据格式
X = hcat(points...)     # d×n 矩阵（特征维度×样本数）
y = Vector{Float64}(values)  # n维向量

# 不是 n×d 矩阵！
```

### 3. 避免超参数优化问题
```julia
# 如果MAPGPOptimizer有问题，使用固定参数
opt = BOpt(f, gp, acquisition, NoModelOptimizer(), bounds...)
```

## 📚 参考资料

- [GitHub 仓库](https://github.com/jbrea/BayesianOptimization.jl)
- [官方文档](https://jbrea.github.io/BayesianOptimization.jl/dev/)
- Brochu et al. (2010): "A Tutorial on Bayesian Optimization"

## 💡 最佳实践

1. **开始简单**: 先用 `NoModelOptimizer` 测试基本功能
2. **数据格式**: 确保 X 是 d×n 格式，y 是向量
3. **边界设置**: 仔细检查 `kernbounds` 的边界关系
4. **采集函数**: `ExpectedImprovement` 通常是最好的起点
5. **调试**: 使用 `verbosity = Progress` 监控优化过程

## 🎯 总结

BayesianOptimization.jl 提供了完整的贝叶斯优化功能，关键是：
- 正确设置数据格式（d×n矩阵）
- 合理配置 `kernbounds`（下界 <= 上界）
- 选择合适的采集函数和模型优化器
- 从简单配置开始，逐步增加复杂性
