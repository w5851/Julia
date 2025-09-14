# Gas_Liquid 模块文件分离说明

## 文件结构重组

为了更好地组织代码结构，我们将原来的 `Function_Gas_Liquid.jl` 文件分为两个部分：

### 1. `Function_Gas_Liquid.jl` - 基础函数模块
**包含内容**：
- 常量和依赖导入
- 基础数学函数（费米子分布函数等）
- 核心物理计算函数
- PNJL模型约束方程和求解器
- 压强计算的基础功能

**主要函数**：
- `get_nodes()` - 积分节点生成
- `fermion()` / `fermion_anti()` - 费米子分布函数
- `calculate_mass()` / `calculate_energy()` - 质量和能量计算
- `calculate_ρ()` / `calculate_ρ_s()` - 密度计算
- `calculate_*_term()` - 各种场项计算
- `calculate_fun_constraint()` - 约束方程
- `solve_fun_constraints()` - 约束方程求解器
- `calculate_pressure()` / `calculate_pressure_wrapper()` - 压强计算
- `calculate_pressure_solved()` - 带求解的压强计算

### 2. `Advanced_Gas_Liquid.jl` - 高阶函数模块
**包含内容**：
- 高阶导数计算功能
- 热力学涨落分析
- 批量处理和扫描功能
- 数据保存和输出功能

**主要函数**：
- `calculate_pressure_derivatives()` - 压强导数计算（通用版本）
- `calculate_pressure_derivatives_efficient()` - 高效导数计算
- `calculate_thermodynamic_fluctuations()` - 热力学涨落计算
- `calculate_derivatives_batch()` - 批量导数计算
- `calculate_fluctuation_ratios_vs_temperature()` - 温度扫描（标准版）
- `calculate_fluctuation_ratios_vs_temperature_advanced()` - 温度扫描（高级版）
- `save_derivatives_results()` - 导数结果保存
- `save_fluctuation_ratios_results()` - 涨落比值结果保存

## 使用方法

### 仅使用基础功能
```julia
include("src/Gas_Liquid/Function_Gas_Liquid.jl")

# 使用基础PNJL模型计算
nodes = get_nodes(256)
couplings = [17.28476, 11.66174, 0.89363, 0.0, 0.00210, -0.00297]
x0 = [1.25, 0.01, 0.35, 0.35]
pressure = calculate_pressure_solved(697.0/hc, 100.0/hc, x0, nodes, couplings)
```

### 使用高阶功能
```julia
include("src/Gas_Liquid/Advanced_Gas_Liquid.jl")  # 自动包含基础模块

# 计算压强导数
pressure_norm, dp1, dp2, dp3, dp4 = calculate_pressure_derivatives_efficient(
    697.0/hc, 100.0/hc, x0, nodes, couplings)

# 计算热力学涨落
kappa1, kappa2, kappa3, kappa4, ratios = calculate_thermodynamic_fluctuations(
    697.0/hc, 100.0/hc, x0, nodes, couplings)
```

## 优势

### 1. **模块化设计**
- 基础功能和高级功能分离
- 更清晰的代码结构
- 便于维护和扩展

### 2. **按需加载**
- 只需要基础计算时，无需加载复杂的高阶功能
- 减少内存占用和加载时间

### 3. **功能分级**
- **基础级**：PNJL模型核心计算
- **高级级**：导数、涨落、批量处理等

### 4. **依赖关系清晰**
- `Advanced_Gas_Liquid.jl` 依赖于 `Function_Gas_Liquid.jl`
- 基础模块可以独立使用

## 注意事项

1. **路径调整**：高阶模块中的输出路径已调整为正确的相对路径
2. **依赖导入**：高阶模块自动包含基础模块的所有功能
3. **向后兼容**：现有脚本只需修改 include 路径即可使用
4. **文档更新**：相关文档和示例需要相应更新

## 迁移指南

### 对于现有脚本
如果原来使用：
```julia
include("src/Gas_Liquid/Function_Gas_Liquid.jl")
```

现在需要根据使用的功能选择：

**仅使用基础功能**：
```julia
include("src/Gas_Liquid/Function_Gas_Liquid.jl")
```

**使用高阶功能**：
```julia
include("src/Gas_Liquid/Advanced_Gas_Liquid.jl")
```

### 对于新开发
- 建议先从基础模块开始
- 需要导数和涨落计算时再引入高阶模块
- 批量计算和扫描功能使用高阶模块

---
**分离完成时间**: 2025年9月14日  
**分离版本**: v1.0  
**维护者**: Rotation_PNJL项目组