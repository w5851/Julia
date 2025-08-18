# PNJL 物理仿真项目 AI 开发代理系统

## 🤖 代理系统简介

本项目采用 AI 驱动的开发方式，通过结构化的提示词系统确保代码质量、架构一致性和文档完整性。

## 📁 代理文件结构

```
agent/
├── prompt.md          # 🎯 主开发提示词 - 核心开发指导原则
├── requirements.md    # 📋 需求文档 - 当前任务和优先级  
├── architecture.md    # 🏗️ 架构设计 - 模块化设计规范
├── api_reference.md   # 📚 API文档 - 所有接口函数参考
├── changelog.md       # 📰 变更日志 - 详细变更记录
└── quick_reference.md # ⚡ 快速参考 - 开发检查清单
```

## 🔄 标准开发工作流

### 开发前 (必读顺序)
1. **📋 `requirements.md`** - 了解当前需求和任务状态
2. **🏗️ `architecture.md`** - 确认架构设计原则  
3. **📚 `api_reference.md`** - 查看相关接口文档
4. **⚡ `quick_reference.md`** - 快速检查清单

### 开发中
- 遵循 `prompt.md` 中的编码规范
- 使用 `quick_reference.md` 中的代码模板
- 参考 `architecture.md` 确保设计一致性

### 开发后 (必须更新)
1. **📋 更新** `requirements.md` - 标记任务完成状态
2. **📚 更新** `api_reference.md` - 新增/修改的接口文档  
3. **📰 记录** `changelog.md` - 详细变更记录
4. **🧪 运行测试** - 确保功能正确性

## 🚨 当前紧急任务

根据 `requirements.md` 中的高优先级需求：

### 1. 数值稳定性修复
```julia
# 问题位置: Function_PNJL_aniso.jl 第58行
log_term = log(value)  # ⚠️ value 可能为负数

# 需要实现安全对数函数
function safe_log(x; min_val=1e-16, handle_negative=:clamp)
```

### 2. 有限差分步长优化  
```julia
# 当前代码
fdm = central_fdm(5, 1)  # 步长自动选择

# 需要支持手动步长控制
fdm = central_fdm(5, 1; step_size=1e-6)
```

## 🎯 架构重构目标

按照 `architecture.md` 设计，将现有代码重构为：

```
Foundation Layer     # 基础工具 (math_utils, constants)
       ↓
Core Computation     # 核心计算 (thermodynamics, solvers) 
       ↓
Physics Models       # 物理模型 (PNJL, gas_liquid, rotation)
       ↓  
Application Layer    # 高级接口 (用户友好API)
```

## 📊 代码质量标准

根据 `prompt.md` 要求：

- **测试覆盖率**: > 80%
- **文档覆盖率**: 100% (所有公开接口)
- **类型稳定性**: 无类型推断警告  
- **性能回归**: < 5%

## 🛠️ 开发工具集成

### Julia 包依赖
- `ForwardDiff.jl` - 自动微分
- `NLsolve.jl` - 非线性方程求解
- `FiniteDifferences.jl` - 有限差分  
- `StaticArrays.jl` - 高性能静态数组
- `BenchmarkTools.jl` - 性能测试

### 推荐开发实践
- 使用 `@code_warntype` 检查类型稳定性
- 使用 `@benchmark` 进行性能测试
- 使用 `Test.jl` 编写完整测试套件

## 📝 文档编写规范

所有函数必须包含完整的 docstring：

```julia
"""
函数功能简述

# 参数
- `param::Type`: 参数描述 (单位，范围)

# 返回值  
- `ReturnType`: 返回值描述

# 物理背景
相关物理概念和数学公式

# 使用示例
```julia
result = function_name(args...)
```

# 注意事项
- 数值稳定性考虑
- 边界条件处理
"""
```

## 🚀 快速开始

对于新的开发者：

1. **阅读核心文档**
   ```bash
   # 按顺序阅读这些文件
   agent/prompt.md          # 了解开发原则
   agent/requirements.md    # 了解当前任务
   agent/architecture.md    # 了解系统设计
   ```

2. **熟悉代码结构**
   ```bash
   # 主要源码文件
   src/Function_PNJL_aniso.jl  # 核心计算函数
   src/Constants_*.jl          # 物理常数定义
   ```

3. **运行测试**
   ```bash
   julia> include("test/runtests.jl")
   ```

## 📈 项目发展路线图

- **第一阶段**: 数值稳定性和基础重构 (当前)
- **第二阶段**: 模块化架构实现  
- **第三阶段**: 高级功能和用户接口
- **第四阶段**: 性能优化和生态系统集成

## 🤝 贡献指南

1. 每次修改前，确保阅读了最新的需求文档
2. 遵循 `architecture.md` 中的设计原则
3. 为所有新功能编写测试和文档
4. 更新相应的代理文档
5. 使用清晰的提交信息格式：`[类型] 简短描述`

## 📞 支持和反馈

- 技术问题：查阅 `api_reference.md`
- 架构疑问：参考 `architecture.md`  
- 任务规划：更新 `requirements.md`
- 变更追踪：查看 `changelog.md`

---

**⭐ 核心理念**: 通过结构化的 AI 代理系统，确保代码质量、架构一致性和持续改进。每一次代码修改都应该让整个项目变得更好、更稳定、更易维护。
