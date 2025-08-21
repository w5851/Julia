## Unified Interface 重构计划

> 目标：将 `unified_physics_public_interface` 中的各接口按职责拆分为独立模块（各司其职），先对单个模块完成最小实现与单元测试，不做模块间组合或集成测试。

---

## 关键说明（你最新的要求）

- `solver_interface` 的核心功能**已实现**于 `src/core/equation_solver.jl`（作为 `EquationSolver` 并导出 `solve_equilibrium_equations`）。
- 对于 `solver_interface`，建议采取以下两种可选操作之一：
	- **重命名**：把 `equation_solver.jl` 文件重命名为 `solver_interface.jl`（保留或调整模块名）；或
	- **re-export**：在包入口或统一接口处把 `EquationSolver.solve_equilibrium_equations` 以 `solver_interface` 的 API 导出或别名化。

	两者区别：重命名会移动/改名文件；re-export 更轻量，只在接口层做符号映射。

---

## 拆分建议（模块列表与最小职责）

1) Solver（优先，已实现）
	 - 文件：`src/core/equation_solver.jl`（当前实现）
	 - 职责：方程组求解的低层实现，输入 residual 或热力学势构建的方程组，返回解向量或抛错。
	 - 建议操作：`重命名` 或 `re-export`（由你选择）。

2) Autodiff 接口（已重命名实现）
	- 文件：`src/core/autodiff_interface.jl`（已创建，基于原 `automatic_differentiation.jl`）
	- 职责：导出 `compute_gradient` / `compute_hessian` 的最小包装（当前直接使用原有实现逻辑，已通过兼容绑定暴露函数）。

3) Types
	 - 文件：`src/core/types.jl`
	 - 职责：集中定义 `PhysicalProperties`, `PhasePoint` 等类型，减少模块间循环依赖。

4) Properties 接口
	 - 文件：`src/core/properties_interface.jl`
	 - 职责：`calculate_physical_properties` 的最小实现（例如计算 pressure = -Ω，并返回 `PhysicalProperties`）。

5) Scan / IO（后置，可选）
	 - 文件：`src/core/scan_interface.jl`, `src/core/io_interface.jl`
	 - 职责：相图扫描与保存，最小实现只需记录成功/失败点并支持简单的 :dat 输出。

---

## 每个模块的「最小合同」示例

- Solver（已实现）
	- 输入：`equation_system::Function, initial_guess::Vector{T}, config::Any`
	- 输出：`Vector{T}`（解向量）
	- 错误：未收敛时抛出异常

- Autodiff
	- 输入：`f::Function, x::Vector{T}, config::Any`
	- 输出：梯度向量

- Properties
	- 输入：`solution_variables::Vector{T}, config::Any, thermodynamic_potential::Function`
	- 输出：`PhysicalProperties` 实例

---

## 最小单元测试策略（每个模块）

- 位置：`test/test_<module>_quick.jl`（例如 `test/test_solver_interface_quick.jl`）
- 内容：
	- Happy path（必做）：使用已知简单函数断言返回期望结果。
	- 错误用例（必做）：设计不可收敛或非法输入，断言抛出异常或按预期失败。
- 运行（PowerShell）：
```powershell
cd 'D:\Desktop\Julia\PNJL_Physics_Simulation'
julia --project=. test/test_solver_interface_quick.jl
```

---

## 推荐开发步骤（每次只做一个模块）

1. 你选择要先做的模块（推荐：Autodiff 或 Types，因 Solver 已就绪）。
2. 我创建对应模块文件（或执行 re-export/重命名），并实现最小函数签名。
3. 我添加 quick 单元测试并在 workspace 中运行验证。
4. 我把修改以补丁形式返回给你（不执行 git 操作）。

---

## 风险与缓解

- 循环依赖：把公共类型抽到 `types.jl` 并让其它模块引用该文件。
- 自动微分兼容性：测试使用 stub 或直接调用现有实现，先不改内部算法。

---

## 请选择下一步（回复其中一项）

- `Solver - re-export`：我将在包入口对现有 `EquationSolver` 做 re-export（最小改动）。
- `Solver - rename`：我将把 `equation_solver.jl` 重命名为 `solver_interface.jl` 并调整 include/export（文件移动）。
- `Autodiff`：我先实现 `autodiff_interface.jl` 并编写 quick 测试。
- `Types`：我先实现 `types.jl`（将数据类型抽出，降低耦合）。
- `Properties`：先实现并测试 `properties_interface.jl`。

回复你选择的项后，我将按计划开始实施并返回补丁与测试结果。