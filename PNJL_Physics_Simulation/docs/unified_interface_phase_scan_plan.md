# 相图扫描接口模块生成计划（CSV-only 精简版）

## 目标
在 `docs` 下提供一份面向实现的、可执行的生成计划，描述一个维度无关的相图扫描（phase-space scanner）接口模块的设计与契约，使用 CSV（CSV.jl + DataFrames.jl 可选，但实现必须在缺少包时回退为原生 CSV 写入）作为唯一持久化格式。

## 设计原则（精简）
- 单一职责：扫描器负责参数空间生成、调度与数据持久化；物理计算由独立函数提供。
- 维度无关：通过 `ParameterAxis`/`ParameterSpace` 抽象任意参数维度。
- 分块/流式 CSV：支持按批次写入 CSV 以控制内存峰值，并写伴随元数据文件（JSON）记录 `grid_shape`、`outputs`、`created_at`、`solver_options` 等。
- 兼容 ForwardDiff：物理函数应能处理 `ForwardDiff.Dual`，或 wrapper 提供数值差分备选。
- 多输出与求解器工作流：支持用户指定任意标量/标量向量类型输出（如 `Omega`, `pressure`, `order_params`），并支持在每个格点基于 `initial_guess` 调用迭代求解器来获得部分输出与求解器元信息。

## 数据结构（精简）
- `ParameterAxis`：name::String, values::Vector{Float64}, unit::Union{String,Nothing}
- `ParameterSpace`：Vector{ParameterAxis}；`generate_grid(axes; as_iterator=true)` 返回 Iterator{NamedTuple}
- `ScanResult`（内存模式）：DataFrame，列包含所有参数列、请求的输出列、以及 solver 元信息列（如 `solver_converged::Bool`, `n_iter::Int`, `residual`）。
- 持久化文件：主 CSV（每行一个格点，列为扁平化后的标量输出），伴随 JSON 元数据文件（例如 `phase_scan.meta.json`）记录 grid layout、outputs 列表、solver_options、生成时间与版本说明。

## 核心 API（精炼草案）
- scan_phase_space(f::Function, axes::Vector{ParameterAxis};
  mode::Symbol = :stream, # :stream 或 :memory
  outfile::String = "output/phase_scan.csv",
  outputs::Vector{Symbol} = [:Omega], # 请求的输出名
  solver::Union{Nothing, Function} = nothing, # 可选：单点求解器，签名参见下文
  initial_guess::Union{Nothing, Function} = nothing, # params -> init
  solver_options::Dict = Dict(),
  parallel::Bool = false,
  batchsize::Int = 128,
  max_points::Int = 1_000_000) -> ScanResult or outfile path

- generate_grid(axes::Vector{ParameterAxis}; as_iterator::Bool=true) -> Iterator{NamedTuple}

- evaluate_and_solve(grid_iterator; solver::Function, initial_guess::Function, solver_options::Dict)
  - 行为：对每个 `params` 调用 `g0 = initial_guess(params)`，随后 `solver(params, g0; solver_options...)`。
  - 要求：`solver` 返回 NamedTuple（包含用户关心的解分量）；同时返回 `solver_meta`：`(converged::Bool, niter::Int, residual::Float64, message::String)`。

- wrap_physics_fn(fn; mapping::Dict) -> f_wrapper：把项目已有函数适配为接受 `NamedTuple` 的接口。

注：实现应允许 `outputs` 同时由直接函数 f 返回和/或由 `solver` 返回；若 `outputs` 中某项在两者都有返回，应以 `solver` 的结果为准（可配置）。

## I/O 约定（CSV-only）
- 主输出为 CSV：每行对应一个格点，列依次为参数列（按 axes 顺序）、请求的输出列、以及 solver 元信息列。
- 向量型输出（长度固定）的每个分量应展开为单独列（例如 `order1, order2`）；长度可变的复杂输出应写入伴随二进制文件（非首选情况），但首版实现可拒绝可变长度输出并返回错误。
- 元数据写入 JSON（同目录，`<outfile>.meta.json`），包含：axes 描述、outputs 列表、scan 参数（batchsize、parallel）、total_points、created_at、tool_version。
- 若 CSV.jl/DataFrames.jl 不可用，实现必须使用轻量原生写入（open/print），保证跨平台可运行。

## 错误处理与边界条件
- NaN/Inf：在 CSV 对应格点写入 `NaN`；统计并在 JSON 元数据中报告计数。
- 点数阈值：若 total_points > max_points（默认 1e6），拒绝执行并返回错误，提示用户减小网格或采用自适应采样。
- 并发写入：禁止多线程直接写同一 CSV；建议工作线程将结果发回主线程，由主线程负责写入。

## 最简实现（目标导向的分析）
目标：用最少增量改动实现功能，使项目能在当前仓库基础上运行 smoke 测试（40×40 网格，输出 CSV 并打印全局最小 Omega），同时支持多输出与 per-point solver 驱动的基本功能。

最简目录改动（最小且安全）
- 新增/修改：`src/core/phase_scan_interface.jl`
  - 必需函数：`ParameterAxis`、`generate_grid`、`wrap_physics_fn`、`scan_phase_space`（实现 stream 与 memory 两种简单模式）、内置轻量 CSV 写入后备。
  - `scan_phase_space` 的最简行为：迭代 grid，按 batch 收集结果（或逐行写入），对每个点：
    1. 若提供 `initial_guess` 与 `solver`：调用 `g0 = initial_guess(params)` 并 `solver(params,g0;solver_options...)`，捕获 `solver_meta`。
    2. 否则调用 `f(params)`（或 wrap_physics_fn 返回）获取直接输出。
    3. 将请求的 `outputs` 从 solver 或 f 中抽取，写入 CSV（缺项写 NaN），并更新内存统计（min Omega）
- 新增演示脚本：`scripts/phase_scan_runner.jl`（已存在），更新为：演示多输出与 solver 模式（若项目有具体 solver，可调用；否则使用示例 solver）。

最小功能集（实现优先级）
1. 支持多输出的 CSV 写入（首要）
2. 逐点调用 `initial_guess` + `solver` 的工作流，并记录 `solver_meta`（次要但必需）
3. 批量/流式写入（避免占用大量内存）
4. 简单并发策略：工作线程仅计算，主线程写文件（可选实现）

实现示例契约（简短）
- `initial_guess(params::NamedTuple) -> NamedTuple`（或任意类型传给 solver）
- `solver(params::NamedTuple, init; solver_options...) -> NamedTuple`，返回解分量键名应与 `outputs` 的 subset 对齐；同时可返回 `converged,niter,residual`（或这些由 wrapper 捕获并封装为 `solver_meta`）

## 我可能缺少的信息（需要你确认或提供以完成代码实现）
1. outputs 具体集合与命名：除了 `:Omega`，常用的是哪些字段？字段名及是否允许向量型？
2. solver 签名与返回格式：仓库内是否已有统一 solver（名称/文件/示例）？期望返回哪些解分量与元信息字段？
3. initial_guess 策略：是否使用邻点 warm-start（例如按行/列顺序重用前一个点的解）？默认是否启用？
4. 收敛标准与容错策略：默认 tol、maxiter、失败是否中断扫描？
5. 性能/并发约束：是否需要多线程或分布式优先？是否允许主线程写入并用 Channels 收集？
6. 向量输出的最大长度或展开规则：是否保证固定长度或需要特殊格式？
7. 测试/验收条件：smoke test 的具体判定（除写文件外是否需要数值阈值或示例最小值位置匹配？）
8. 兼容性：是否必须支持 ForwardDiff 对所有物理函数？还是仅在需要时用户启用？

## Smoke test（精简）
- 要求：在 40×40（1600 点）网格上完成扫描，将结果写为 `output/phase_scan_sample.csv`，并在控制台输出全局最小 `Omega` 的位置与值。
- 运行示例（PowerShell）：

```powershell
julia --project=. scripts/phase_scan_runner.jl
```

## 依赖（精简）
- 可选：`CSV.jl`, `DataFrames.jl`（便于格式化与内存操作）
- 必须：实现必须在缺少这些包时回退为原生 CSV 写入（open/print）。

## 后续任务（建议）
- 基于此精简计划实现 `src/core/phase_scan_interface.jl` 的最小版本并运行 smoke test。实现完成后，可逐步添加并行、稳健的 CSV 分块写入、以及对大规模数据的后处理方案。

---

文件结束。
