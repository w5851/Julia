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

### SolverInterface 新签名（仓库内统一 solver 工厂）

仓库内的 `SolverInterface.solve_equilibrium_equations` 已扩展以支持可选的 `solution_names` 参数，用于直接把解向量展开为 `NamedTuple`，便于 scanner 使用。签名如下：

```julia
function solve_equilibrium_equations(equation_system, initial_guess::Vector{T}, config::Any; solution_names::Union{Nothing, Vector{Symbol}}=nothing) where T<:Real
    # 若 solution_names 给定且长度匹配，返回 NamedTuple（字段由 solution_names 指定），并保留一个 :zero 字段
    # 否则返回 Vector{T}（原行为）
end
```

行为要点：
- 当调用方传入 `solution_names`（例如 `[:phi_u, :phi_d, :phi_s]`）且长度与解向量一致时，返回的值为一个 NamedTuple，包含按 `solution_names` 命名的字段；同时包含 `:zero` 字段保存原始向量（便于下游再次访问向量）。
- 若未传入 `solution_names`，函数维持向后兼容，仍然返回 `Vector{T}`。

示例：在 `scan_phase_space` 中使用（封装为符合 scanner 约定的 `solver(params, init; kwargs...)`）

```julia
using SolverInterface

function my_solver(params, init; kwargs...)
  config = Dict(:T => params.T, :mu => params.mu)
  # 传入 solution_names 以便直接获得命名字段
  res = SolverInterface.solve_equilibrium_equations((x,c)->equation_system(x,c), init, config; solution_names=[:phi_u])
  # res 现在是 NamedTuple，比如 (phi_u = 1.23, zero = [1.23])
  return merge(res, (converged=true, niter=1, residual=0.0))
end

scan_phase_space(f, axes; solver = my_solver, solution_names = [:phi_u])
```

说明：示例中 `my_solver` 使用了 `solve_equilibrium_equations(...; solution_names=...)`，使返回值直接包含可用于 CSV 列名的字段；这样 scanner 在构建 header 与写入数据时更可靠且无需额外映射。

- wrap_physics_fn(fn; mapping::Dict) -> f_wrapper：把项目已有函数适配为接受 `NamedTuple` 的接口。

注：实现应允许 `outputs` 同时由直接函数 f 返回和/或由 `solver` 返回；若 `outputs` 中某项在两者都有返回，应以 `solver` 的结果为准（可配置）。

## 解向量后处理与输出展开（postprocess 与 header 探测）

- 新增参数（API 草案）
  - `postprocess::Union{Nothing, Function} = nothing`  
    - 签名约定：`postprocess(solution, params) -> Dict{Symbol,Any}` 或 `NamedTuple`  
    - 语义：在 solver 成功返回解向量（或 `outvals`）后调用，用于从解向量计算派生的物理量（如 `:pressure`, `:entropy`）。返回键会并入最终输出列。
  - `solution_names::Union{Nothing,Vector{Symbol}} = nothing`  
    - 若 solver 返回单一向量（例如 `:zero`），可由调用方提供元素名称列表；否则默认展开为 `zero_1, zero_2, ...`。
  - `header_probe::Symbol = :first_point`（可选值 `:first_point` 或 `:first_batch`）
    - 语义：决定流式（stream）写入时如何确定完整 header（列名及向量展开列数）。
    - `:first_point`：对第一个格点先执行求解 + postprocess（探测列结构），写 header，然后写该点及后续点。  
    - `:first_batch`：先缓存 `probe_batchsize` 个点（默认 16），基于这批点构建 header，再写出 header + 缓存数据，然后继续流式写入（对不一致列/长度做检测）。
  - `probe_batchsize::Int = 16` （仅当 `header_probe == :first_batch` 有效）

- Solver / 解向量返回规范（约定）
  1. 优先：solver 返回 `NamedTuple`，其中标量字段直接作为列名写入；若包含 `:zero`（或其他向量字段），则根据 `solution_names` 或默认 `zero_1..zero_N` 展开为多列。  
  2. 其次：solver 返回 `Dict`，按键名映射列名。  
  3. 否则：若返回单个标量或向量，默认作为 `:Omega` 或 `:zero` 处理并展开（向量需提供 `solution_names` 否则使用 `zero_#` 名称）。  
  4. 若不同点间同一输出的向量长度不一致：首版实现拒绝并抛错（用户可更改为 JSON 编码或伴随二进制文件的后续实现）。

- postprocess 行为
  - 在 solver 成功后立即调用 `postprocess(solution, params)`（若提供），将返回的键并入 `outvals`，并写入 CSV。  
  - 若 `postprocess` 返回可变长度向量，同样按“拒绝可变长度”的规则处理。  
  - `postprocess` 可用于计算压强、熵、以及任何需要解向量作为输入的派生量。这样实现的职责划分清晰：solver 求解序参量，postprocess 计算热力学量。

- Header 与流式写入细节（实现者说明）
  - 选用 `header_probe = :first_point` 时，在写入任何 CSV 行前，对第一个格点执行 solver + postprocess（如果这一步失败，应根据 `fill_failed_with_last` 策略处理；若无可行结果，则无法决定列结构，报错或改用 `:first_batch`）。写好 header 后继续写第一个点及其余点。  
  - 选用 `header_probe = :first_batch` 时，缓冲前 `probe_batchsize` 点以决定列名与向量长度；缓冲后写 header 并将缓冲数据写出，若缓冲中出现不一致（例如向量长度冲突），在首版直接报错并停止，提示用户调整网格或提供 `solution_names`。  
  - header 包含：参数列（按 axes 顺序），解向量分量列（已展开），postprocess 产生的标量列，最后为 `solver_converged, n_iter, residual` 等 meta 列。

- meta.json 建议字段（应在 `<outfile>.meta.json` 中写入）
  - `axes`：每个 axis 的 name 与长度  
  - `outputs`：声明的 outputs 列（含展开细节）  
  - `solution_names`：若对向量做了命名展开，则记录名称列表与长度  
  - `total_points`、`created_at`、`tool_version`  
  - `failed_points`：若存在失败点，记录坐标/索引或参数值列表（可选）

- 最小调用示例（文档示例）
  - 假设已有 solver 与 postprocess：
    ```julia
    scan_phase_space(f, axes;
        solver = my_solver,
        initial_guess = my_init,
        postprocess = (solution, params) -> Dict(:pressure => compute_pressure(solution), :entropy => compute_entropy(solution)),
        outputs = [:Omega, :pressure, :entropy],
        solution_names = [:phi_u, :phi_d, :phi_s],
        header_probe = :first_point)
    ```
  - 语义：每点先由 `my_solver` 求解序参量（可返回 NamedTuple 或包含 `:zero` 向量），随后 `postprocess` 用该解计算 `:pressure` 与 `:entropy` 并写入 CSV。

- 错误/边界提示（行为规范）
  - 若 `postprocess` 抛异常，按 `fill_failed_with_last` 策略处理（默认填充或写 NaN），并把失败信息记录到 `meta.json`。  
  - 若 header 探测失败或出现向量长度不一致，首版将报错并提示使用更小的 probe 或提供 `solution_names`。

### 明确的返回值契约与兼容策略（必读）

为避免像 `Omega` 全部为 NaN 的情况，scanner 要求并推荐如下契约：

- 优先约定（推荐，强制化文档化）：
  1. 优先让 `solver` 返回 `NamedTuple`，并把该 NamedTuple 的键直接作为输出列名（例如 `(:phi_u, :phi_d, :pressure)`）。若包含向量字段（例如 `:zero`），应同时提供 `solution_names` 或直接返回解向量的命名字段（例如 `:phi_u, :phi_d, :phi_s`）。
  2. `postprocess` 若提供，也应返回 `NamedTuple` 或 `Dict`，其键会被并入最终输出列。

- 向后兼容（scanner 支持的回退）：
  1. 若 `solver` 返回单一向量（未命名），调用方必须提供 `solution_names` 以便展开为列；否则 scanner 将无法创建稳定列名并会抛错或写 NaN（视实现策略）。
  2. 若 `solver` 返回的 NamedTuple/Dict 中缺少 `:Omega`，scanner 不会自动臆断 Omega 的含义，除非调用方在 `postprocess` 中显式返回 `:Omega`，或明确约定 `solution_names` 的第一个分量作为 Omega（此行为应在调用时明确开启/接受）。

- Omega 自动填充策略（可选行为，需在调用时确认）：
  - scanner 可提供一个布尔参数（例如 `fill_Omega_from_solution=true`），在 `:Omega` 缺失时自动用 `solvec[1]` 或 `outvals[:zero][1]` 填充；推荐仅在 caller 明确同意时启用以避免语义歧义。

- warm-start 格式约定：
  - scanner 会把上一次成功的“解向量”以 `prev_solution_vector` 的形式保存并传给下一点作为 `g0`（若 `initial_guess` 返回 `nothing`）。为避免类型错配，`solver` 的 `init` 参数应期望一个与 `solvec` 相容的类型（通常是 `Vector{Float64}`）。不要使用 `outvals`（包含 Omega、pressure 等）直接作为初解，scanner 已在实现中区分 `prev_solution_vector` 与 `prev_solution_outvals`。

- meta.json 建议补充字段（必写或强烈建议）：
  - `solution_names`：当向量被展开时写入数组（用于后处理与列名再现）。
  - `fill_Omega_from_solution`：若启用自动填充，记录该布尔值。
  - `failed_points`：记录失败点（参数或索引）与失败原因。

- 快速 wrapper 示例（当旧 solver 只返回向量时，用 wrapper 返回 NamedTuple）：

```julia
# 假设旧 solver 返回 Vector{Float64}
function solver_wrapper(params, init; kwargs...)
    vec = old_solver(params, init; kwargs...)
    # 指定分量名
    return (phi_u = vec[1], phi_d = vec[2], zero = vec, converged=true)
end
```

文档结论：推荐调用方或仓库内的 solver 工厂统一返回 NamedTuple；scanner 仍保持对 `solution_names` 的兼容支持，但生产环境应采用 NamedTuple 约定以保证列名稳定与可复现。

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
A:包括"phi_u","phi_d"等与解向量(序参量)相关的字段；"pressure","rho"等与热力学量相关的字段
2. solver 签名与返回格式：仓库内是否已有统一 solver（名称/文件/示例）？期望返回哪些解分量与元信息字段？
A:solver的生成由core文件夹下的solver_interface中定义的接口函数生成，期望返回解向量(一般是[phi_u,phi_d,phi_s,Phi1,Phi2]，但后续有可能变动)
3. initial_guess 策略：是否使用邻点 warm-start（例如按行/列顺序重用前一个点的解）？默认是否启用？
A:是，注意换行/列时不能使用上一个行/列的最后一个点作为初解，因为它们在相空间上并非零点；默认启用
4. 收敛标准与容错策略：默认 tol、maxiter、失败是否中断扫描？
A:收敛标准默认，失败不中断扫描，但需要返回错误信息(至少需要知道哪个点失败了)，失败的这个点可以用上一个点的信息填充（或者跳过这个点）
5. 性能/并发约束：是否需要多线程或分布式优先？是否允许主线程写入并用 Channels 收集？
A:可以暂时不需要多线程，在之后的开发中再做打算
6. 向量输出的最大长度或展开规则：是否保证固定长度或需要特殊格式？
A:不固定长度，因为用户有时候需要输出压强，有时候不需要，因此需要输出的向量应该作为接口函数的输入让用户决定
7. 测试/验收条件：smoke test 的具体判定（除写文件外是否需要数值阈值或示例最小值位置匹配？）
A:暂时不需要
8. 兼容性：是否必须支持 ForwardDiff 对所有物理函数？还是仅在需要时用户启用？
A:物理函数应该作为该模块的输入传入，因此不需要考虑对自动微分的兼容

## 实现默认值（我将采用，除非你另行指定）

以下默认值将用于最小可行实现与 smoke test：

- solver 调用入口：假定注入式 `solver(params, init; kwargs...)`，实现中通过参数 `solver` 接受用户传入函数；若未提供则在 runner 中使用示例 solver。
- solver 返回结构：默认假定 solver 返回一个包含解向量的 NamedTuple（例如 `result.zero`），解向量内的分量顺序/含义由上层约定；额外的标量或热力学量（如压强）可以由额外的后处理函数或 wrapper 从解向量计算并加入输出列。

- outputs 序列化规则：定长向量展开为多列（例如 `order1, order2`）；首版实现拒绝可变长度向量并返回错误（可在未来版本通过 JSON 编码或伴随二进制文件支持）。
- initial_guess 与 warm-start：按行优先（行内前点作为 warm-start），默认启用；换行时不携带上一行末点的解，改用该行第一个点的默认初解或上一行的首点解作为初解。
- 收敛与容错默认值：使用求解器（例如 NLsolve）的默认收敛参数；求解失败时不终止全局扫描，默认用最近成功点的解填充该点（若无可用则写 NaN），并在 CSV 中将 `solver_converged=false`，在 meta.json 中记录失败点列表与错误信息；提供可选参数切换为“跳过该点”或“写 NaN”（不填充）。
- smoke test 默认 outputs：`[:Omega, :phi_u, :phi_d]`（可在 runner 中覆盖）。
- ForwardDiff 策略：scanner 不强制对物理函数做 AD 兼容性处理；若需自动微分，要求调用方传入接受 Dual 的物理函数或在 wrapper 中显式处理。
- smoke test 验收：检查 `output/phase_scan_sample.csv` 存在且行数为 1600，并在控制台打印全局最小 Omega 及其参数位置；不做数值阈值匹配。


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
