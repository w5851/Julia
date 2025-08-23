# SolverInterface API 文档

概览
- 模块: `SolverInterface`（文件：`src/core/solver_interface.jl`）
- 目的: 封装对 NLsolve 的简易调用，约定使用 in-place 残差函数签名并在需要时把解向量 `result.zero` 通过 `pack` 转为 NamedTuple 返回，简化上层 scanner/工厂 的交互。

导出符号
- `solve_equilibrium_equations` — 主要入口，接受 in-place 残差函数并返回解或带命名字段的 NamedTuple。

接口契约
- solve_equilibrium_equations(equation_system, initial_guess::Vector{T}, config; pack=nothing)
  - 输入:
    - `equation_system`：必须为 in-place 残差函数，签名为 `equation_system(res::AbstractVector, x::AbstractVector, config)`。
      - `equation_system` 在内部应将残差写入 `res`（并可返回 `res`）。
    - `initial_guess::Vector{T}`：初始猜测向量，类型参数 `T<:Real`。
    - `config::Any`：传递给 `equation_system` 的参数容器（常见为 `NamedTuple` 或 `Dict`）。
    - 关键字参数:
      - `pack`（可选）：若提供，应为 `pack(::AbstractVector) -> NamedTuple` 的函数；当求解成功后，`solve_equilibrium_equations` 会调用 `pack(result.zero)` 并返回一个合并了 `zero` 字段的 NamedTuple（例如 `merge(nt, (zero = result.zero,))`）。
  - 输出:
    - 若未提供 `pack`，返回 `result.zero`（Vector）。
    - 若提供 `pack`，返回 `NamedTuple`（含解的命名字段以及 `:zero` 字段保存原始向量）。
  - 错误与异常:
    - 若 NLsolve 未收敛，会抛出错误：`error("求解未收敛: residual_norm=..., iterations=...")`。
    - 调用者可用 try/catch 捕获并执行自定义策略（例如记录失败点、返回 NaN、或重试）。

实现细节与约定
- 内部包装: 把传入的 `equation_system(res,x,config)` 包装为 `residual_func!(F,x)` 以适配 `nlsolve` 的接口。
- 自动微分: 当前 `nlsolve` 调用使用 `autodiff=:forward`（ForwardDiff）。因此，`equation_system` 在实现时须对 `res` 的类型敏感（在需要时使用 `zeros(eltype(x), n)` 分配残差以兼容 Dual）。
- pack/field 对齐: 推荐 `pack` 由 `EquationFactory.create_equation_system_inplace` 提供，保证 `pack(result.zero)` 返回的字段名与工厂返回的 `names` 一致，从而 CSV 列与解字段对齐。

示例
```julia
# 假设工厂返回 F_param!, names, param_names, pack, unpack
F_param!, names, param_names, pack, unpack = create_equation_system_inplace(...)
params = (T=1.0, mu=0.0)
init = zeros(length(names))
# 直接传入工厂返回的 in-place 残差函数和 pack
res = solve_equilibrium_equations(F_param!, init, params; pack=pack)
# res 是 NamedTuple，包含字段对齐的解和 :zero
```

兼容性与迁移提示
- 旧代码如果仍期望 `equation_system(x, config)` 返回残差向量，可用 `bind_params` 或自己包装成一个满足三参数签名的闭包：
```julia
# 把 F_param! 和 params 绑定为不带 params 的 in-place 函数
F = bind_params(F_param!, params)
# 或者显式包装为 (res,x,config)->F_param!(res,x,params)
```
- 若你的 solver/扫描器以前依赖 `solution_names`，现在应改为接受 `pack` 函数（由工厂提供），以确保字段命名的一致性。

扩展建议
- 暴露 `nlsolve` 的可选参数（如收敛容忍、最大迭代次数、线搜索策略）作为 `solve_equilibrium_equations` 的可选关键字参数，以便在不同问题上调优。
- 为失败点提供可配置返回模式（抛错 / 记录 NaN / 返回部分结果），方便 scanner 策略实现。

维护者注记：如需支持并行扫描或线程安全，请确认传入的 `config`/`params` 是线程安全的不可变值（例如 `NamedTuple`）并在 doc 中记录相关限制。
