# EquationFactory API 文档

概览
- 模块: `EquationFactory`（文件：`src/core/equation_factory.jl`）
- 目的: 提供“方程组工厂函数”，把按子方程/变量声明的方程集合编译为可被求解器直接调用的 in-place 残差函数，并提供 `pack`/`unpack` 工具以保证解字段与 CSV 列名对齐。

导出符号
- `create_equation_system_inplace` — 工厂函数，生成 in-place 原型残差及相关工具。
- `bind_params` — 把带参数的 in-place 原型绑定成不带 params 的短闭包，方便传给只接受 `(res,x)` 签名的接口。

接口契约（契约式说明）
- create_equation_system_inplace(funcs_with_vars...; global_vars=nothing, param_names=nothing)
  - 输入:
    - `funcs_with_vars...`：若干项，每项为二元元组 `(func, vars)`。
      - `func`：一个函数，约定为 `func(x..., p...)`，其中前面的 `x...` 对应子方程所需的未知量（按 `vars` 中的顺序），后面的 `p...` 可选，为按 `param_names` 顺序展开的外部参数。
      - `vars`：未知量名字序列（每个元素可为 `Symbol`、`String` 或数字索引）；建议使用 `Symbol`，例如 `(:x,:y)`。
    - 关键字参数:
      - `global_vars`（可选）：显式指定全局未知量顺序（覆盖从子方程推断的顺序）。
      - `param_names`（可选）：外部参数名序列（例如 `(:T,:mu)`）；若不提供，工厂认为无外部参数。
  - 输出（返回五元组）:
    1. `F_param!(res::AbstractVector, X::AbstractVector, params)`：in-place 残差函数。计算并写入 `res`，返回 `res`。
       - `params` 通常为 `NamedTuple`，函数内部会按 `param_names` 顺序使用 `getfield(params, name)` 展开成位置参数传给 `func`。
       - 为了兼容自动微分（ForwardDiff），实现中对 `res` 的分配与 `X` 的元素访问做了注意（例如用 `zeros(eltype(x), n)` 在调用方创建残差向量）。
    2. `names`：全局未知量名字的 `Tuple{Symbol,...}`，对应 `X` 向量的索引顺序。
    3. `param_names`：参数名字的 `Tuple{Symbol,...}`（或空 Tuple）。
    4. `pack(::AbstractVector)`：把解向量转换成 `NamedTuple`，字段顺序与 `names` 保持一致（用于把 `result.zero` 转为可用字段名的结构，便于写 CSV / postprocess）。
    5. `unpack(::NamedTuple)`：把 NamedTuple 的值转换回普通 Vector（用于测试或逆操作）。

示例（基本用法）
```julia
# 子方程函数签名：func(x..., p...)
f(x,y,T,mu) = x + y - T
g(y,z,T,mu) = y - z + mu
h(x,z,T,mu) = x - z

F_param!, names, param_names, pack, unpack = create_equation_system_inplace((f, (:x,:y)), (g, (:y,:z)), (h, (:x,:z)); param_names=(:T,:mu))

# 调用示例（params 为 NamedTuple）
params = (T=1.0, mu=0.0)
res = zeros(eltype([0.0]), length(names)) # 或 zeros(eltype(x0), n)
X = [0.5, 0.5, 0.5]
F_param!(res, X, params)

# 把解向量打包为命名字段
nt = pack(X)  # NamedTuple(:x=>..., :y=>..., :z=>...)
```

辅助函数 `bind_params`
- 签名: `bind_params(F_param!, params)`
- 返回: 一个函数 `(res, X) -> F_param!(res, X, params)`。
- 用途: scanner 在每个扫描点可创建局部 `params` 并用 `bind_params` 生成短闭包传给只接受 `(res,x)` 签名的旧式 solver。

错误与异常
- 如果 `vars` 中声明的变量未包含在全局变量顺序中，会抛出错误并指出缺失变量名。
- 若 `funcs_with_vars` 为空，会抛出错误。

性能与自动微分注意事项
- 工厂返回的是 in-place 原型以最大限度减少分配。
- 实现对小长度的 `vars` 做了手写展开以避免创建临时 Tuple 的分配，从而对热点循环友好。
- 当配合 `ForwardDiff` 使用时，调用方应使用 `zeros(eltype(x), n)` 创建残差向量以兼容 Dual 类型。

兼容性建议
- 优先用 `Symbol` 指定变量名。若需要显式控制索引顺序，请使用 `global_vars`。
- 在 scanner 中推荐按“方案 B”在每个扫描点创建短生命周期的 `params`（NamedTuple），并用 `bind_params`/或直接把 `F_param!` 与 `params` 一起传入支持三参数签名的 solver。


----
文档维护者提示：若工厂 API 需要扩展（例如支持可选雅可比回调、稀疏模式或并行安全），请在此处记录向后兼容的变更策略。
