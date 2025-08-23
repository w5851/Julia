# Debug script for EquationFactory
# 打印 create_equation_system_inplace 的返回值与示例调用

using Printf

# 引入工厂实现
include(joinpath(@__DIR__, "..", "src", "core", "equation_factory.jl"))
using .EquationFactory

# 定义简单子方程，仅接受局部需要的参数
f(x, y, T) = x - T*0.1
g(y, x, mu) = y - mu*0.01

println("Creating equation system using create_equation_system_inplace...")
F_param!, names, param_names, pack, unpack = create_equation_system_inplace(
    (f, (:x, :y), (:T,)),
    (g, (:y, :x), (:mu,))
)

println("Types:")
println(" F_param! => ", typeof(F_param!))
println(" names => ", typeof(names), " -> ", names)
println(" param_names => ", typeof(param_names), " -> ", param_names)
println(" pack => ", typeof(pack))
println(" unpack => ", typeof(unpack))

# 测试 pack/unpack
X = [1.23, 4.56]
nt = pack(X)
println("pack([1.23,4.56]) => ", nt, " (typeof: ", typeof(nt), ")")
println("unpack(nt) => ", unpack(nt))

# 测试 F_param! 调用
res = zeros(2)
params = (T=100.0, mu=10.0)
println("Calling F_param!(res, X, params) with X=", X, " params=", params)
F_param!(res, X, params)
println("res => ", res)

# 测试 bind_params
F = bind_params(F_param!, params)
res2 = zeros(2)
F(res2, X)
println("After bind_params F(res2, X) => ", res2)

println("Done.")
