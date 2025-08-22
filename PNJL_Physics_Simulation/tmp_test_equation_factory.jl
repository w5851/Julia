# 临时测试脚本：加载 EquationFactory 并运行小示例
f(x, y) = x^2 + y^2 - 1      # 单位圆
g(y, z) = y - z              # 直线 y = z
h(x, z) = x + z - 2          # 平面 x + z = 2

include("src/core/equation_factory.jl")

F!, names, pack, unpack = EquationFactory.create_equation_system_inplace((f, (:x,:y)), (g, (:y,:z)), (h, (:x,:z)))
println("names=", names)

X = [0.6, 0.8, 0.8]
res = zeros(3)
F!(res, X)
println("res=", res)
println("pack(X)=", pack(X))
println("unpack(pack(X))=", unpack(pack(X)))
