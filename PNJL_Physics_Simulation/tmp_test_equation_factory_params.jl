# 临时测试：验证带 params 的工厂和点级闭包绑定（方案B）
include("src/core/equation_factory.jl")
using .EquationFactory

# 定义一个简单线性方程组（有解析解），依赖未知量 x,y 和参数 T, mu
# 方程：
#   f: x + y - T == 0
#   g: y - z + mu == 0
#   h: x - z == 0
# 解法（解析解）： z = (T + mu)/2, x = z, y = z - mu

f(x,y,T,mu) = x + y - T
g(y,z,T,mu) = y - z + mu
h(x,z,T,mu) = x - z

# 注意：这里函数签名为 func(x..., p...), 所以把参数写在最后
F_param!, names, param_names, pack, unpack = EquationFactory.create_equation_system_inplace((f, (:x,:y)), (g, (:y,:z)), (h, (:x,:z)); param_names=(:T,:mu))
println("names=", names, " param_names=", param_names)

include("src/core/solver_interface.jl")
using .SolverInterface

# 两组不同的参数点
params1 = (T=1.0, mu=0.0)
params2 = (T=0.36, mu=0.2)

# 初始猜测
x0 = [0.5, 0.5, 0.5]

for (i, p) in enumerate((params1, params2))
    println("--- point ", i, ": ", p)
    # 每点创建短闭包用于采样残差
    F!(res, X) = F_param!(res, X, p)
    res0 = zeros(3)
    F!(res0, [0.6,0.8,0.8])
    println("sample residual=", res0)

    # 为兼容现有 SolverInterface 定义一个返回 Vector 的 equation_system(x, config)
    # 这里直接构造 in-place 版本的 equation_system(res, x, cfg)
    function equation_system(res, x, cfg)
        F_param!(res, x, cfg)
        return res
    end

    try
    # 直接把工厂返回的 in-place F_param! 传给求解器，并把 pack 函数交给求解器用于生成 NamedTuple
    result = SolverInterface.solve_equilibrium_equations(F_param!, x0, p; pack=pack)
        println("result=", result)
    catch e
        println("solver error: ", e)
    end
end
