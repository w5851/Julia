using PNJLPhysicsSimulation

# 计算 f(x)=x1^2 + x2^2 在点 [1.0,2.0] 的梯度，期望 [2.0,4.0]
g = PNJLPhysicsSimulation.compute_gradient(x->(x[1]^2 + x[2]^2), [1.0,2.0])
println("grad=", g)
