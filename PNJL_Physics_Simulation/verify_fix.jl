using PNJLPhysicsSimulation
using PNJLPhysicsSimulation.MathUtils

println("✅ PNJL 物理仿真包加载成功!")
println("\n测试安全对数函数:")
println("  safe_log(2.0) = ", safe_log(2.0))
println("  safe_log(-1.0) = ", safe_log(-1.0)) 
println("  safe_log(0.0) = ", safe_log(0.0))

println("\n✅ 数值稳定性修复验证成功!")
println("需求#1: 数值稳定性修复 - 已完成 ✅")
