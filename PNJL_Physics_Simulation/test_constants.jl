#!/usr/bin/env julia

println("正在测试包加载...")
try
    using PNJLPhysicsSimulation
    println("✅ 包成功加载！")
    
    # 测试常数访问
    using PNJLPhysicsSimulation.PhysicalConstants
    using PNJLPhysicsSimulation.PNJLConstants
    
    println("hc = ", hc)
    println("Nc = ", Nc)
    println("rho0 = ", rho0)
    println("T0 = ", T0)
    
catch e
    println("❌ 错误: ", e)
    showerror(stdout, e)
    println()
end
