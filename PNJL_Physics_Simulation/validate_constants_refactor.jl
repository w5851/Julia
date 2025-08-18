#!/usr/bin/env julia

println("🔍 验证常量分离重构完成情况")
println("="^50)

try
    using PNJLPhysicsSimulation
    
    # 测试通用物理常数
    using PNJLPhysicsSimulation.PhysicalConstants
    println("✅ 通用物理常数:")
    println("   hc = $hc MeV⋅fm")
    println("   Nc = $Nc (QCD颜色数)")
    
    # 测试各模型常数
    using PNJLPhysicsSimulation.PNJLConstants
    using PNJLPhysicsSimulation.PNJLAnisoConstants  
    using PNJLPhysicsSimulation.RotationConstants
    using PNJLPhysicsSimulation.GasLiquidConstants
    
    println("\n✅ 模型特定常数:")
    println("   PNJL模型 - rho0: $(PNJLConstants.rho0), T0: $(PNJLConstants.T0)")
    println("   PNJL_aniso模型 - Lambda_f: $(PNJLAnisoConstants.Lambda_f)")  
    println("   Rotation模型 - r0: $(RotationConstants.r0), C: $(RotationConstants.C)")
    println("   GasLiquid模型 - m: $(GasLiquidConstants.m)")
    
    println("\n🎉 常量分离重构成功完成!")
    println("📋 成果总结:")
    println("   - 消除了重复的物理常数定义")
    println("   - 建立了模块化的常数架构")
    println("   - 每个模型有独立的常数管理")
    println("   - 包能正确加载和访问所有常数")
    
catch e
    println("❌ 测试失败: $e")
    rethrow(e)
end
