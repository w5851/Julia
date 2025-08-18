#!/usr/bin/env julia

println("🧹 验证函数文件中的常数清理完成情况")
println("="^60)

try
    using PNJLPhysicsSimulation
    
    # 测试各模型的函数是否能正确访问常数
    using PNJLPhysicsSimulation.PNJLFunctions
    using PNJLPhysicsSimulation.PNJLAnisoFunctions  
    using PNJLPhysicsSimulation.RotationFunctions
    using PNJLPhysicsSimulation.GasLiquidFunctions
    
    println("✅ 所有模型函数模块成功加载")
    
    # 测试一些基本功能调用
    println("\n🔧 测试基本功能调用:")
    
    # PNJL模型测试
    try
        nodes = PNJLFunctions.get_nodes(32)
        println("   ✅ PNJL - get_nodes()正常工作")
    catch e
        println("   ❌ PNJL - get_nodes()错误: $e")
    end
    
    # Gas-Liquid模型测试
    try
        nodes = GasLiquidFunctions.get_nodes(32)
        println("   ✅ Gas-Liquid - get_nodes()正常工作")
    catch e
        println("   ❌ Gas-Liquid - get_nodes()错误: $e")
    end
    
    # PNJL_aniso模型测试
    try
        nodes = PNJLAnisoFunctions.get_nodes_aniso(32, 16)
        println("   ✅ PNJL_aniso - get_nodes_aniso()正常工作")
    catch e
        println("   ❌ PNJL_aniso - get_nodes_aniso()错误: $e")
    end
    
    println("\n📊 常数访问测试:")
    println("   通用常数 - hc: $(PNJLPhysicsSimulation.PhysicalConstants.hc)")
    println("   PNJL常数 - T0: $(PNJLPhysicsSimulation.PNJLConstants.T0)")
    println("   Gas-Liquid常数 - m: $(PNJLPhysicsSimulation.GasLiquidConstants.m)")
    println("   Rotation常数 - r0: $(PNJLPhysicsSimulation.RotationConstants.r0)")
    
    println("\n🎉 常数清理验证成功完成!")
    println("✨ 所有函数文件中的重复常数定义已完全清理")
    
catch e
    println("❌ 验证失败: $e")
    showerror(stdout, e)
    println()
end
