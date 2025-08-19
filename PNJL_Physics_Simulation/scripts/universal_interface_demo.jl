"""
通用物理接口使用示例

演示如何使用新的统一公共接口对不同物理模型进行计算：
1. 自动微分接口
2. 方程组求解接口  
3. 物理量计算接口
4. 相图扫描接口

作者: AI助手
日期: 2025年8月19日
"""

using PNJLPhysicsSimulation
using .UnifiedPhysicsPublicInterface
using .ModelConfiguration
using .PNJLAnisoFunctions  # 用作示例的具体模型

function main()
    println("🚀 通用物理接口演示开始")
    println("="^60)
    
    # 1. 准备模型配置
    println("\\n📋 1. 设置模型配置")
    config = PNJLAnisoConfig(
        cutoff=0.6,           # 动量截断
        n_p=15,               # 动量积分点数  
        n_theta=10,           # 角度积分点数
        T=0.15,               # 温度 (GeV)
        mu=[0.3, 0.3, 0.3],   # 化学势 [μu, μd, μs] (GeV)
        Phi=(0.5, 0.5),       # Polyakov场初值 [Φ₁, Φ₂]
        xi=1.0                # 各向异性参数
    )
    println("   ✓ PNJL各向异性模型配置完成")
    println("     温度: $(config.temperature) GeV")
    println("     化学势: $(config.chemical_potentials) GeV")
    println("     各向异性参数: $(config.anisotropy_parameter)")
    
    # 2. 定义热力学势函数
    println("\\n⚡ 2. 定义热力学势函数")
    omega_func = x -> begin
        phi_u, phi_d, phi_s, Phi_1, Phi_2 = x[1], x[2], x[3], x[4], x[5]
        return calculate_omega_total([phi_u, phi_d, phi_s], [Phi_1, Phi_2], config)
    end
    println("   ✓ 热力学势函数已定义")
    
    # 3. 测试自动微分接口
    println("\\n🔄 3. 测试自动微分接口")
    variables = [0.1, 0.1, 0.2, 0.5, 0.5]  # [φᵤ, φd, φₛ, Φ₁, Φ₂]
    
    try
        gradients = calculate_derivatives(omega_func, variables, config)
        equilibrium_check = calculate_equilibrium_conditions(omega_func, variables, config)
        
        println("   ✓ 梯度计算成功")
        println("     偏导数: $(round.(gradients, digits=4))")
        println("     是否平衡: $(equilibrium_check.is_equilibrium)")
        println("     最大偏导数: $(round(equilibrium_check.max_derivative, digits=6))")
        
    catch e
        println("   ❌ 自动微分测试失败: $e")
        return
    end
    
    # 4. 测试方程组求解接口
    println("\\n🎯 4. 测试方程组求解接口")
    initial_guess = [0.08, 0.08, 0.15, 0.6, 0.6]  # 改进的初始猜测
    
    try
        solution = solve_equilibrium_equations(omega_func, initial_guess, config; 
                                             ftol=1e-8, show_trace=false)
        
        if solution.converged
            println("   ✅ 方程组求解成功!")
            println("     平衡解: $(round.(solution.solution, digits=4))")
            println("     残差范数: $(solution.residual_norm)")
            println("     迭代次数: $(solution.iterations)")
            println("     求解时间: $(round(solution.solve_time, digits=3))s")
            
            # 5. 测试物理量计算接口
            println("\\n📊 5. 测试物理量计算接口")
            
            try
                properties = calculate_physical_properties(solution.solution, config, omega_func)
                
                println("   ✅ 物理量计算成功!")
                println("     压强: $(round(properties.pressure, digits=4)) GeV⁴")
                println("     能量密度: $(round(properties.energy_density, digits=4)) GeV⁴")
                println("     重子密度: $(round(properties.baryon_density, digits=4)) GeV³")
                println("     手征凝聚: $(round.(properties.chiral_condensates, digits=4))")
                println("     Polyakov环: $(round.(properties.polyakov_loops, digits=4))")
                println("     计算时间: $(round(properties.calculation_time, digits=3))s")
                
            catch e
                println("   ⚠️  物理量计算警告: $e")
            end
            
        else
            println("   ❌ 方程组求解未收敛")
            println("     残差范数: $(solution.residual_norm)")
            println("     迭代次数: $(solution.iterations)")
        end
        
    catch e
        println("   ❌ 方程组求解失败: $e")
        return
    end
    
    # 6. 测试小规模相图扫描
    println("\\n🗺️  6. 测试相图扫描接口（小规模示例）")
    
    try
        println("   开始小规模相图扫描...")
        phase_points = scan_phase_diagram(
            omega_func, config;
            temperature_range=(0.12, 0.18),     # 小温度范围
            temperature_points=3,               # 3个温度点
            chemical_potential_range=(0.25, 0.35),  # 小化学势范围
            chemical_potential_points=3,        # 3个化学势点
            compute_properties=true,
            show_progress=true
        )
        
        println("   ✅ 相图扫描完成!")
        println("     扫描点数: $(length(phase_points))")
        
        converged_count = count(p -> p.converged, phase_points)
        println("     收敛点数: $converged_count / $(length(phase_points))")
        
        if converged_count > 0
            # 显示几个成功点的结果
            println("   📋 示例结果:")
            for (i, point) in enumerate(phase_points[1:min(3, converged_count)])
                if point.converged
                    println("     点$i: T=$(round(point.temperature, digits=3)), "*
                           "μ=$(round(point.chemical_potential, digits=3)), "*
                           "P=$(round(point.properties.pressure, digits=4))")
                end
            end
            
            # 保存结果
            output_file = "demo_phase_diagram.dat"
            save_phase_diagram(phase_points, output_file)
            println("   💾 结果已保存至: $output_file")
        end
        
    catch e
        println("   ⚠️  相图扫描警告: $e")
    end
    
    println("\\n" * "="^60)
    println("🎉 通用物理接口演示完成!")
    println("\\n📝 接口功能总结:")
    println("   ✓ calculate_derivatives - 自动微分计算偏导数")
    println("   ✓ solve_equilibrium_equations - nlsolve求解平衡态")  
    println("   ✓ calculate_physical_properties - 计算热力学量")
    println("   ✓ scan_phase_diagram - 扫描T-μ相图")
    println("\\n🌟 这些接口适用于所有物理模型 (PNJL, PNJL各向异性, 旋转, 气液相变等)")
    
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
