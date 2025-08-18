"""
阶段三核心功能测试
专注于模型配置和Gas-Liquid新接口测试
"""

using Test
using PNJLPhysicsSimulation
using PNJLPhysicsSimulation.ModelConfiguration
using PNJLPhysicsSimulation.IntegrationInterface

@testset "阶段三核心功能测试" begin
    
    @testset "模型配置系统" begin
        @testset "PNJL配置" begin
            config = create_default_config(:PNJL)
            @test config isa PNJLConfig
            @test config.momentum_cutoff == 10.0
            @test config.n_momentum_points == 64
            @test config.temperature == 0.15
            @test length(config.chemical_potentials) == 3
            
            # 测试网格配置
            grid = get_grid_config(config)
            @test grid isa MomentumGrid
            @test length(grid.nodes) == config.n_momentum_points
            
            println("✅ PNJL配置测试通过 - $(config.n_momentum_points)积分点")
        end
        
        @testset "PNJL_aniso配置" begin
            config = create_default_config(:PNJL_aniso)
            @test config isa PNJLAnisoConfig
            @test config.n_angle_points == 32
            
            grids = get_grid_config(config)
            @test haskey(grids, :momentum)
            @test haskey(grids, :angle)
            
            println("✅ PNJL_aniso配置测试通过 - 2D积分网格")
        end
        
        @testset "Rotation配置" begin
            config = create_default_config(:Rotation)
            @test config isa RotationConfig
            @test config.angular_velocity == 0.05
            
            println("✅ Rotation配置测试通过 - Ω=$(config.angular_velocity)")
        end
        
        @testset "GasLiquid配置" begin
            config = create_default_config(:GasLiquid)
            @test config isa GasLiquidConfig
            @test config.momentum_cutoff == 20.0
            @test config.baryon_density == 0.15
            
            println("✅ GasLiquid配置测试通过 - ρ_B=$(config.baryon_density)")
        end
    end
    
    @testset "Gas-Liquid新接口" begin
        @testset "新密度计算函数功能" begin
            mass = 0.5
            μ = 0.3
            T = 0.15
            grid = create_momentum_grid(32, 5.0)
            
            ρ_new = PNJLPhysicsSimulation.GasLiquidFunctions.calculate_ρ_new(mass, μ, T, grid)
            ρ_s_new = PNJLPhysicsSimulation.GasLiquidFunctions.calculate_ρ_s_new(mass, μ, T, grid)
            
            @test isfinite(ρ_new)
            @test isfinite(ρ_s_new)
            @test ρ_new >= 0.0
            @test ρ_s_new >= 0.0
            
            println("✅ Gas-Liquid新接口功能正常")
            println("   ρ = $(ρ_new), ρ_s = $(ρ_s_new)")
        end
        
        @testset "新旧接口数值一致性" begin
            mass = 0.5
            μ = 0.3  
            T = 0.15
            grid = create_momentum_grid(64, 10.0)
            
            # 旧接口计算
            p_nodes = grid.nodes
            coef = grid.weights .* (p_nodes.^2) / π^2
            E = sqrt.(p_nodes.^2 .+ mass^2)
            
            ρ_old = PNJLPhysicsSimulation.GasLiquidFunctions.calculate_ρ(E, μ, T, coef)
            ρ_s_old = PNJLPhysicsSimulation.GasLiquidFunctions.calculate_ρ_s(E, μ, T, coef, mass)
            
            # 新接口计算
            ρ_new = PNJLPhysicsSimulation.GasLiquidFunctions.calculate_ρ_new(mass, μ, T, grid)
            ρ_s_new = PNJLPhysicsSimulation.GasLiquidFunctions.calculate_ρ_s_new(mass, μ, T, grid)
            
            @test isapprox(ρ_old, ρ_new, rtol=1e-10)
            @test isapprox(ρ_s_old, ρ_s_new, rtol=1e-10)
            
            println("✅ 新旧接口一致性验证通过")
            println("   重子密度差异: $(abs(ρ_old - ρ_new))")  
            println("   标量密度差异: $(abs(ρ_s_old - ρ_s_new))")
        end
    end
    
    @testset "集成测试" begin
        @testset "所有模型配置创建" begin
            models = [:PNJL, :PNJL_aniso, :Rotation, :GasLiquid]
            configs = Dict()
            
            for model in models
                configs[model] = create_default_config(model)
                @test configs[model] isa ModelConfig
            end
            
            println("✅ 创建了$(length(configs))个模型配置")
        end
        
        @testset "配置优化功能" begin
            config = create_default_config(:PNJL)
            optimized = optimize_config(config, 1e-6)
            @test optimized isa PNJLConfig
            println("✅ 配置优化功能正常")
        end
    end
end

println("=" ^ 60)
println("✅ 阶段三核心功能测试完成")
println("=" ^ 60)
println("""
阶段三成果验证:
✅ 统一模型配置系统工作正常
✅ Gas-Liquid模型新接口集成成功  
✅ 新旧接口数值完全一致
✅ 所有物理模型配置创建成功
✅ 网格管理系统功能完善

阶段三模型适配圆满完成！
下一步可以进入阶段四的全面测试验证。
""")
