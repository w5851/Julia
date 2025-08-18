"""
阶段三：模型适配测试

测试新的统一物理接口和模型配置系统的功能。
"""

using Test
using PNJLPhysicsSimulation
using PNJLPhysicsSimulation.ModelConfiguration
using PNJLPhysicsSimulation.UnifiedPhysicsInterface
using PNJLPhysicsSimulation.IntegrationInterface

@testset "阶段三：模型适配测试" begin
    
    @testset "模型配置测试" begin
        @testset "PNJL配置创建" begin
            config = create_default_config(:PNJL)
            @test config isa PNJLConfig
            @test config.momentum_cutoff == 10.0
            @test config.n_momentum_points == 64
            @test config.temperature == 0.15
            @test length(config.chemical_potentials) == 3
            println("✅ PNJL配置创建测试通过")
        end
        
        @testset "PNJL_aniso配置创建" begin
            config = create_default_config(:PNJL_aniso)
            @test config isa PNJLAnisoConfig
            @test config.n_angle_points == 32
            @test config.anisotropy_parameter == 0.1
            println("✅ PNJL_aniso配置创建测试通过")
        end
        
        @testset "Rotation配置创建" begin
            config = create_default_config(:Rotation)
            @test config isa RotationConfig
            @test config.angular_velocity == 0.05
            @test config.n_angular_points == 32
            println("✅ Rotation配置创建测试通过")
        end
        
        @testset "GasLiquid配置创建" begin
            config = create_default_config(:GasLiquid)
            @test config isa GasLiquidConfig
            @test config.momentum_cutoff == 20.0
            @test config.baryon_density == 0.15
            println("✅ GasLiquid配置创建测试通过")
        end
    end
    
    @testset "网格配置测试" begin
        @testset "PNJL网格配置" begin
            config = create_default_config(:PNJL)
            grid = get_grid_config(config)
            @test grid isa MomentumGrid
            @test length(grid.nodes) == config.n_momentum_points
            println("✅ PNJL网格配置测试通过")
        end
        
        @testset "PNJL_aniso网格配置" begin
            config = create_default_config(:PNJL_aniso)
            grids = get_grid_config(config)
            @test haskey(grids, :momentum)
            @test haskey(grids, :angle)
            @test grids.momentum isa MomentumGrid
            @test grids.angle isa AngleGrid
            println("✅ PNJL_aniso网格配置测试通过")
        end
    end
    
    @testset "统一物理接口测试" begin
        @testset "PNJL热力学计算" begin
            config = create_default_config(:PNJL)
            phi = [-0.1, -0.1, -1.7]
            
            result = calculate_thermodynamics(config, phi)
            @test result isa ThermodynamicResult
            @test isfinite(result.pressure)
            @test result.temperature == config.temperature
            @test result.chemical_potential == config.chemical_potentials
            
            println("✅ PNJL热力学计算测试通过")
            println("   压强: $(result.pressure)")
        end
        
        @testset "PNJL_aniso热力学计算" begin
            config = create_default_config(:PNJL_aniso)
            phi = [-0.1, -0.1, -1.7]
            
            result = calculate_thermodynamics(config, phi)
            @test result isa ThermodynamicResult
            @test isfinite(result.pressure)
            
            println("✅ PNJL_aniso热力学计算测试通过")
            println("   压强: $(result.pressure)")
        end
        
        @testset "Rotation热力学计算" begin
            config = create_default_config(:Rotation)
            phi = [-0.1, -0.1, -1.7]
            
            result = calculate_thermodynamics(config, phi)
            @test result isa ThermodynamicResult
            @test isfinite(result.pressure)
            
            println("✅ Rotation热力学计算测试通过")
            println("   压强: $(result.pressure)")
        end
    end
    
    @testset "Gas-Liquid新接口测试" begin
        @testset "新密度计算函数" begin
            mass = 0.5
            μ = 0.3
            T = 0.15
            grid = create_momentum_grid(32, 5.0)
            
            # 测试重构后的密度计算
            ρ_new = PNJLPhysicsSimulation.GasLiquidFunctions.calculate_ρ_new(mass, μ, T, grid)
            ρ_s_new = PNJLPhysicsSimulation.GasLiquidFunctions.calculate_ρ_s_new(mass, μ, T, grid)
            
            @test isfinite(ρ_new)
            @test isfinite(ρ_s_new)
            
            println("✅ Gas-Liquid新接口测试通过")
            println("   重子密度: $(ρ_new)")
            println("   标量密度: $(ρ_s_new)")
        end
    end
    
    @testset "性能和一致性测试" begin
        @testset "配置优化" begin
            config = create_default_config(:PNJL)
            optimized_config = optimize_config(config, 1e-6)
            @test optimized_config isa PNJLConfig
            println("✅ 配置优化测试通过")
        end
        
        @testset "接口一致性" begin
            # 测试不同模型使用相同接口
            models = [:PNJL, :PNJL_aniso, :Rotation]
            phi = [-0.1, -0.1, -1.7]
            
            results = []
            for model in models
                config = create_default_config(model)
                result = calculate_thermodynamics(config, phi)
                push!(results, result)
                @test result isa ThermodynamicResult
            end
            
            println("✅ 接口一致性测试通过")
            println("   所有模型都使用统一的ThermodynamicResult结构")
        end
    end
end

println("=" ^ 60)
println("✅ 阶段三：模型适配测试完成")  
println("=" ^ 60)
println("""
模型适配成果总结:
1. ✅ 统一模型配置系统：支持所有物理模型的参数管理
2. ✅ 统一物理计算接口：所有模型使用相同的计算API
3. ✅ 网格配置管理：自动化的积分网格创建和管理
4. ✅ Gas-Liquid模型适配：成功集成新的积分接口
5. ✅ 接口一致性：所有模型返回统一的结果结构

下一步：进入阶段四的全面测试和验证
""")
