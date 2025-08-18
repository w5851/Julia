"""
积分接口模块的单元测试

测试IntegrationInterface模块的基本功能，包括：
- 积分方法和网格的创建
- 基础积分函数的正确性
- 物理专用积分函数的数值精度
- 错误处理和边界情况
"""

using Test
using PNJLPhysicsSimulation
using PNJLPhysicsSimulation.IntegrationInterface

@testset "IntegrationInterface Tests" begin
    
    @testset "Grid Creation" begin
        # 测试动量网格创建
        @testset "MomentumGrid Creation" begin
            grid = create_momentum_grid(64, 20.0)
            
            @test grid isa MomentumGrid
            @test length(grid.nodes) == 64
            @test length(grid.weights) == 64
            @test grid.cutoff == 20.0
            @test grid.domain == (0.0, 20.0)
            
            # 检查节点是否在正确范围内
            @test all(0.0 ≤ node ≤ 20.0 for node in grid.nodes)
            @test all(weight > 0.0 for weight in grid.weights)
        end
        
        # 测试角度网格创建
        @testset "AngleGrid Creation" begin
            grid = create_angle_grid(16)
            
            @test grid isa AngleGrid
            @test length(grid.nodes) == 16
            @test length(grid.weights) == 16
            @test grid.domain == (-1.0, 1.0)
            
            # 检查节点范围
            @test all(-1.0 ≤ node ≤ 1.0 for node in grid.nodes)
            @test all(weight > 0.0 for weight in grid.weights)
        end
        
        # 测试乘积网格
        @testset "ProductGrid Creation" begin
            p_grid = create_momentum_grid(32, 10.0)
            t_grid = create_angle_grid(8)
            product_grid = create_product_grid(p_grid, t_grid)
            
            @test product_grid isa ProductGrid
            @test length(product_grid.grids) == 2
            @test product_grid.grids[1] isa MomentumGrid
            @test product_grid.grids[2] isa AngleGrid
        end
    end
    
    @testset "Integration Methods" begin
        # 测试基础积分功能
        @testset "Basic Integration" begin
            grid = create_momentum_grid(64, 1.0)
            method = GaussLegendreIntegration()
            
            # 测试简单函数积分
            # ∫₀¹ x dx = 1/2
            result = integrate(method, grid, x -> x)
            @test isapprox(result, 0.5, rtol=1e-10)
            
            # ∫₀¹ x² dx = 1/3
            result = integrate(method, grid, x -> x^2)
            @test isapprox(result, 1/3, rtol=1e-10)
            
            # ∫₀¹ exp(x) dx = e - 1
            result = integrate(method, grid, x -> exp(x))
            @test isapprox(result, exp(1.0) - 1.0, rtol=1e-10)
        end
        
        # 测试二维积分
        @testset "2D Integration" begin
            p_grid = create_momentum_grid(32, 1.0)
            t_grid = create_angle_grid(16)
            method = GaussLegendreIntegration()
            
            # ∫₀¹ ∫₋₁¹ xy dx dy = 0 (奇函数)
            result = integrate_2d(method, p_grid, t_grid, (x, y) -> x * y)
            @test isapprox(result, 0.0, atol=1e-12)
            
            # ∫₀¹ ∫₋₁¹ x² dx dy = 2/3
            result = integrate_2d(method, p_grid, t_grid, (x, y) -> x^2)
            @test isapprox(result, 2.0/3.0, rtol=1e-10)
        end
        
        # 测试向量化积分
        @testset "Vectorized Integration" begin
            grid = create_momentum_grid(64, 1.0)
            method = GaussLegendreIntegration()
            
            functions = [x -> x, x -> x^2, x -> x^3]
            expected = [0.5, 1/3, 0.25]
            
            results = integrate_vectorized(method, grid, functions)
            
            @test length(results) == 3
            for (result, expect) in zip(results, expected)
                @test isapprox(result, expect, rtol=1e-10)
            end
        end
    end
    
    @testset "Physics Integration Functions" begin
        # 测试物理积分函数的基本功能
        @testset "Thermal Integral Structure" begin
            # 设置测试参数
            masses = [0.1, 0.1, 0.5]  # 简化的夸克质量
            mu = [0.3, 0.3, 0.3]      # 化学势
            T = 0.15                   # 温度
            Phi1, Phi2 = 0.5, 0.5     # Polyakov loop参数
            
            grid = create_momentum_grid(64, 10.0)
            method = GaussLegendreIntegration()
            
            # 测试函数是否正常执行（不测试具体数值，因为物理模型复杂）
            result = omega_thermal_integral(masses, mu, T, Phi1, Phi2, grid, method)
            
            @test result isa Float64
            @test isfinite(result)
            
            # 测试温度依赖性：高温时积分应该更大
            T_high = 0.3
            result_high = omega_thermal_integral(masses, mu, T_high, Phi1, Phi2, grid, method)
            @test abs(result_high) > abs(result)  # 高温贡献应该更大
        end
        
        @testset "Vacuum Energy Integral" begin
            masses = [0.1, 0.1, 0.5]
            grid = create_momentum_grid(64, 10.0)
            method = GaussLegendreIntegration()
            
            result = vacuum_energy_integral(masses, grid, method)
            
            @test result isa Float64
            @test isfinite(result)
            # 真空能量的符号取决于计算方式和规范化，这里只检查有限性
            
            # 测试质量依赖性：更大的质量应该导致不同的真空能量
            heavy_masses = [1.0, 1.0, 1.0]
            result_heavy = vacuum_energy_integral(heavy_masses, grid, method)
            @test result_heavy ≠ result
        end
    end
    
    @testset "Error Handling" begin
        # 测试异常情况的处理
        @testset "Invalid Grid Parameters" begin
            # 测试不匹配的节点和权重长度
            @test_throws ArgumentError MomentumGrid([1.0, 2.0], [1.0], (0.0, 2.0), 2.0)
            @test_throws ArgumentError AngleGrid([1.0, 2.0], [1.0])
        end
        
        @testset "Problematic Integrands" begin
            grid = create_momentum_grid(32, 1.0)
            method = GaussLegendreIntegration()
            
            # 测试发散积分的处理（应该不会崩溃）
            result = integrate(method, grid, x -> x > 0.5 ? Inf : 0.0)
            @test isfinite(result) || isnan(result)  # 应该有合理的处理
        end
    end
    
    @testset "Consistency with Legacy Functions" begin
        # 这里为未来的一致性测试预留空间
        # 当实施阶段二时，将添加与原有函数的对比测试
        @test true  # 占位符测试
    end
end

# 性能基准测试（可选运行）
function run_benchmark_tests()
    using BenchmarkTools
    
    println("Integration Interface Performance Benchmarks")
    println("=" ^ 50)
    
    grid = create_momentum_grid(128, 20.0)
    method = GaussLegendreIntegration()
    
    # 基础积分性能
    println("Basic integration benchmark:")
    @btime integrate($method, $grid, x -> x^2 * exp(-x))
    
    # 物理积分性能
    masses = [0.1, 0.1, 0.5]
    mu = [0.3, 0.3, 0.3]
    T = 0.15
    Phi1, Phi2 = 0.5, 0.5
    
    println("Physics integration benchmark:")
    @btime omega_thermal_integral($masses, $mu, $T, $Phi1, $Phi2, $grid, $method)
    
    println("Vacuum energy integration benchmark:")
    @btime vacuum_energy_integral($masses, $grid, $method)
end

# 运行基准测试的提示
println("""
基础测试已完成。如需运行性能基准测试，请执行：
    include("test/test_integration_interface.jl")
    run_benchmark_tests()
""")
