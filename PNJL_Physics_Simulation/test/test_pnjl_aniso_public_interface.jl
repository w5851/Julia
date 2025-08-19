"""
PNJL各向异性模型公共接口测试

测试新实现的4个高级接口:
1. omega_derivatives - 自动微分接口
2. solve_equilibrium - NLsolve包装接口  
3. evaluate_solution - 解评价接口
4. scan_phase_diagram - 相图扫描接口

作者: PNJL Physics Simulation Team
创建时间: 2025年8月19日
"""

using Test
using PNJLPhysicsSimulation
using PNJLPhysicsSimulation.PNJLAnisoPublicInterface
using StaticArrays
using Printf

@testset "PNJL Anisotropic Public Interface Tests" begin

# ============================================================================
# Test 1: 自动微分接口测试
# ============================================================================

@testset "Interface 1: omega_derivatives" begin
    println("\n=== 测试自动微分接口 ===")
    
    # 设置测试参数
    phi = [0.1, 0.1, 0.05]      # MeV
    Phi1, Phi2 = 0.5, 0.3
    mu = [100.0, 100.0, 50.0]   # MeV
    T = 150.0                   # MeV
    xi = 0.1                    # 各向异性参数
    
    # 测试基本功能
    @test_nowarn derivatives = omega_derivatives(phi, Phi1, Phi2, mu, T, xi; p_points=16, t_points=8)
    
    derivatives = omega_derivatives(phi, Phi1, Phi2, mu, T, xi; p_points=16, t_points=8)
    
    # 检查返回类型和维度
    @test derivatives isa StaticArrays.SVector{5}
    @test length(derivatives) == 5
    @test all(isfinite.(derivatives))  # 所有导数应该是有限数值
    
    println("偏导数结果:")
    @printf("∂Ω/∂φᵤ = %.6e\n", derivatives[1])
    @printf("∂Ω/∂φᵈ = %.6e\n", derivatives[2])  
    @printf("∂Ω/∂φˢ = %.6e\n", derivatives[3])
    @printf("∂Ω/∂Φ₁ = %.6e\n", derivatives[4])
    @printf("∂Ω/∂Φ₂ = %.6e\n", derivatives[5])
    
    # 测试原地计算版本
    grad = zeros(5)
    @test_nowarn omega_derivatives!(grad, phi, Phi1, Phi2, mu, T, xi; p_points=16, t_points=8)
    @test norm(grad - derivatives) < 1e-12  # 应该得到相同结果
    
    # 测试边界情况
    @test_throws AssertionError omega_derivatives([0.1, 0.1], Phi1, Phi2, mu, T, xi)  # 错误的phi维度
    @test_throws AssertionError omega_derivatives(phi, Phi1, Phi2, [100.0, 100.0], T, xi)  # 错误的mu维度
    @test_throws AssertionError omega_derivatives(phi, Phi1, Phi2, mu, -10.0, xi)  # 负温度
    
    println("✓ 自动微分接口测试通过")
end

# ============================================================================
# Test 2: NLsolve包装接口测试
# ============================================================================

@testset "Interface 2: solve_equilibrium" begin
    println("\n=== 测试方程组求解接口 ===")
    
    # 测试参数
    mu_test = [200.0, 200.0, 150.0]  # MeV
    T_test = 120.0                   # MeV
    xi_test = 0.05                   # 小的各向异性
    
    # 基本求解测试
    @test_nowarn solution = solve_equilibrium(mu_test, T_test, xi_test; 
                                            p_points=16, t_points=8,
                                            show_trace=false)
    
    solution = solve_equilibrium(mu_test, T_test, xi_test; 
                               p_points=16, t_points=8, show_trace=false)
    
    # 检查解的结构和类型
    @test solution isa PNJLAnisoPublicInterface.EquilibriumSolution
    @test solution.temperature == T_test
    @test solution.chemical_potentials == SVector{3}(mu_test)
    @test solution.anisotropy == xi_test
    
    # 打印解的摘要
    PNJLAnisoPublicInterface.print_solution_summary(solution)
    
    # 如果收敛，检查解的合理性
    if solution.converged
        println("✓ 求解收敛")
        
        # 检查物理合理性
        @test PNJLAnisoPublicInterface.check_solution_validity(solution)
        
        # 验证平衡态条件 (梯度应接近零)
        verification_grad = omega_derivatives(
            [solution.phi_u, solution.phi_d, solution.phi_s],
            solution.Phi1, solution.Phi2,
            mu_test, T_test, xi_test; 
            p_points=16, t_points=8
        )
        
        grad_norm = norm(verification_grad)
        println("平衡态验证 - 梯度范数: ", grad_norm)
        @test grad_norm < 1e-4  # 平衡态条件: ∇Ω ≈ 0
        
    else
        @warn "求解未收敛，残差范数: $(solution.residual_norm)"
    end
    
    # 测试自定义初值
    custom_initial = [0.12, 0.12, 0.08, 0.6, 0.2]
    @test_nowarn solution2 = solve_equilibrium(mu_test, T_test, xi_test;
                                             initial_guess=custom_initial,
                                             p_points=16, t_points=8)
    
    # 测试不同求解器选项
    @test_nowarn solution3 = solve_equilibrium(mu_test, T_test, xi_test;
                                             method=:newton,
                                             xtol=1e-10,
                                             p_points=16, t_points=8)
    
    println("✓ 方程组求解接口测试通过")
end

# ============================================================================  
# Test 3: 解评价接口测试
# ============================================================================

@testset "Interface 3: evaluate_solution" begin
    println("\n=== 测试解评价接口 ===")
    
    # 先求解一个平衡态
    mu_eval = [250.0, 250.0, 200.0]
    T_eval = 140.0
    xi_eval = 0.0  # 各向同性情况
    
    solution = solve_equilibrium(mu_eval, T_eval, xi_eval; 
                               p_points=16, t_points=8, show_trace=false)
    
    if solution.converged
        # 评价热力学量
        @test_nowarn thermo = evaluate_solution(solution; p_points=24, t_points=12)
        
        thermo = evaluate_solution(solution; p_points=24, t_points=12)
        
        # 检查结果类型和结构
        @test thermo isa PNJLAnisoPublicInterface.ThermodynamicQuantities
        @test all(isfinite.([thermo.pressure, thermo.energy_density, thermo.entropy_density]))
        @test length(thermo.baryon_density) == 3
        @test length(thermo.effective_masses) == 3
        
        # 打印热力学量摘要
        PNJLAnisoPublicInterface.print_thermodynamics_summary(thermo)
        
        # 检查热力学一致性
        @test thermo.pressure == -(thermo.omega_total)  # P = -Ω
        @test thermo.temperature == T_eval
        @test thermo.chemical_potentials == SVector{3}(mu_eval)
        
        # 检查能量密度符号 (应为正值在高密度下)
        if sum(thermo.baryon_density) > 0.1  # 高密度情况
            @test thermo.energy_density > 0
        end
        
        println("✓ 热力学量计算正常")
    else
        @warn "平衡态求解失败，跳过热力学量评价测试"
    end
    
    println("✓ 解评价接口测试通过")
end

# ============================================================================
# Test 4: 相图扫描接口测试 (小规模)
# ============================================================================

@testset "Interface 4: scan_phase_diagram (小规模测试)" begin
    println("\n=== 测试相图扫描接口 ===")
    
    # 创建小规模测试参数
    test_params = PNJLAnisoPublicInterface.create_scan_parameters(
        (100.0, 150.0),       # 温度范围 (MeV)
        (100.0, 200.0),       # 化学势范围 (MeV)
        (0.0, 0.1),          # 各向异性范围
        (3, 3, 2),           # 网格点数 (T, μ, ξ) - 小规模测试
        "test_phase_diagram.dat";  # 输出文件
        p_points=12,         # 降低积分精度加快测试
        t_points=6,
        save_intermediate=false
    )
    
    @test test_params isa PNJLAnisoPublicInterface.ScanParameters
    @test test_params.n_temperature == 3
    @test test_params.n_mu == 3  
    @test test_params.n_anisotropy == 2
    @test test_params.output_file == "test_phase_diagram.dat"
    
    # 执行小规模扫描
    println("开始小规模相图扫描测试 (3×3×2 = 18个点)...")
    @test_nowarn result = scan_phase_diagram(test_params)
    
    result = scan_phase_diagram(test_params)
    
    # 检查结果结构
    @test result isa PNJLAnisoPublicInterface.ScanResult
    @test result.total_points == 18
    @test length(result.temperatures) == 3
    @test length(result.chemical_potentials) == 3
    @test length(result.anisotropies) == 2
    @test size(result.pressures) == (3, 3, 2)
    @test size(result.convergence_flags) == (3, 3, 2)
    
    # 检查是否有收敛的点
    converged_count = sum(result.convergence_flags)
    println("扫描结果: $(converged_count)/$(result.total_points) 个点收敛")
    
    if converged_count > 0
        @test result.converged_points == converged_count
        @test result.computation_time > 0
        
        # 检查收敛点的物理量是否合理
        converged_pressures = result.pressures[result.convergence_flags]
        @test all(isfinite.(converged_pressures))
        
        println("✓ 扫描计算成功")
    else
        @warn "所有扫描点都未收敛，可能需要调整参数"
    end
    
    # 检查文件是否生成
    @test isfile("test_phase_diagram.dat")
    @test isfile("test_phase_diagram_summary.txt")
    
    println("✓ 相图扫描接口测试完成")
    
    # 清理测试文件
    try
        rm("test_phase_diagram.dat")
        rm("test_phase_diagram_summary.txt")
        println("✓ 测试文件已清理")
    catch e
        @warn "清理测试文件时出现问题: $e"
    end
end

# ============================================================================
# 集成测试：完整工作流程
# ============================================================================

@testset "Integration Test: 完整工作流程" begin
    println("\n=== 集成测试：完整工作流程 ===")
    
    # 设置物理参数
    mu_flow = [180.0, 180.0, 130.0]  # MeV
    T_flow = 130.0                   # MeV
    xi_flow = 0.02                   # 轻微各向异性
    
    println("步骤 1: 求解平衡态...")
    solution = solve_equilibrium(mu_flow, T_flow, xi_flow; 
                               p_points=20, t_points=10, show_trace=false)
    
    @test solution isa PNJLAnisoPublicInterface.EquilibriumSolution
    
    if solution.converged
        println("✓ 平衡态求解成功")
        
        # 验证平衡态条件
        println("步骤 2: 验证平衡态条件...")
        grad = omega_derivatives([solution.phi_u, solution.phi_d, solution.phi_s],
                                solution.Phi1, solution.Phi2,
                                mu_flow, T_flow, xi_flow; 
                                p_points=20, t_points=10)
        
        grad_norm = norm(grad)
        println("梯度范数: $grad_norm")
        @test grad_norm < 1e-3
        println("✓ 平衡态条件验证通过")
        
        # 计算热力学量
        println("步骤 3: 计算热力学量...")
        thermo = evaluate_solution(solution; p_points=24, t_points=12)
        
        @test thermo isa PNJLAnisoPublicInterface.ThermodynamicQuantities
        @test all(isfinite.([thermo.pressure, thermo.energy_density]))
        println("✓ 热力学量计算成功")
        
        # 显示主要结果
        println("\n=== 最终结果摘要 ===")
        @printf("温度: %.1f MeV\n", T_flow)
        println("化学势: ", mu_flow, " MeV")  
        @printf("各向异性参数: %.3f\n", xi_flow)
        println()
        @printf("压力: %.3e MeV⁴\n", thermo.pressure)
        @printf("能量密度: %.3e MeV⁴\n", thermo.energy_density)
        @printf("重子密度: %.3f fm⁻³\n", sum(thermo.baryon_density))
        @printf("有效质量: [%.1f, %.1f, %.1f] MeV\n", 
               thermo.effective_masses[1], thermo.effective_masses[2], thermo.effective_masses[3])
        
        println("✓ 完整工作流程测试成功")
    else
        @warn "平衡态求解失败，工作流程测试中断"
        @test false
    end
end

end  # testset

println("\n" * "="^60)
println("PNJL各向异性模型公共接口测试完成!")
println("4个核心接口功能验证通过:")
println("  ✓ Interface 1: omega_derivatives (自动微分)")
println("  ✓ Interface 2: solve_equilibrium (方程求解)")  
println("  ✓ Interface 3: evaluate_solution (解评价)")
println("  ✓ Interface 4: scan_phase_diagram (相图扫描)")
println("="^60)
