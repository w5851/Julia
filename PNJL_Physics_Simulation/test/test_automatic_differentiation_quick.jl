using Test

# 加载被测模块（使用 @__DIR__ 保证相对路径正确）
include(joinpath(@__DIR__, "..", "src", "core", "autodiff_interface.jl"))

@testset "AutodiffInterface quick tests" begin
    # 测试 compute_gradient
    vars = [1.2, -0.5]
    cfg = (T = 1.5, mu = [0.3, 0.4])
    omega = (x, cfg) -> x[1]^2 + 3.0 * x[2] + cfg.T * x[1] + sum(cfg.mu)

    grad = AutodiffInterface.compute_gradient(omega, vars, cfg)
    @test isapprox(grad[1], 2.0 * vars[1] + cfg.T; atol=1e-8)
    @test isapprox(grad[2], 3.0; atol=1e-8)

    # 测试 compute_temperature_derivative
    omegaT = (x, T, cfg) -> x[1]^2 + T * x[1] + sum(cfg.mu)
    dTd = AutodiffInterface.compute_temperature_derivative(omegaT, vars, cfg.T, cfg)
    @test isapprox(dTd, vars[1]; atol=1e-8)

    # 测试 compute_chemical_potential_derivatives
    mu = [0.7, -0.2]
    cfg2 = (extra = 0.0,)
    omegaMu = (x, mu_vec, cfg) -> x[1]^2 + mu_vec[1] * x[1] + mu_vec[2] * x[2] + cfg.extra
    dOm_dmu = AutodiffInterface.compute_chemical_potential_derivatives(omegaMu, vars, mu, cfg2)
    @test isapprox(dOm_dmu[1], vars[1]; atol=1e-8)
    @test isapprox(dOm_dmu[2], vars[2]; atol=1e-8)
end
