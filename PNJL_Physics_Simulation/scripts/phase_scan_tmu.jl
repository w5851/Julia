# Phase-scan wrapper: t-μ 扫描（mu 为标量）
# 运行示例： julia --project=. scripts/phase_scan_tmu.jl

# 加载核心接口模块和物理模块
include(joinpath(@__DIR__, "..", "src", "core", "autodiff_interface.jl"))
using .AutodiffInterface
include(joinpath(@__DIR__, "..", "src", "core", "phase_scan_interface.jl"))
using .PhaseScanInterface: ParameterAxis, scan_phase_space

include(joinpath(@__DIR__, "..", "src", "models", "pnjl_aniso", "functions.jl"))
using .PNJLAnisoFunctions: create_aniso_grids, pressure_wrapper_modern, calculate_core_modern

using NLsolve, Printf

# 默认网格（粗）
p_grid, t_grid = create_aniso_grids(8, 8)
xi_default = 0.0

# 初始猜测函数（5 个自由度：phi_u,phi_d,phi_s,Phi1,Phi2）
initial_guess = params -> zeros(5)

# solver wrapper: 接受 params::NamedTuple (包含 :T 和 :mu)，返回 NamedTuple 包含 :Omega 和 :zero
function solver_wrapper(params, g0=nothing; solver_options...)
    T = params.T/197.33
    mu_scalar = params.mu/197.33
    # 将标量 mu 扩展为三味向量
    mu_vec = [mu_scalar, mu_scalar, mu_scalar]

    # 初始猜测
    X0 = g0 === nothing ? zeros(5) : g0
    # convert types to match T
    X0_typed = convert.(promote_type(eltype(X0), typeof(T)), X0)

    try
        res = nlsolve(x -> calculate_core_modern(x, mu_vec, T, p_grid, t_grid, xi_default), X0_typed; autodiff=:forward)
        converged_flag = converged(res)
        zero = res.zero
        omega = pressure_wrapper_modern(zero, mu_vec, T, p_grid, t_grid, xi_default)
        return (Omega = omega, zero = zero, converged = converged_flag, niter = res.iterations, residual = res.residual_norm)
    catch e
        @printf("solver exception at T=%.3f mu=%.3f : %s\n", T, mu_scalar, sprint(showerror,e))
        return (Omega = NaN, zero = nothing, converged = false, niter = 0, residual = Inf)
    end
end

# small grid for smoke test
T_values = collect(range(100.0, stop=150.0, length=2))
mu_values = collect(range(0.0, stop=50.0, length=2))
axes = [ParameterAxis("T", T_values), ParameterAxis("mu", mu_values)]

outfile = joinpath(@__DIR__, "..", "output", "phase_scan_tmu.csv")

# run scan
scan_phase_space(s -> begin
        # PhaseScanInterface passes a NamedTuple params; just forward to solver via wrapper
        # we expect scan to call solver when provided; here we pass the solver argument instead
        end,
    axes; solver=solver_wrapper, initial_guess=initial_guess, outfile=outfile, outputs=[:Omega])

@printf("phase scan finished, results in %s\n", outfile)
