# Runner for PhaseScanInterface
# 使用方式： julia --project=. scripts/phase_scan_runner.jl

using Dates
include(joinpath(@__DIR__, "..", "src", "core", "phase_scan_interface.jl"))
using .PhaseScanInterface

# 占位物理函数（或可替换为项目中的函数）
function placeholder_omega(params::NamedTuple)
    T = haskey(params, :T) ? params.T : 300.0
    mu = haskey(params, :mu) ? params.mu : 0.0
    return sin(T/100.0) + cos(mu/200.0)
end

# 生成两个轴并运行
T_axis = ParameterAxis("T", collect(range(0.0, stop=400.0, length=40)))
mu_axis = ParameterAxis("mu", collect(range(0.0, stop=500.0, length=40)))
axes = [T_axis, mu_axis]

outfile = joinpath(@__DIR__, "..", "output", "phase_scan_from_core.csv")

println("开始运行 phase scan runner...")
scan_phase_space(placeholder_omega, axes; outfile=outfile)
println("runner 完成。 输出: $outfile")
