using Test
using Dates

# Tests for PhaseScanInterface
include(joinpath(@__DIR__, "..", "src", "core", "phase_scan_interface.jl"))
include(joinpath(@__DIR__, "..", "src", "core", "solver_interface.jl"))
using .PhaseScanInterface
using .SolverInterface

@testset "PhaseScanInterface basic" begin
    # ParameterAxis
    a = ParameterAxis("T", [1.0, 2.0])
    @test a.name == "T"
    @test length(a.values) == 2

    # generate_grid
    b = ParameterAxis("mu", [10.0, 20.0])
    axes = [a, b]
    it = generate_grid(axes)
    vals = collect(it)
    @test length(vals) == 4
    @test isa(vals[1], NamedTuple)
    @test haskey(vals[1], :T) && haskey(vals[1], :mu)
    @test vals[1].T == 1.0 && vals[1].mu == 10.0

    # wrap_physics_fn: different signatures
    f1 = params -> params.T + params.mu
    w1 = wrap_physics_fn(f1)
    @test w1((T=1.0, mu=2.0)) == 3.0

    f2 = (T, mu) -> T + mu
    w2 = wrap_physics_fn(f2)
    @test w2((T=4.0, mu=1.0)) == 5.0

end

@testset "PhaseScanInterface integration (scan_phase_space)" begin
    # ensure output dir
    outdir = joinpath(@__DIR__, "..", "output")
    try
        mkpath(outdir)
    catch
    end
    outfile = joinpath(outdir, "test_phase_scan_interface.csv")
    metafile = outfile * ".meta.json"
    # remove existing
    try
        isfile(outfile) && rm(outfile)
        isfile(metafile) && rm(metafile)
    catch
    end

    # dummy solver implemented via EquationFactory -> SolverInterface (complete chain)
    include(joinpath(@__DIR__, "..", "src", "core", "equation_factory.jl"))
    using .EquationFactory

    # define simple linear sub-equations compatible with factory
    f(x, y, T, mu) = x - T*0.1
    # param_names used by factory are (:T, :mu), so function must accept params in that order
    g(y, x, T, mu) = y - mu*0.01

    # declare local param_names per sub-equation: f needs (:T,), g needs (:mu,)
    F_param!, names_f, param_names_f, pack_f, unpack_f = EquationFactory.create_equation_system_inplace(
        (f, (:x, :y), (:T,)),
        (g, (:y, :x), (:mu,))
    )

    @test param_names_f == (:T, :mu)

    function dummy_solver(params, init; kwargs...)
        g0 = init isa AbstractVector ? init : [0.0, 0.0]
        # call solver with factory-produced in-place function and pack
        return solve_equilibrium_equations(F_param!, g0, params; pack=pack_f)
    end

    # postprocess: compute Omega from solution vector or named fields
    function my_postprocess(solution, params)
        if isa(solution, AbstractVector)
            x = solution[1]; y = solution[2]
        elseif isa(solution, NamedTuple) && haskey(solution, :x)
            x = solution.x; y = solution.y
        elseif isa(solution, Dict) && haskey(solution, :x)
            x = solution[:x]; y = solution[:y]
        else
            return Dict(:Omega => NaN)
        end
        return Dict(:Omega => x^2 + y^2)
    end

    axes = [ParameterAxis("T", [100.0, 200.0]), ParameterAxis("mu", [10.0])]

    scan_phase_space((p)->0.0, axes; solver = dummy_solver, outfile = outfile, outputs = [:x, :y, :Omega], header_probe = :first_point, postprocess = my_postprocess)

    @test isfile(outfile)
    @test isfile(metafile)

    # check header contains parameter names and output fields
    open(outfile, "r") do io
        header = readline(io)
        @test contains(header, "T")
        @test contains(header, "mu")
        @test contains(header, "x")
        @test contains(header, "y")
    end

end
