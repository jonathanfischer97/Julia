begin 
    using Plots; #theme(:juno)
    # using Compose
    using Catalyst
    using DifferentialEquations, ModelingToolkit
    using Statistics
    # using Peaks
    using FindPeaks1D
    using Evolutionary, FFTW
    using Random
    using Distributions
    using DataFrames#, DataFrameMacros
    using CSV
    # using Unitful
    # using Unitful: µM, M, nm, µm, s, μs, Na, L, 𝐍
    using StaticArrays
    using BenchmarkTools, Profile, ProgressMeter
    # using Cthulhu
    # using JET
    # using MultivariateStats, UMAP, TSne, StatsPlots
    # using GlobalSensitivity, QuasiMonteCarlo
    using LinearAlgebra
    # using BifurcationKit, Setfield, ForwardDiff, Parameters; const BK = BifurcationKit
    using ColorSchemes, Plots.PlotMeasures
    # using OrderedCollections
    # using Combinatorics
    # using LazySets, Polyhedra
    default(lw = 2, size = (1000, 600), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)
    # plotlyjs()
    # import CairoMakie as cm 
    # gr()
    # push!(LOAD_PATH, "../../UTILITIES")

    include("../../UTILITIES/EvolutionaryOverloads.jl")

    # import the Catalyst model "fullrn"
    include("../../UTILITIES/ReactionNetwork.jl")

    # import the cost function and other evaluation functions
    include("../../UTILITIES/EvaluationFunctions.jl")
    # using .EvaluationFunctions

    # import the genetic algorithm and associated functions
    include("../../UTILITIES/GA_functions.jl")

    include("../../UTILITIES/TestBenchPlotUtils.jl")

    # include("../../UTILITIES/UnitTools.jl")


    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))
end








function testbench(constraints::ConstraintType, prob::ODEProblem)
    test_gaproblem = GAProblem(constraints, prob)
    Random.seed!(1234)
    test_results = run_GA(test_gaproblem; population_size = 1000, iterations = 5, show_trace=true)

    if length(test_results.fitvals) == 0
        @info "No oscillatory points found."
    else 
        test_results_df = make_ga_dataframe(test_results, prob)
        return test_results_df
    end
end

function run_test()
    tspan = (0., 2000.0)
    fullrn = make_fullrn()
    ogprob = ODEProblem(fullrn, [], tspan, [])
    de = modelingtoolkitize(ogprob)

    ogprobjac = ODEProblem(de, [], tspan, jac=true)

    param_constraints = define_parameter_constraints(; karange = (1e-3, 1e2), kbrange = (1e-3, 1e3), kcatrange = (1e-3, 1e3), dfrange = (1e2, 2e4))
    ic_constraints = define_initialcondition_constraints(; Lrange = (1e-1, 1e2), Krange = (1e-2, 1e2), Prange = (1e-2, 1e2), Arange = (1e-1, 1e2))

    allconstraints = AllConstraints(param_constraints, ic_constraints)

    ogprobjac = remake(ogprobjac, u0 = [[100., 0.2, 0.2, 4.64]; zeros(12)])

    test_results_df = testbench(param_constraints, ogprobjac)
end


run_test()

# test_results_df = testbench(allconstraints, ogprobjac)

