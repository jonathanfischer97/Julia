begin 
    using Plots; #theme(:juno)
    using Catalyst
    using DifferentialEquations
    using Statistics
    using Peaks
    using Evolutionary, FFTW
    using Random
    using Distributions
    using DataFrames
    # using CSV
    # using Unitful
    # using Unitful: µM, M, nm, µm, s, μs, Na, L, 𝐍
    # using StaticArrays
    using BenchmarkTools, Profile, ProgressMeter
    # using Cthulhu
    # using JET
    # using MultivariateStats, UMAP, TSne, StatsPlots
    # using GlobalSensitivity, QuasiMonteCarlo
    using LinearAlgebra
    # using BifurcationKit, Setfield, ForwardDiff, Parameters; const BK = BifurcationKit
    # using ColorSchemes, Plots.PlotMeasures
    # using OrderedCollections
    # using Combinatorics
    # using LazySets, Polyhedra
    # default(lw = 2, size = (1000, 600), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)
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

    include("../../UTILITIES/ODEProbMaker.jl")

    include("../../UTILITIES/UnitTools.jl")


    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))
end

fullrn = make_fullrn()
ogprob = ODEProblem(fullrn, [], (0.,100.0), [])

param_constraints = define_parameter_constraints(karange = (1e-3, 1e2), kbrange = (1e-2, 1e3), kcatrange = (1e-2, 1e3), dfrange = (1e3, 1e5))

#! TESTING GA FUNCTIONALITY
test_gaproblem = GAProblem(param_constraints, ogprob)

test_results = run_GA(test_gaproblem; population_size = 5000, iterations = 1, parallelization=:serial)

testind = test_results.trace[2].metadata["staterecord"][1].ind
testsol = solve(remake(ogprob, p=testind), saveat = 0.01, save_idxs = 1)
CostFunction(testsol)

#TODO fix the evaluation of the first generation in initial_state
#* have ruled out the in-place objective function, as well as my replacing the map with a for loop
#* doesn't seem to be an issue with the trace recording either 
#* something wrong at the time of evaluation within initial_state


eval_param_fitness(test_results.ind[2], ogprob)



plotsol(row) = plotsol(row, test_results, ogprob)

plotsol(1)


test_fitness(row) = eval_param_fitness(test_results.ind[row], ogprob)
test_fitness(500)

reogprob = remake(ogprob, p=test_results.ind[1])
testsol = solve(reogprob, saveat = 0.01, save_idxs = 1)
CostFunction(testsol)


plot(testsol)


fitind = [0.005904610505030443, 131.13364934592303, 314.4503575261949, 
        3.3658232279193547, 3.0772188176974846, 1.6497701059703715, 
        0.22809186527536418, 16.53422149313229, 0.8938599354855014, 
        3.34373724825129, 0.6270815692033677, 101.721977264349, 44216.07505733261]

fitprob = remake(ogprob, p=fitind)
fitsol = solve(fitprob, saveat = 0.01, save_idxs = 1)
CostFunction(fitsol)

