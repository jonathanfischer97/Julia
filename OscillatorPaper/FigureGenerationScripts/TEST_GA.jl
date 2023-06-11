begin 
    using Plots 
    using Catalyst
    using DifferentialEquations
    using Statistics
    using Peaks
    using Evolutionary, FFTW
    using Random
    using Distributions
    using DataFrames
    # using Unitful
    # using Unitful: µM, nm, s
    using StaticArrays
    using ProgressMeter
    # using BenchmarkTools, Profile, ProgressMeter
    # using Cthulhu
    # using JET
    # using MultivariateStats, UMAP, TSne, StatsPlots
    # using GlobalSensitivity, QuasiMonteCarlo
    using LinearAlgebra
    # using BifurcationKit, Setfield, ForwardDiff, Parameters; const BK = BifurcationKit
    using ColorSchemes, Plots.PlotMeasures
    # using OrderedCollections
    using Combinatorics
    # using LazySets, Polyhedra
    default(lw = 2, size = (1000, 600))
    # plotlyjs()
    # gr()
    # using Base.Threads


    include("../../UTILITIES/EvolutionaryOverloads.jl")

    # import the Catalyst model "fullrn"
    include("../../UTILITIES/ReactionNetwork.jl")

    # import the cost function and other evaluation functions
    include("../../UTILITIES/EvaluationFunctions.jl")

    # import the genetic algorithm and associated functions
    include("../../UTILITIES/GA_functions.jl")

    include("../../UTILITIES/ODEProbMaker.jl")

    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))
end


param_constraints = define_parameter_constraints()

# ic_constraints = define_initialcondition_constraints()


testprob = make_ODEProb()

paramgaproblem = GAProblem(param_constraints, testprob)

# testsol = solve(testprob, saveat = 0.1)
# plot(testsol)

testrun = run_GA(paramgaproblem; population_size = 5000)

# plotsol(row) = plotsol(row, testrun, testprob; vars = [1])

# plotsol(1)


# testsol = solve(testprob, saveat = 0.1)[1,:]
# plot(testsol)

# argmaxima(testsol, 10)




extract_solution(row) = extract_solution(row, testrun, testprob; vars = [1])

dfsol = extract_solution(6);
plot(dfsol)

findmaxima(dfsol[1,:], 1)

dfsolfreqs = getFrequencies(dfsol[1,:])

plot(dfsolfreqs, title = "FFT of Solution", xlabel = "Frequency", ylabel = "Magnitude")

argmaxima(dfsolfreqs, 10)

CostFunction(dfsol)