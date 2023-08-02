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

    include("../../UTILITIES/ODEProbMaker.jl")

    include("../../UTILITIES/UnitTools.jl")


    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))
end


fullrn = make_fullrn()
ogprob = ODEProblem(fullrn, [], (0.,1000.0), [])

# ogsol = solve(ogprob, Rosenbrock23(), save_idxs = 1)
# plot(ogsol)


#* Optimization of initial conditions to produce data for CSV
param_constraints = define_parameter_constraints(ogprob)

param_gaproblem = GAProblem(param_constraints, ogprob)

param_record = run_GA(param_gaproblem; population_size = 50000)


param_names = [Symbol(x.name) for x in param_constraints.ranges]

csv_record = DataFrame()
for (i,name) in enumerate(param_names)
    csv_record[!, name] = [x[i] for x in param_record.ind]
end
csv_record.period = param_record.per
csv_record.amplitude = param_record.amp
csv_record

CSV.write("OscillatorPaper/FigureGenerationScripts/optimized_parameters.csv", csv_record)
