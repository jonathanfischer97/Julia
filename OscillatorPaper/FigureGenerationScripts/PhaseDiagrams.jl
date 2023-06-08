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
    using BenchmarkTools, Profile, ProgressMeter
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
    default(lw = 2, size = (1000, 600), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 8px)
    plotlyjs()
    using CairoMakie
    # gr()
    # using Base.Threads


    include("../../UTILITIES/EvolutionaryOverloads.jl")

    # import the Catalyst model "fullrn"
    include("../../UTILITIES/ReactionNetwork.jl")

    # import the cost function and other evaluation functions
    include("../../UTILITIES/EvaluationFunctions.jl")

    # import the genetic algorithm and associated functions
    include("../../UTILITIES/GA_functions.jl")

    # import the ODE problem generator. Sets `fullprob` as a constant
    include("../../UTILITIES/ODEProbMaker.jl")

    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))
end

#* Solve model for arbitrary oscillatory parameters and initial conditions
sol = solve(fullprob, saveat=0.1)

#* Plot the phase diagrams with time as the third dimension
Plots.plot(sol, idxs=(0, 1, 3), colorbar = false, xlabel = "Time", ylabel = "L", zlabel = "K", title = "Phase Diagrams")#, camera = (-100, 30))


#! Makie plotting 
begin
    t = sol.t
    x = sol[1]
    y = sol[2,:]
    z = sol[3,:]
end


f = Figure()
ax = Axis(f[1,1])
lines!(ax, t, x, color = :blue, linewidth = 3, linestyle = :dash)
CairoMakie.scatter!(ax, t, x, color = range(0,1,length=length(t)), markersize = range(5, 15, length=length(t)), colormap = :viridis)
f