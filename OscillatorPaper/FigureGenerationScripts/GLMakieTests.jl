begin 
    using Plots
    using GLMakie; GLMakie.activate!()
    using Catalyst
    using DifferentialEquations
    using Statistics
    using Peaks
    using Evolutionary, FFTW
    using Random
    using Distributions
    using DataFrames, DataFrameMacros
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
ogsol = solve(ogprob, saveat = 0.1)
Plots.plot(ogsol)

Base.@kwdef mutable struct Lorenz
    dt::Float64 = 0.01
    σ::Float64 = 10
    ρ::Float64 = 28
    β::Float64 = 8/3
    x::Float64 = 1
    y::Float64 = 1
    z::Float64 = 1
end

function step!(l::Lorenz)
    dx = l.σ * (l.y - l.x)
    dy = l.x * (l.ρ - l.z) - l.y
    dz = l.x * l.y - l.β * l.z
    l.x += l.dt * dx
    l.y += l.dt * dy
    l.z += l.dt * dz
    Point3f(l.x, l.y, l.z)
end

attractor = Lorenz()

function ODE_movie(sol)
    points = Observable(Point3f[]) # Signal that can be used to update plots efficiently
    colors = Observable(Int[])

    set_theme!(theme_black())

    fig, ax, l = lines(points, color = colors,
        colormap = :inferno, transparency = true, 
        axis = (; type = Axis3, protrusions = (0, 0, 0, 0), 
                viewmode = :fit, limits = (1, 4, 0, 1, 0, 2)))

    record(fig, "odetest.mp4", 1:120) do frame
        for i in eachindex(ogsol.t)
            # update arrays inplace
            push!(points[], Point3f(sol.u[i][1], sol.u[i][2], sol.u[i][4]))
            push!(colors[], frame)
        end
        ax.azimuth[] = 1.7pi + 0.3 * sin(2pi * frame / 120) # set the view angle of the axis
        notify(points); notify(colors) # tell points and colors that their value has been updated
        l.colorrange = (0, frame) # update plot attribute directly
    end
end
ODE_movie(ogsol)