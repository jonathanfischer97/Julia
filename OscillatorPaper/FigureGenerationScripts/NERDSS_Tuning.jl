begin 
    using Plots 
    using Catalyst
    using DifferentialEquations
    using Statistics
    # using Peaks
    # using Evolutionary, FFTW
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
    # using Combinatorics
    # using LazySets, Polyhedra
    default(lw = 2, size = (1000, 600), dpi = 200)
    # plotlyjs()
    import CairoMakie as cm 
    gr()


    # import the Catalyst model "fullrn"
    include("../../UTILITIES/ReactionNetwork.jl")

    # import the ODE problem generator. Sets `fullprob` as a constant
    include("../../UTILITIES/ODEProbMaker.jl")

    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))
end


const fullrn = make_fullrn()

format_equations(fullrn)

psym = [:ka1 => 5.453534109021441e-05 #ka1, 1
    :kb1 => 0.0643048008980449 #kb1, 2
    :kcat1 => 286.6382253193995 #kcat1, 3
    :ka2 => 0.0033210781343476934 #ka2, 4
    :kb2 => 0.39569337786534897 #kb2, 5
    :ka3 => 8.23096854112845e-05 #ka3, 6
    :kb3 => 0.5393197910059361 #kb3, 7
    :ka4 => 0.0001089706534070095 #ka4, 8
    :kb4 => 0.2897657637531564 #kb4, 9
    :ka7 => 0.0003802716418207893 #ka7, 10
    :kb7 => 0.0028126177315505618 #kb7, 11
    :kcat7 => 1.2733781341040291 #kcat7, 12
    :DF => 2650.5049775102034] #DF, 13
p = [x[2] for x in psym]
usym = [:L => 634.6315289074139, :K => 47.42150049582334, :P => 239.66312964177104,  :A => 838.7724702072363, :Lp => 790.5014385747756, :LpA => 0.0, :LK => 0.0, 
        :LpP => 0.0, :LpAK => 0.0, :LpAP => 0.0, :LpAKL => 0.0, :LpAPLp => 0.0, :AK => 0.0, :AP => 0.0, :AKL => 0.0, :APLp => 0.0]

nerdssprob = make_ODEProb(; psym = psym, usym = usym, tspan = (0., 200.))

nerdsol = solve(nerdssprob, Rosenbrock23(), saveat = 0.01)

# plot(nerdsol, idxs = [1,5,2,3,4], color = [:blue :gold :red :purple :green])


pl = plot(nerdsol, idxs = (1,3), xlabel = "PIP2 (uM)", ylabel = "Synaptojanin (µM)", title = "Tuning the period through V/A", label = "Nominal")

for i in 1:5
    p[end] = p[end] * 0.8;
    @info "Running with V/A = $(p[end])"
    new_nerdsol = solve(remake(nerdssprob; p = p), Rosenbrock23(), saveat = 0.01);
    plot!(pl, new_nerdsol, idxs = (1,3), xlims = (300, 1300),label = "V/A scaled by 0.8^$i", ls = :dash)
end
display(pl)
savefig(pl, "OscillatorPaper/FigureGenerationScripts/ProgressReportFigures/NERDSS_Tuning.png")
