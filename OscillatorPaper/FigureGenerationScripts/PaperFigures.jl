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
    using Unitful: µM, nm, s
    using BenchmarkTools, Profile, ProgressMeter
    # using MultivariateStats, UMAP, TSne, StatsPlots
    # using GlobalSensitivity, QuasiMonteCarlo
    using LinearAlgebra
    # using BifurcationKit, Setfield, ForwardDiff, Parameters; const BK = BifurcationKit
    using ColorSchemes, Plots.PlotMeasures
    using OrderedCollections
    using Combinatorics
    # using LazySets, Polyhedra
    default(lw = 2, size = (1000, 600))
end

#! PLOTTING BACKEND ##
plotlyjs()

"""Full oscillator model"""
fullrn = @reaction_network fullrn begin
    @parameters ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 ka7 kb7 kcat7 y
    @species L(t) Lp(t) K(t) P(t) A(t) LpA(t) LK(t) LpP(t) LpAK(t) LpAP(t) LpAKL(t) LpAPLp(t) AK(t) AP(t) AKL(t) APLp(t)
    # ALIASES: L = PIP, Lp = PIP2, K = Kinase, P = Phosphatase, A = AP2 
    # reactions between the same binding interfaces will have the same rate constant no matter the dimensionality or complex
    (ka1,kb1), L + K <--> LK # L binding to kinase
    kcat1, LK --> Lp + K # L phosphorylation by kinase into Lp
    (ka2,kb2), Lp + A <--> LpA # Lp binding to AP2 adaptor
    (ka3,kb3), LpA + K <--> LpAK # Membrane-bound adaptor binding to kinase
    (ka1*y,kb1), LpAK + L <--> LpAKL # 2D reaction: Membrane-bound kinase binds to L with greater affinity as determined by y (V/A)
    kcat1, LpAKL --> Lp + LpAK # L phosphorylation by kinase into Lp, same as 3D: first order reactions aren't dependent on dimensionality 
    (ka7,kb7), Lp + P <--> LpP # Lp binding to phosphatase
    kcat7, LpP --> L + P # L dephosphorylation by phosphatase
    (ka4,kb4), LpA + P <--> LpAP # Membrane-bound adaptor binding to phosphatase 
    (ka7*y,kb7), Lp + LpAP <--> LpAPLp # 2D reaction: Membrane-bound phosphatase binds to Lp with greater affinity as determined by y (V/A)
    kcat7, LpAPLp --> L + LpAP # L dephosphorylation by phosphatase, same as 3D: first order reactions aren't dependent on dimensionality

    #previously excluded reactions, all possible combinations possible in vitro
    (ka2,kb2), Lp + AK <--> LpAK
    (ka2*y,kb2), Lp + AKL <--> LpAKL
    (ka2,kb2), Lp + AP <--> LpAP
    (ka2*y,kb2), Lp + APLp <--> LpAPLp
    (ka3,kb3), A + K <--> AK
    (ka4,kb4), A + P <--> AP
    (ka3,kb3), A + LK <--> AKL
    (ka4,kb4), A + LpP <--> APLp
    (ka3*y,kb3), LpA + LK <--> LpAKL
    (ka4*y,kb4), LpA + LpP <--> LpAPLp
    (ka1,kb1), AK + L <--> AKL #binding of kinase to lipid
    kcat1, AKL --> Lp + AK #phosphorylation of lipid
    (ka7,kb7), AP + Lp <--> APLp #binding of phosphatase to lipid
    kcat7, APLp --> L + AP #dephosphorylation of lipid
end  

#! Solve model for arbitrary oscillatory parameters and initial conditions
begin
#parameter list
    const param_names = ["ka1", "kb1", "kcat1", "ka2", "kb2", "ka3", "kb3", "ka4", "kb4", "ka7", "kb7", "kcat7", "V/A"]

    rawpvals = [7.345728551517871, 758.1068961118469, 1.2291947166253903, 0.009153774616923016, 0.5009816004035071, 0.6179849986163015, 0.001, 0.32594456847574915, 0.001, 0.07113844033009084, 0.001, 2.521341233885747, 19975.20816768003]
    psym = [:ka1 => rawpvals[1], :kb1 => rawpvals[2], :kcat1 => rawpvals[3], :ka2 => rawpvals[4], :kb2 => rawpvals[5], :ka3 => rawpvals[6], :kb3 => rawpvals[7], :ka4 => rawpvals[8],
            :kb4 => rawpvals[9], :ka7 => rawpvals[10], :kb7 => rawpvals[11], :kcat7 => rawpvals[12], :y => rawpvals[13]]
    p = [x[2] for x in psym]
        
    #initial condition list
    usym = [:L => 0.0, :Lp => 3.0, :K => 0.5, :P => 0.3, :A => 2.0, :LpA => 0.0, :LK => 0.0, 
            :LpP => 0.0, :LpAK => 0.0, :LpAP => 0.0, :LpAKL => 0.0, :LpAPLp => 0.0, :AK => 0.0, :AP => 0.0, 
            :AKL => 0.0, :APLp => 0.0]
    u0 = [x[2] for x in usym]

    #timespan for integration
    tspan = (0., 100.)
    #solve the reduced ODEs
    prob = ODEProblem(fullrn, u0, tspan, p)
    sol = solve(prob, save_idxs=(1:6)) #solve adaptively
end


#* plot the results and save
legendlabels = ["PIP" "PIP2" "Kinase" "Phosphatase" "AP2" "PIP2-AP2"]
plot(sol, xlabel= "Time (s)", ylabel = "Concentration (µM)", label=legendlabels, grid= false, guidefontsize=16, tickfontsize=12, legendfontsize=12, thickness_scaling=2.5)
savefig("OscillatorPaper/FigureGenerationScripts/Figures/Figure1-timeseries.png")



plotattr()



# Sample data for demonstration purposes
oscillators = ["Membrane-localization oscillator", "Autocatalytic protein oscillator", "Repressilator"]
max_frequencies = [10.0, 8.0, 3.5]  # Example maximum frequencies for the oscillators
tunability_indices = [0.8, 0.3, 0.1]  # Example tunability indices for the oscillators

function plot_max_frequencies(oscillators, max_frequencies)
    p1 = bar(oscillators, max_frequencies, ylabel="Max Frequency (Hz)",fillcolor=[:red,:green,:blue], rotation=15,  legend=false)
    return p1
end

function plot_tunability_indices(oscillators, tunability_indices)
    p2 = bar(oscillators, tunability_indices, ylabel="Tunability Index",fillcolor=[:red,:green,:blue], rotation=15,legend=false)
    return p2
end

p1 = plot_max_frequencies(oscillators, max_frequencies)
p2 = plot_tunability_indices(oscillators, tunability_indices)

plot(p1, p2, layout=(1, 2), size=(1000, 600), title=["Max Frequencies" "Tunability Indices"],bottom_margin = 16px, left_margin = 18px)
savefig("OscillatorPaper/FigureGenerationScripts/Figures/Figure9-compare-oscillators.png")