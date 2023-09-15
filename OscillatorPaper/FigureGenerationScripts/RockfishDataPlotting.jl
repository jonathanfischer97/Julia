begin 
    using Plots; #theme(:juno)
    using StatsPlots
    using Catalyst
    using OrdinaryDiffEq, ModelingToolkit
    using Statistics
    using Evolutionary, FFTW
    using Random
    using Distributions
    using DataFrames, DataFrameMacros
    using CSV
    using StaticArrays
    using BenchmarkTools, ProgressMeter

    # using JET

    using LinearAlgebra

    # using Setfield
    
    using ColorSchemes, Plots.PlotMeasures
    default(lw = 2, size = (1000, 600), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)

    #* import the overloads for Evolutionary.jl
    include("../../UTILITIES/EvolutionaryOverloads.jl")

    #* import the Catalyst model "fullrn"
    include("../../UTILITIES/ReactionNetwork.jl")

    #* import the cost function and other evaluation functions
    include("../../UTILITIES/EvaluationFunctions.jl")
    # using .EvaluationFunctions

    #* import the genetic algorithm and associated functions
    include("../../UTILITIES/GA_functions.jl")

    #* import the plotting functions
    include("../../UTILITIES/TestBenchPlotUtils.jl")


    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))
end


df100 = CSV.read("OscillatorPaper/FigureGenerationScripts/ROCKFISH_DATA/4Fixed/PopSize_1000/SummaryResults/Summary_DF=100.0.csv", DataFrame)
df1000 = CSV.read("OscillatorPaper/FigureGenerationScripts/ROCKFISH_DATA/4Fixed/PopSize_1000/SummaryResults/Summary_DF=1000.0.csv", DataFrame)
df10000 = CSV.read("OscillatorPaper/FigureGenerationScripts/ROCKFISH_DATA/4Fixed/PopSize_1000/SummaryResults/Summary_DF=10000.0.csv", DataFrame)


using GLMakie; GLMakie.activate!()

function plot_4fixed_makie(df, Lval = 100.)
    df = df[df.L .== Lval, :]

    xlog = log10.(df.K) 
    ylog = log10.(df.P) 
    zlog = log10.(df.A)


    #* make the plot
    fig = Figure(resolution = (1000, 600))

    #* make the scatter plot
    # Identify non-NaN indices and values
    nonan_indices = findall(!isnan, df[:, :average_period])
    # nonan_amplitudes = df[:, :average_amplitude][nonan_indices]
    nonan_numpoints = df[:, :num_oscillatory_points][nonan_indices]

    nonan_periods = df[:, :average_period][nonan_indices]

    # Normalize sizes for non-NaN values
    sizes = fill(0.1, size(df, 1))
    # sizes[nonan_indices] = ((nonan_amplitudes .- minimum(nonan_amplitudes)) ./ (maximum(nonan_amplitudes) - minimum(nonan_amplitudes))) ./ 2 
    sizes[nonan_indices] = ((nonan_numpoints .- minimum(nonan_numpoints)) ./ (maximum(nonan_numpoints) - minimum(nonan_numpoints))) ./ 3 


    # Normalize periods for non-NaN values
    norm_periods = fill(NaN, size(df, 1))

    norm_periods[nonan_indices] = (nonan_periods .- minimum(nonan_periods)) ./ (maximum(nonan_periods) - minimum(nonan_periods)) 

    # Create the figure and axis
    
    ax = Axis3(fig[1:3,2]; aspect=:data, perspectiveness=0.5, title="K vs. P vs. A for L=$(round(Lval))uM, DF=10000.", xlabel = "K", ylabel = "P", zlabel = "A")

    # Scatter plot for non-NaN values
    hm = meshscatter!(ax, xlog, ylog, zlog; markersize=sizes, ssao=true, color=df.average_period, colormap=:thermal, transparency=true, nan_color=:gray,
                        diffuse = Vec3f(0.5, 0.5, 0.5), specular = Vec3f(0.3, 0.3, 0.3), shininess = 100f0, ambient = Vec3f(0.1), shading=true, alpha=0.7)

    # Colorbar and labels
    Colorbar(fig[2, 4], hm, label="Period (s)", height=Relative(2.0))
    colgap!(fig.layout, 6)

    fig
end

plot_4fixed_makie(df10000, 46.4158883361278)


unique(df100.L)

nonan_indices = findall(!isnan, df10000[:, :average_period])


df10000[nonan_indices,:]


testdf = CSV.read("OscillatorPaper/FigureGenerationScripts/ROCKFISH_DATA/4Fixed/PopSize_1000/4FixedICRawSets/DF=10000.0/L=0.1_K=0.01_P=0.03_A=2.15.csv", DataFrame)

plot_everything(testdf, make_ODE_problem(); jump=1)
