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


"""
    read_csvs_in_directory(directory_path::String)

Read all CSV files in a given directory into an array of DataFrames.

# Arguments
- `directory_path::String`: The path to the directory containing the CSV files.

# Returns
- `Array{DataFrame, 1}`: An array of DataFrames, each corresponding to a CSV file.
"""
function read_csvs_in_directory(directory_path::String)
    # Initialize an empty array to store DataFrames
    dfs = DataFrame[]
    
    # Loop over each file in the directory
    for file_name in readdir(directory_path)
        # Check if the file is a CSV file
        if occursin(".csv", file_name)
            # Full path to the CSV file
            full_file_path = joinpath(directory_path, file_name)
            
            # Read the CSV file into a DataFrame
            df = DataFrame(CSV.File(full_file_path))
            
            # Append the DataFrame to the array
            push!(dfs, df)
        end
    end
    
    return dfs
end

# Example usage
dfs_array = read_csvs_in_directory("OscillatorPaper/FigureGenerationScripts/ROCKFISH_DATA/3Fixed/PopSize_10000/K_P_A/SummaryResults")


using GLMakie; GLMakie.activate!()
using GeometryBasics: Rect3f

# Assuming dfs_array is your array of DataFrames
# For demonstration, let's assume you've read one DataFrame from dfs_array
df = dfs_array[3]

# Log-transform the fixed values for plotting
df_log = transform(df, :K => ByRow(log10), :P => ByRow(log10), :A => ByRow(log10))

# Marker size normalization and scaling
marker_size = df.num_oscillatory_points ./ maximum(df.num_oscillatory_points) * 0.3

# Log-transform for better color mapping
color_values = log10.(df.average_period)

# Create the meshscatter plot
fig, ax, pltobj = meshscatter(
    vec([(df_log.K_log10[i], df_log.P_log10[i], df_log.A_log10[i]) for i in 1:length(df_log.K)]);
    color = color_values,
    # marker = Rect3f(Vec3f(-0.5), Vec3f(1)),
    markersize = marker_size,
    colormap = :viridis,
    figure = (;
        resolution = (1200, 800)),
    axis = (;
        type = Axis3,
        xticklabels = df.K,
        xlabel = "Log(K)",
        ylabel = "Log(P)",
        zlabel = "Log(A)",
        aspect = :data)
)

# Display the plot
fig



function create_meshscatter_plot(df::DataFrame)
    # Log-transform the fixed values for plotting
    df_log = transform(df, :K => ByRow(log10), :P => ByRow(log10), :A => ByRow(log10))

    # Marker size normalization and scaling
    marker_size = df.num_oscillatory_points ./ maximum(df.num_oscillatory_points) * 0.3

    # Log-transform for better color mapping
    color_values = log10.(df.average_period)

    # Create the meshscatter plot
    fig, ax, pltobj = meshscatter(
        vec([(df_log.K_log10[i], df_log.P_log10[i], df_log.A_log10[i]) for i in 1:length(df_log.K)]);
        color = color_values,
        markersize = marker_size,
        colormap = (:Egypt, 1.0),
        figure = (;
            resolution = (1200, 800)),
        axis = (;
            type = Axis3,
            xtickformat = values -> [string(round(10^x; digits=2)) for x in values], #* converts the log10 values back to the original values
            ytickformat = values -> [string(round(10^x; digits=2)) for x in values],
            ztickformat = values -> [string(round(10^x; digits=2)) for x in values],
            xlabel = "K",
            ylabel = "P",
            zlabel = "A",
            aspect = :data)
    )
    # Display the plot
    return fig
end

create_meshscatter_plot(df)




df100 = CSV.read("OscillatorPaper/FigureGenerationScripts/ROCKFISH_DATA/4Fixed/PopSize_1000/SummaryResults/Summary_DF=100.0.csv", DataFrame)
df1000 = CSV.read("OscillatorPaper/FigureGenerationScripts/ROCKFISH_DATA/4Fixed/PopSize_1000/SummaryResults/Summary_DF=1000.0.csv", DataFrame)
df10000 = CSV.read("OscillatorPaper/FigureGenerationScripts/ROCKFISH_DATA/4Fixed/PopSize_1000/SummaryResults/Summary_DF=10000.0.csv", DataFrame)

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
