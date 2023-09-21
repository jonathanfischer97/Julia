begin 
    using Plots; #theme(:juno)
    using Catalyst
    using OrdinaryDiffEq, ModelingToolkit
    # using DiffEqCallbacks
    using Statistics

    using Random
    using Distributions
    using DataFrames#, DataFrameMacros
    using CSV
    using StaticArrays
    using BenchmarkTools, ProgressMeter

    # using JET

    using LinearAlgebra

    using MultivariateStats
    using Clustering

    # using Setfield
    
    using ColorSchemes, Plots.PlotMeasures
    default(lw = 2, size = (1000, 600), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)



    using OscTools 


    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))

    using FFTW
    FFTW.set_num_threads(18)
end

function test4fixedGA()
    #* Set up the default GA problem
    ga_problem = GAProblem()

    #* Fixed some constraints 
    set_fixed_constraints!(ga_problem.constraints, [:DF, :K, :P, :A])

    #* Assign the fixed values 
    set_fixed_values!(ga_problem.constraints, 1000., 1.0, 1.0, 3.16)

    #* Set seed 
    Random.seed!(1234)

    #* Generate the initial population
    population = generate_population(ga_problem.constraints, 10000)

    #* Run the GA
    run_GA(ga_problem, population)
end

ga_result = test4fixedGA()





#* Clustering test 
dfarray = read_csvs_in_directory("/home/local/WIN/jfisch27/Desktop/Julia/OscillatorPaper/FigureGenerationScripts/ROCKFISH_DATA/4Fixed/PopSize_10000/4FixedICRawSets/DF=1000.0")


function get_optimal_k(df::DataFrame, max_k::Int; exclude_cols::Vector{Symbol} = Symbol[])
    distortions = elbow_method(df, max_k, exclude_cols = exclude_cols)
    return argmin(distortions)
end

 
clusterarray = [kmeans(df, 3, exclude_cols = [:fit, :per, :amp, :DF, :L, :K, :P, :A]) for df in dfarray]

df = dfarray[1]
optimal_k = get_optimal_k(df, 10, exclude_cols = [:fit, :per, :amp, :DF, :L, :K, :P, :A])





result = kmeans(df, 3, exclude_cols = [:fit, :per, :amp, :DF, :L, :K, :P, :A])

pca_mat = df_to_matrix(df, [:fit, :per, :amp, :DF, :L, :K, :P, :A])

pca_mat_transposed = transpose(pca_mat)

pca_model = fit(PCA, pca_mat_transposed, maxoutdim=2)
reduced_data = MultivariateStats.transform(pca_model, pca_mat_transposed)
reduced_centroids = MultivariateStats.transform(pca_model, result.centers)

#* Extract the coordinates of the reduced data
x_coords = reduced_data[1, :]
y_coords = reduced_data[2, :]

#* Extract the coordinates of the reduced centroids
centroid_x_coords = reduced_centroids[1, :]
centroid_y_coords = reduced_centroids[2, :]

#* Create the scatter plot for data points
scatter(x_coords, y_coords, zcolor=result.assignments, colorbar=true, label="", xlabel="PC1", ylabel="PC2", title="PCA Scatter Plot with K-means Clusters")

#* Overlay the centroids
scatter!(centroid_x_coords, centroid_y_coords, marker=:x, label="Centroids", color=:black)







#*FFTW testing 

ogprobjac = make_ODE_problem()
sol = solve_odeprob(ogprobjac, [6, 9, 10, 11, 12, 15, 16])
Amem_sol = map(sum, sol.u)

@btime rfft($Amem_sol)
rfft_plan = plan_rfft(Amem_sol)

@btime rfft_plan * $Amem_sol

@btime dct($Amem_sol)

dct_plan = plan_dct(Amem_sol)

@btime dct_plan * $Amem_sol

dct_plan = plan_dct!(Amem_sol)

@btime dct_plan * $Amem_sol


