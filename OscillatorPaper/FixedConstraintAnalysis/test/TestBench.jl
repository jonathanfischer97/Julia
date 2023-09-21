begin 
    using Plots; #theme(:juno)
    using Catalyst
    using OrdinaryDiffEq, ModelingToolkit
    using DiffEqCallbacks
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



#* Clustering test 
dfarray = read_csvs_in_directory("/home/local/WIN/jfisch27/Desktop/Julia/OscillatorPaper/FigureGenerationScripts/ROCKFISH_DATA/4Fixed/PopSize_10000/4FixedICRawSets/DF=1000.0")

function df_to_matrix(df::DataFrame, exclude_cols::Vector{Symbol})
    return Matrix(df[:, Not(exclude_cols)])
end

function kmeans(df::DataFrame, k::Int; exclude_cols::Vector{Symbol} = Symbol[])
    data_matrix = df_to_matrix(df, exclude_cols)
    return kmeans(data_matrix, k)
end

function elbow_method(df::DataFrame, max_k::Int; exclude_cols::Vector{Symbol} = Symbol[])
    data_matrix = df_to_matrix(df, exclude_cols)
    return elbow_method(data_matrix, max_k)
end

function elbow_method(X::Matrix{Float64}, max_k::Int)
    distortions = Vector{Float64}(undef, max_k)
    for k in 1:max_k
        result = kmeans(X, k)
        distortions[k] = result.totalcost
    end
    return distortions
end

function get_optimal_k(df::DataFrame, max_k::Int; exclude_cols::Vector{Symbol} = Symbol[])
    distortions = elbow_method(df, max_k, exclude_cols = exclude_cols)
    return argmin(distortions)
end

 
clusterarray = [kmeans(df, 3, exclude_cols = [:fit, :per, :amp, :DF, :L, :K, :P, :A]) for df in dfarray]

df = dfarray[1]
optimal_k = get_optimal_k(df, 10, exclude_cols = [:fit, :per, :amp, :DF, :L, :K, :P, :A])
cluster = clusterarray[1]

a = assignments(cluster)
c = counts(cluster)
m = cluster.centers
n = nclusters(cluster)

using MultivariateStats



pca_mat = df_to_matrix(df, [:fit, :per, :amp, :DF, :L, :K, :P, :A])

pca_model = fit(PCA, pca_mat, maxoutdim=2)
reduced_data = MultivariateStats.transform(pca_model, pca_mat)

x_coords = reduced_data[:, 1]
y_coords = reduced_data[:, 2]

scatter(x_coords, y_coords, zcolor=df.per, colorbar=true)


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


