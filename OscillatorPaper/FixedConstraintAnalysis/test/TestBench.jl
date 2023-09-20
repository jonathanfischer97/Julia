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

import Clustering: kmeans
function kmeans(df::DataFrame, k::Int; exclude_cols::Vector{Symbol} = Symbol[])
    included_cols = setdiff(propertynames(df), exclude_cols)
    data_matrix = Matrix(df[:, included_cols])
    return kmeans(data_matrix, k)
end

function elbow_method(df::DataFrame, max_k::Int; exclude_cols::Vector{Symbol} = Symbol[])
    included_cols = setdiff(propertynames(df), exclude_cols)
    data_matrix = Matrix(df[:, included_cols])
    return elbow_method(data_matrix, max_k)
end

function elbow_method(X::Matrix{Float64}, max_k::Int)
    distortions = Vector{Float64}(undef, max_k)
    for k in 1:max_k
        result = kmeans(X, k)
        distortions[k] = totalcost(result)
    end
    return distortions
end

 
clusterarray = [kmeans(df, 4, exclude_cols = [:fit, :per, :amp, :DF, :L, :K, :P, :A]) for df in dfarray]

df = dfarray[1]
cluster = clusterarray[1]

a = assignments(cluster)
c = counts(cluster)
m = cluster.centers
n = nclusters(cluster)

plot(df.ka1, df.kb1, marker_z = a, color=:lightrainbow)