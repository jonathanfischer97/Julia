begin 
    # using Plots; #theme(:juno)
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
    
    # using ColorSchemes, Plots.PlotMeasures
    # default(lw = 2, size = (1000, 600), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)



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

population_to_matrix(ga_result)



#* Clustering test 
dfarray = read_csvs_in_directory("/home/local/WIN/jfisch27/Desktop/Julia/OscillatorPaper/FigureGenerationScripts/ROCKFISH_DATA/4Fixed/PopSize_10000/4FixedICRawSets/DF=1000.0")


function get_optimal_k(df::DataFrame, max_k::Int; exclude_cols::Vector{Symbol} = Symbol[])
    distortions = elbow_method(df, max_k, exclude_cols = exclude_cols)
    return argmin(distortions)
end

 
clusterarray = [kmeans(df, 3, exclude_cols = [:fit, :per, :amp, :DF, :L, :K, :P, :A]) for df in dfarray]

df = dfarray[1]
optimal_k = get_optimal_k(df, 10, exclude_cols = [:fit, :per, :amp, :DF, :L, :K, :P, :A])


df_mat = df_to_matrix(df, [:fit, :per, :amp, :DF, :L, :K, :P, :A])
result = kmeans(df_mat, 3)


pca_mat = df_to_matrix(df, [:fit, :per, :amp, :DF, :L, :K, :P, :A])

# pca_mat_transposed = transpose(pca_mat)

pca_model = fit(PCA, pca_mat, maxoutdim=3)
reduced_data = MultivariateStats.transform(pca_model, pca_mat)
reduced_centroids = MultivariateStats.transform(pca_model, result.centers)

#* Extract the coordinates of the reduced data
x_coords = reduced_data[1, :]
y_coords = reduced_data[2, :]
z_coords = reduced_data[3, :]

#* Extract the coordinates of the reduced centroids
centroid_x_coords = reduced_centroids[1, :]
centroid_y_coords = reduced_centroids[2, :]
centroid_z_coords = reduced_centroids[3, :]

#* Create the scatter plot for data points
scatter(x_coords, y_coords, z_coords, zcolor=result.assignments, colorbar=true, label="", xlabel="PC1", ylabel="PC2", title="PCA Scatter Plot with K-means Clusters")

#* Overlay the centroids
scatter!(centroid_x_coords, centroid_y_coords, centroid_z_coords, marker=:x, label="Centroids", color=:black)


# Function to identify fixed columns in a DataFrame
function identify_fixed_columns(df)
    fixed_cols = Symbol[]
    for col in propertynames(df)
        if length(unique(df[!, col])) == 1
            push!(fixed_cols, col)
        end
    end
    return fixed_cols
end


using GLMakie; GLMakie.activate!()

function plot_3D_PCA_clusters(df)
    # Step 1: Identify fixed columns
    fixed_cols = identify_fixed_columns(df)
    
    # Step 2: Perform K-means clustering
    result = kmeans(df, 3, exclude_cols = fixed_cols)
    
    # Step 3: Perform PCA
    pca_mat = df_to_matrix(df, fixed_cols)
    pca_model = fit(PCA, pca_mat, maxoutdim=3)
    reduced_data = MultivariateStats.transform(pca_model, pca_mat)
    reduced_centroids = MultivariateStats.transform(pca_model, result.centers)
    
    # Step 4: Extract 3D coordinates
    x_coords = reduced_data[1, :]
    y_coords = reduced_data[2, :]
    z_coords = reduced_data[3, :]

    centroid_x = reduced_centroids[1, :]
    centroid_y = reduced_centroids[2, :]
    centroid_z = reduced_centroids[3, :]
    
    # Step 5: Plotting
    fig = Figure(;resolution = (1000, 600))
    ax = Axis3(fig[1,1]; aspect=:data, perspectiveness = 0.5)
    GLMakie.scatter!(ax, x_coords, y_coords, z_coords, label="", xlabel="PC1", ylabel="PC2", title="PCA Scatter Plot with K-means Clusters")
    # GLMakie.scatter!(ax, centroid_x, centroid_y, centroid_z, marker=:x, label="Centroids", color=:black)

    
    return fig
end

plot_3D_PCA_clusters(dfarray[1])




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


