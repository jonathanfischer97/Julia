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
    
    using ColorSchemes, Plots.PlotMeasures
    # default(lw = 2, size = (1000, 600), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)


    using OscTools 

    import GLMakie: Figure, Axis3, meshscatter, meshscatter!, Colorbar, colgap!, rowgap!, Vec3f, Relative, save
    import GLMakie as gm

    using Base.Threads
    using Distances


    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))

    using FFTW
    FFTW.set_num_threads(18)
end

# function test4fixedGA()
#     #* Set up the default GA problem
#     ga_problem = GAProblem()

#     #* Fixed some constraints 
#     set_fixed_constraints!(ga_problem.constraints, [:DF, :K, :P, :A])

#     #* Assign the fixed values 
#     set_fixed_values!(ga_problem.constraints, 1000., 1.0, 1.0, 3.16)

#     #* Set seed 
#     Random.seed!(1234)

#     #* Generate the initial population
#     population = generate_population(ga_problem.constraints, 10000)

#     #* Run the GA
#     run_GA(ga_problem, population)
# end

ga_result = test4fixedGA()

population_to_matrix(ga_result)



#* Summary dataframes 
dfarray = read_csvs_in_directory("/home/local/WIN/jfisch27/Desktop/Julia/OscillatorPaper/FixedConstraintAnalysis/ROCKFISH_DATA/3Fixed/PopSize_15000/kcat1_kb2_P/SummaryResults")


"""
    plot_fixed_makie(df::DataFrame)
Plots 3D scatter plot of the dataframe for a single DF value using Makie. 
"""
function plot_fixed_makie(df::DataFrame)

    xlog = log10.(df[:,1]) 
    ylog = log10.(df[:,2]) 
    zlog = log10.(df[:,3])


    xname, yname, zname = names(df)[1:3]
    dfval = df.DF[1]

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
    sizes[nonan_indices] = ((nonan_numpoints .- minimum(nonan_numpoints)) ./ (maximum(nonan_numpoints) - minimum(nonan_numpoints))) ./ 2 

    # Normalize periods for non-NaN values
    norm_periods = fill(NaN, size(df, 1))
    norm_periods[nonan_indices] = (nonan_periods .- minimum(nonan_periods)) ./ (maximum(nonan_periods) - minimum(nonan_periods)) 

    # Create the figure and axis
    ax = Axis3(fig[1,1]; aspect=:data, perspectiveness=0.5, title="$xname vs. $yname vs $zname at DF = $dfval", 
                                                                xlabel = xname, ylabel = yname, zlabel = zname, xlabelfont = :bold, ylabelfont = :bold, zlabelfont=:bold)
    

    # Scatter plot for non-NaN values
    hm = meshscatter!(ax, xlog, ylog, zlog; markersize=sizes, ssao=true, color=df.average_period, colormap=:thermal, transparency=true, nan_color=:gray,
                        diffuse = Vec3f(0.1), specular = Vec3f(0.2), shininess = 100f0, ambient = Vec3f(0.1), shading=true, alpha=0.9, overdraw=true, backlight=1.0f0)


    #* Projection on the x-y plane (z = minimum of zlog)
    #Make Z-matrix for contour plot based on the period values, mapping each period value to the x-y coordinates

    # Use actual scatter plot data limits for defining heatmap position
    zmin = minimum(zlog[nonan_indices])
    # zmax = maximum(zlog[nonan_indices])
    xmax = maximum(xlog[nonan_indices])
    ymax = maximum(ylog[nonan_indices])

    # Projection on the x-y plane (z = minimum of zlog)
    xyplane, xzplane, yzplane = generate_z_matrices(df)

    # Define x, y, and z positions for the heatmap. They should match the actual scatter plot data points.
    x_positions = sort(unique(xlog[nonan_indices]))
    y_positions = sort(unique(ylog[nonan_indices]))
    z_positions = sort(unique(zlog[nonan_indices]))

    # Plot the heatmap for the x-y plane at z = zmin
    gm.heatmap!(ax, x_positions, y_positions, xyplane; transformation = (:xy, zmin + (zmin/10.)), colormap=:thermal, nan_color=:gray, alpha=0.9)

    # Plot the heatmap for the x-z plane at y = ymin
    gm.heatmap!(ax, x_positions, z_positions, xzplane; transformation = (:xz, ymax + (ymax/3.)), colormap=:thermal, nan_color=:gray, alpha=0.9)

    # Plot the heatmap for the y-z plane at x = xmin
    gm.heatmap!(ax, y_positions, z_positions, yzplane; transformation = (:yz, xmax + (xmax/3.)), colormap=:thermal, nan_color=:gray, alpha=0.9)

    # Colorbar and labels
    Colorbar(fig[1, 2], hm; label="Period (s)", labelfont=:bold)
    colgap!(fig.layout, 6)

    gm.hidespines!(ax)

    fig
end


fig = plot_fixed_makie(dfarray[3])
save("test_heat_projection.png", fig)

df = dfarray[3]
x = df[:,1]
x = reshape(x, 5, 5, 5)
x = reduce(+, x, dims=2)
x = dropdims(x, dims=3)
y = df[:,2]
z = df[:,3]
periods = df[:, :average_period]

z = reshape(periods, 5, 5, 5)
z = reduce(+, z, dims=3)
z = dropdims(z, dims=3)

gm.contour(z; transformation=(:xy, minimum(z)))
gm.heatmap(z)


# Step 1: Group the DataFrame by kcat1 and kb2
grouped = groupby(df, [:kcat1, :kb2])

# Step 2: Compute the average of :average_period for each group
reduced_df = combine(grouped, :average_period => mean)

# Extract unique values of kcat1 and kb2 for matrix dimensions
unique_kcat1 = unique(reduced_df.kcat1)
unique_kb2 = unique(reduced_df.kb2)

# Step 3: Transform the reduced :average_period values into a matrix
z_matrix = Matrix{Float64}(undef, length(unique_kcat1), length(unique_kb2))

for i in eachindex(unique_kcat1)
    for j in eachindex(unique_kb2)
        mask = (reduced_df.kcat1 .== unique_kcat1[i]) .& (reduced_df.kb2 .== unique_kb2[j])
        if any(mask)
            z_matrix[i, j] = reduced_df[mask, :average_period_mean][1]
        else
            z_matrix[i, j] = NaN
        end
    end
end

gm.contourf(z_matrix)



function generate_z_matrices(df::DataFrame)
    independent_vars = propertynames(df)[1:3]

    # Accessing raw arrays
    raw_data = Dict(var => df[!, var] for var in independent_vars)
    avg_period = df[!, :average_period]

    matrices = Vector{Matrix{Float64}}(undef, 3)

    idx = 1
    for i in 1:3
        for j in (i+1):3
            x_var, y_var = independent_vars[i], independent_vars[j]
            
            # Directly applying log10 transformation
            x_values = log10.(raw_data[x_var])
            y_values = log10.(raw_data[y_var])

            # Extract unique values for matrix dimensions
            unique_x = unique(x_values)
            unique_y = unique(y_values)

            x_map = Dict(val => index for (index, val) in enumerate(unique_x))
            y_map = Dict(val => index for (index, val) in enumerate(unique_y))

            # Initialize z_matrix and count_matrix
            z_matrix = fill(0.0, length(unique_x), length(unique_y))
            count_matrix = zeros(Int, length(unique_x), length(unique_y))

            # Iterate through the raw arrays
            for k in eachindex(x_values)
                x_idx = x_map[x_values[k]]
                y_idx = y_map[y_values[k]]

                # Accumulating values
                if !isnan(avg_period[k])
                    z_matrix[x_idx, y_idx] += avg_period[k]
                    count_matrix[x_idx, y_idx] += 1
                end
            end

            # Finalize the averaging process
            z_matrix .= z_matrix ./ count_matrix
            # replace!(z_matrix, NaN => 0.0)  # Replace NaNs with 0s

            # Store the z_matrix
            matrices[idx] = z_matrix
            idx += 1
        end
    end

    return (xy=matrices[1], xz=matrices[2], yz=matrices[3])
end



# Call the function
xy, xz, yz = generate_z_matrices(df)








"""
    plot_3fixed_makie(dfarray::Vector{DataFrame})
Plots three side by side 3D scatter plots of the dataframe for 3 values of DF using Makie.
"""
function plot_3fixed_makie(dfarray::Vector{DataFrame})

    xname, yname, zname = names(dfarray[1])[1:3]
    dfvals = [df.DF[1] for df in dfarray]


    #* make the figure
    fig = Figure(resolution = (1600, 600))

    #* make 3 axes, one for each DF value
    axs = [Axis3(fig[1,i]; aspect = :data, perspectiveness=0.5, title="$xname vs. $yname vs $zname at DF = $(dfvals[i])", xlabel = xname, ylabel = yname, zlabel = zname,            
                                    xtickformat = values -> [string(round(10^x; digits=2)) for x in values], #* converts the log10 values back to the original values
                                    ytickformat = values -> [string(round(10^x; digits=2)) for x in values],
                                    ztickformat = values -> [string(round(10^x; digits=2)) for x in values]) for i in 1:3]

    for (i,ax) in enumerate(axs)
        df = dfarray[i]

        xlog = log10.(df[:,1]) 
        ylog = log10.(df[:,2]) 
        zlog = log10.(df[:,3])


        #* make the scatter plot
        # Identify non-NaN indices and values
        nonan_indices = findall(!isnan, df[:, :average_period])

        nan_indices = findall(isnan, df[:, :average_period])
        # nonan_amplitudes = df[:, :average_amplitude][nonan_indices]
        nonan_numpoints = df[:, :num_oscillatory_points][nonan_indices]

        nonan_periods = df[:, :average_period][nonan_indices]

        # Normalize sizes for non-NaN values
        sizes = fill(0.1, size(df, 1))
        # sizes[nonan_indices] = ((nonan_amplitudes .- minimum(nonan_amplitudes)) ./ (maximum(nonan_amplitudes) - minimum(nonan_amplitudes))) ./ 2 
        sizes[nonan_indices] = ((nonan_numpoints .- minimum(nonan_numpoints)) ./ (maximum(nonan_numpoints) - minimum(nonan_numpoints))) ./ 2.5 

        # Normalize periods for non-NaN values
        norm_periods = fill(NaN, size(df, 1))
        norm_periods[nonan_indices] = (nonan_periods .- minimum(nonan_periods)) ./ (maximum(nonan_periods) - minimum(nonan_periods)) 

        # Scatter plot for non-NaN values
        pl = meshscatter!(ax, xlog[nonan_indices], ylog[nonan_indices], zlog[nonan_indices]; markersize=sizes[nonan_indices], ssao=true, color=df.average_period[nonan_indices], colormap=:thermal, transparency=true, nan_color=:gray,
                            diffuse = Vec3f(0.5, 0.5, 0.5), specular = Vec3f(0.3, 0.3, 0.3), shininess = 100f0, ambient = Vec3f(0.1), shading=true, alpha=0.7)

        # Scatter with NaN values, high transparency
        meshscatter!(ax, xlog[nan_indices], ylog[nan_indices], zlog[nan_indices]; markersize=0.1, ssao=true, colormap=:thermal, transparency=true,
                            diffuse = Vec3f(0.5, 0.5, 0.5), specular = Vec3f(0.3, 0.3, 0.3), shininess = 100f0, ambient = Vec3f(0.1), shading=true, alpha=0.1)

        # Define shared attributes
        # shared_sizes = sizes[nonan_indices]
        # shared_colors = df.average_period[nonan_indices]

        #* Adding 2D projections
        # Projection on the x-y plane (z = minimum of zlog)
        # scatter!(ax, xlog[nonan_indices], ylog[nonan_indices], fill(minimum(zlog), length(nonan_indices)); 
        #             color=shared_colors, marker=:circle, alpha=0.5)

        # # Projection on the x-z plane (y = minimum of ylog)
        # scatter!(ax, xlog[nonan_indices], fill(minimum(ylog), length(nonan_indices)), zlog[nonan_indices]; 
        #             markersize=shared_sizes, color=shared_colors, marker=:circle, alpha=0.5)

        # # Projection on the y-z plane (x = minimum of xlog)
        # scatter!(ax, fill(minimum(xlog), length(nonan_indices)), ylog[nonan_indices], zlog[nonan_indices]; 
        #             markersize=shared_sizes, color=shared_colors, marker=:circle, alpha=0.5)

        # Colorbar and labels
        Colorbar(fig[2, i], pl, label="Period (s)", vertical = false)#, height=Relative(1.0))
    end

    # Colorbar and labels
    colgap!(fig.layout, 6)
    rowgap!(fig.layout, 1)

    fig
end


fig = plot_3fixed_makie(dfarray)
save("test_glmakie.png", fig)


#< Loop through all the summary dataframes for each fixed combination and plot them
for dir in readdir("/home/local/WIN/jfisch27/Desktop/Julia/OscillatorPaper/FixedConstraintAnalysis/ROCKFISH_DATA/3Fixed/PopSize_15000", join=true)
    for subdir in readdir(dir, join=true)
        if basename(subdir) == "SummaryResults" && length(readdir(subdir)) == 3
            # println("Plotting $subdir")
            dfarray = read_csvs_in_directory(subdir)
            try
                fig = plot_3fixed_makie(dfarray)
                save("SUMMARY_SCATTERPLOTS/$(basename(dir))_summary_plot.png", fig)
            catch
                println("Error plotting $dir")
            end
        end
    end
end
#> END ###










#< CLUSTERING #####

rawdf_array = read_csvs_in_directory("OscillatorPaper/FixedConstraintAnalysis/ROCKFISH_DATA/3Fixed/PopSize_15000/kcat1_kb4_A/RawData/DF=1000.0")
combined_rawdf = vcat(rawdf_array...)

# log10.(combined_rawdf)



"""
    silhouette_score(X::AbstractMatrix{Float64}, labels::Vector{Int}, sample_size::Int=100)
Computes the silhouette score for a sample of data and labels.
"""
function silhouette_score(X::AbstractMatrix{Float64}, labels::Vector{Int}, sample_size::Int=100)
    # Sample data and corresponding labels
    idx = rand(1:size(X, 2), sample_size)
    sampled_X = X[:, idx]
    sampled_labels = labels[idx]
    
    dist_matrix = pairwise(Euclidean(), sampled_X, sampled_X)
    silhouettes = Clustering.silhouettes(sampled_labels, dist_matrix)
    
    return mean(silhouettes)
end

# Determine optimal cluster count using silhouette method
function optimal_kmeans_clusters(df::DataFrame, max_k::Int, exclude_cols::Vector{Symbol} = [])
    # Convert DataFrame to Matrix and exclude the specified columns
    data_matrix = df_to_matrix(df, exclude_cols)
    
    best_k = 2
    best_score = -Inf

    # Iterate through k values and compute silhouette score
    for k in 2:max_k
        result = kmeans(data_matrix, k)
        score = silhouette_score(data_matrix, assignments(result))
        if score > best_score
            best_score = score
            best_k = k
        end
    end
    
    return best_k
end




"""
    optimal_kmeans_clusters(data_matrix::AbstractMatrix{Float64}, max_k::Int)
Determine optimal cluster count using silhouette method with multithreading.
"""
function optimal_kmeans_clusters(data_matrix::AbstractMatrix{Float64}, max_k::Int)
    best_k = 2
    best_score = -Inf

    # Create an array to store the scores for each k
    scores = zeros(max_k - 1)

    # Iterate through k values and compute silhouette score in parallel
    @threads for k in 2:max_k
        result = kmeans(data_matrix, k)
        score = silhouette_score(data_matrix, assignments(result))
        scores[k - 1] = score # Store the score in the array
    end
    
    # Find the maximum score and the corresponding k value
    best_score, best_k = findmax(scores)
    best_k += 1 # Adjust the index to match the k value
    
    return best_k
end

"""
    get_optimal_clusters(df::DataFrame, max_k::Int, exclude_cols::Vector{Symbol} = [])
Wrapper function for optimal_kmeans_clusters that converts a DataFrame to a Matrix, and returns the optimal cluster count.
"""
function get_optimal_clusters(df::DataFrame, max_k::Int, exclude_cols::Vector{Symbol} = [])
    data_matrix = df_to_matrix(df, exclude_cols)
    return optimal_kmeans_clusters(data_matrix, max_k)
end


"""
    plot_3D_PCA_clusters(df::DataFrame)
Plots a 3D scatter plot of the log-transformed dataframe using PCA and K-means clustering.
"""
function plot_3D_PCA_clusters(df::DataFrame)
    # gm.activate!()
    # Step 1: Identify fixed columns
    fixed_cols = identify_fixed_columns(df)
    # @show fixed_cols

    # Step 2: Log transform the data prior to clustering, otherwise the clustering will be dominated by the largest values
    df = log10.(df)
    dfmat = df_to_matrix(df, fixed_cols)
    
    # Step 3: Dynamically determine the optimal number of clusters
    optimal_clusters = get_optimal_clusters(dfmat, 10) # Assuming a maximum of 10 clusters for demonstration
    
    
    # Step 4: Perform K-means clustering using optimal cluster count
    result = kmeans(dfmat, optimal_clusters)
   
    
    # Step 5: Perform PCA
    pca_model = fit(PCA, dfmat, maxoutdim=3)
    reduced_data = MultivariateStats.transform(pca_model, dfmat)
    reduced_centroids = MultivariateStats.transform(pca_model, result.centers)
    
    
    # Step 6: Extract 3D coordinates using new bases from PCA
    x_coords = (reduced_data[1, :])
    y_coords = (reduced_data[2, :])
    z_coords = (reduced_data[3, :])

    # Step 7: Map the centroids to the new bases
    centroid_x = reduced_centroids[1, :]
    centroid_y = reduced_centroids[2, :]
    centroid_z = reduced_centroids[3, :]
   
    
    # Step 8: Plotting
    fig = Figure(;resolution = (1000, 600))
    ax = Axis3(fig[1,1]; aspect=:data, perspectiveness = 0.5, title="PCA Scatter Plot with K-means Clusters (Log transformed data)")
    pc_variances = principalvars(pca_model)
    pc_percentages = pc_variances / sum(pc_variances) * 100
    ax.xlabel = "PC1 ($(round(pc_percentages[1], digits=2))%)"
    ax.ylabel = "PC2 ($(round(pc_percentages[2], digits=2))%)"
    ax.zlabel = "PC3 ($(round(pc_percentages[3], digits=2))%)"
    # ax = Axis3(fig[1,1]; aspect=:data, perspectiveness = 0.5, xlabel="PC1", ylabel="PC2", zlabel="PC3", title="PCA Scatter Plot with K-means Clusters (Log transformed data)")
    gm.scatter!(ax, x_coords, y_coords, z_coords, color = result.assignments, alpha=0.3)
    gm.scatter!(ax, centroid_x, centroid_y, centroid_z, marker=:x, label="Centroids", color=:red, markersize = 50, overdraw=true)

    return fig
end

fig = plot_3D_PCA_clusters(combined_rawdf)

save("test_pca_cluster.png", fig)







using Clustering, MultivariateStats
df = combined_rawdf
# Step 1: Identify fixed columns
fixed_cols = identify_fixed_columns(df)

# Step 2: Log transform the data prior to clustering, otherwise the clustering will be dominated by the largest values
df = log10.(df)
dfmat = df_to_matrix(df, fixed_cols)

# Step 3: Dynamically determine the optimal number of clusters
optimal_clusters = optimal_kmeans_clusters(dfmat, 10) # Assuming a maximum of 10 clusters for demonstration
@show optimal_clusters

# Step 4: Perform K-means clustering using optimal cluster count
result = kmeans(dfmat, optimal_clusters)
@show result.centers

# Step 5: Perform PCA
pca_model = fit(PCA, dfmat, maxoutdim=3)
reduced_data = MultivariateStats.transform(pca_model, dfmat)
reduced_centroids = MultivariateStats.transform(pca_model, result.centers)
@show reduced_centroids

# Step 6: Extract 3D coordinates using new bases from PCA
x_coords = (reduced_data[1, :])
y_coords = (reduced_data[2, :])
z_coords = (reduced_data[3, :])

# Step 7: Map the centroids to the new bases
centroid_x = reduced_centroids[1, :]
centroid_y = reduced_centroids[2, :]
centroid_z = reduced_centroids[3, :]
@show centroid_x, centroid_y, centroid_z

# Step 8: Plotting
fig = Figure(;resolution = (1000, 600))

ax = Axis3(fig[1,1]; aspect=:data, perspectiveness = 0.5, title="PCA Scatter Plot with K-means Clusters (Log transformed data)")
pc_variances = principalvars(pca_model)
pc_percentages = pc_variances / sum(pc_variances) * 100
ax.xlabel = "PC1 ($(round(pc_percentages[1], digits=2))%)"
ax.ylabel = "PC2 ($(round(pc_percentages[2], digits=2))%)"
ax.zlabel = "PC3 ($(round(pc_percentages[3], digits=2))%)"

gm.scatter!(ax, x_coords, y_coords, z_coords, color = result.assignments, alpha=0.1)
gm.scatter!(ax, centroid_x, centroid_y, centroid_z, marker=:x, label="Centroids", color=:black, markersize = 10)

fig









#< Loadings biplot

# Plot
import CairoMakie as cm
"""
    plot_loadings_biplot_2D(df::DataFrame)
Plots a 2D biplot of the loadings for the first two principal components.
"""
function plot_loadings_biplot_2D(df::DataFrame)
    cm.activate!()
    # Step 1: Identify fixed columns
    fixed_cols = identify_fixed_columns(df)

    # Step 2: Log transform the data prior to clustering, otherwise the clustering will be dominated by the largest values
    df = log10.(df)
    dfmat = df_to_matrix(df, fixed_cols)

    # Step 3: Perform PCA and extract loadings
    pca_model = fit(PCA, dfmat, maxoutdim=2)
    pcaloadings = loadings(pca_model)
    # Extract loadings for the first two principal components
    pc1_loadings = pcaloadings[:, 1]
    pc2_loadings = pcaloadings[:, 2]

    # Plot
    fig = cm.Figure(;resolution = (1000, 600))
    ax = cm.Axis(fig[1,1], title="PCA Loadings Biplot")

    # Get magnitude of loadings
    magnitude = sqrt.(pc1_loadings .^ 2 .+ pc2_loadings .^ 2)
    # colors = cm.colormap(:viridis, magnitude)

    quivplot = cm.quiver!(ax, zeros(size(pc1_loadings,)), zeros(size(pc2_loadings)), pc1_loadings, pc2_loadings, color = magnitude, colormap=:plasma, alpha=0.5)

    # Add labels and annotate arrows
    var_labels = names(df[:, Not(fixed_cols)])
    for (i, label) in enumerate(var_labels)
        cm.text!(ax, pc1_loadings[i], pc2_loadings[i]; text = label)#, halign=:center, valign=:center, rotation=atan(pc2_loadings[i] / pc1_loadings[i]) * 180 / pi)
    end

    # Label axes with percentage of variance explained
    pc_variances = principalvars(pca_model)
    pc_percentages = pc_variances / sum(pc_variances) * 100
    ax.xlabel = "PC1 ($(round(pc_percentages[1], digits=2))%)"
    ax.ylabel = "PC2 ($(round(pc_percentages[2], digits=2))%)"

    cm.Colorbar(fig[1,2], quivplot, label="Loadings Magnitude")

    fig
end

biplot = plot_loadings_biplot(combined_rawdf)
cm.save("test_biplot2D.png", biplot)


"""
    plot_loadings_biplot_3D(df::DataFrame)
Plots a 3D biplot of the loadings for the first three principal components.
"""
function plot_loadings_biplot_3D(df::DataFrame)
    # Step 1: Identify fixed columns
    fixed_cols = identify_fixed_columns(df)

    # Step 2: Log transform the data prior to clustering, otherwise the clustering will be dominated by the largest values
    df = log10.(df)
    dfmat = df_to_matrix(df, fixed_cols)

    # Step 3: Perform PCA and extract loadings
    pca_model = fit(PCA, dfmat, maxoutdim=3)
    pcaloadings = loadings(pca_model)
    # Extract loadings for the first three principal components
    pc1_loadings = pcaloadings[:, 1]
    pc2_loadings = pcaloadings[:, 2]
    pc3_loadings = pcaloadings[:, 3]

    # Plot
    gm.with_theme(gm.theme_dark()) do
        fig = Figure(;resolution = (1000, 600))
        ax = Axis3(fig[1,1], title="PCA Loadings Biplot")

        # Get magnitude of loadings
        magnitude = sqrt.(pc1_loadings .^ 2 .+ pc2_loadings .^ 2 .+ pc3_loadings .^ 2)

        origin = (zeros(size(pc1_loadings,)), zeros(size(pc2_loadings)), zeros(size(pc3_loadings)))

        quivplot = gm.arrows!(ax, origin..., pc1_loadings, pc2_loadings, pc3_loadings; color = magnitude, colormap=:inferno, alpha=0.1, arrowsize=0.1)

        # Add labels and annotate arrows
        var_labels = names(df[:, Not(fixed_cols)])
        for (i, label) in enumerate(var_labels)
            gm.text!(ax, pc1_loadings[i], pc2_loadings[i], pc3_loadings[i]; text = label, overdraw=true, font = :bold,fontsize=20, color=:white)#, halign=:center, valign=:center, rotation=atan(pc2_loadings[i] / pc1_loadings[i]) * 180 / pi)
        end

        # Label axes with percentage of variance explained
        pc_variances = principalvars(pca_model)
        pc_percentages = pc_variances / sum(pc_variances) * 100
        ax.xlabel = "PC1 ($(round(pc_percentages[1], digits=2))%)"
        ax.ylabel = "PC2 ($(round(pc_percentages[2], digits=2))%)"
        ax.zlabel = "PC3 ($(round(pc_percentages[3], digits=2))%)"

        Colorbar(fig[1,2], quivplot, label="Loadings Magnitude")

        fig
    end
end

biplot = plot_loadings_biplot_3D(combined_rawdf)
save("test_biplot3D.png", biplot)












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


