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

    import GLMakie: Figure, Axis3, meshscatter, meshscatter!, Colorbar, colgap!, rowgap!, Vec3f, Relative, save, with_theme
    import GLMakie as gm

    using Base.Threads
    using Distances


    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))

    using FFTW
    FFTW.set_num_threads(18)
end



ga_result = test4fixedGA()
population_to_matrix(ga_result)






#* Summary dataframes 
dfarray = read_csvs_in_directory("/home/local/WIN/jfisch27/Desktop/Julia/OscillatorPaper/FixedConstraintAnalysis/ROCKFISH_DATA/3Fixed/PopSize_15000/L_K_P/SummaryResults")

for df in dfarray
    df.period_range = df.maximum_period .- df.minimum_period
end

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
    gm.heatmap!(ax, x_positions, y_positions, xyplane; transformation = (:xy, zmin + (zmin/3.)), colormap=:thermal, nan_color=:gray, alpha=0.9)

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




"""
    generate_z_matrices(df::DataFrame)
Generates the z-matrices for the x-y, x-z, and y-z planes for a given DataFrame, for either contour or heatmaps.
"""
function generate_z_matrices(df::DataFrame, zvar=:average_period)
    independent_vars = propertynames(df)[1:3]

    # Accessing raw arrays
    raw_data = Dict(var => df[!, var] for var in independent_vars)
    zvec = df[!, zvar]

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
                if !isnan(zvec[k])
                    z_matrix[x_idx, y_idx] += zvec[k]
                    count_matrix[x_idx, y_idx] += 1
                end
            end

            # Finalize the averaging process
            z_matrix .= z_matrix ./ count_matrix

            # Store the z_matrix
            matrices[idx] = z_matrix
            idx += 1
        end
    end

    return (xy=matrices[1], xz=matrices[2], yz=matrices[3])
end

# Call the function
xy, xz, yz = generate_z_matrices(dfarray[1])

function generate_z_matrices(xvec, yvec, zvec, colorvec)

    # Directly applying log10 transformation to the columns and subset unique
    x_values = log10.(xvec) |> unique
    y_values = log10.(yvec) |> unique
    z_values = log10.(zvec) |> unique


    # Initialize a 3D matrix
    data_cube = fill(0.0, length(x_values), length(y_values), 3)
    
    # Populate the 3D matrix




    # Compute averages along each dimension of the 3D matrix
    xy = dropdims(mean(data_cube, dims=3), dims=3)
    xz = dropdims(mean(data_cube, dims=2), dims=2)
    yz = dropdims(mean(data_cube, dims=1), dims=1)

    return (xy=xy, xz=xz, yz=yz)
end

xy, xz, yz = generate_z_matrices(dfarray[1])




"""
    get_new_label(oldlabel)
Converts the old label to a new label for plotting.
"""
function get_new_label(oldlabel)
    labelmap = Dict("ka1" => "k⁽ᴸᴷ⁾", "kb1" => "kᵣ⁽ᴸᴷ⁾", "kcat1" => "kcat⁽ᴸᴷ⁾", 
                    "ka2" => "k⁽ᴸᴬ⁾", "kb2" => "kᵣ⁽ᴸᴬ⁾",
                    "ka3" => "k⁽ᴬᴷ)", "kb3" => "kᵣ⁽ᴬᴷ⁾",
                    "ka4" => "k⁽ᴬᴾ⁾", "kb4" => "kᵣ⁽ᴬᴾ⁾",
                    "ka7" => "k⁽ᴸᴾ⁾", "kb7" => "kᵣ⁽ᴸᴾ⁾", "kcat7" => "kcat⁽ᴸᴾ⁾",
                    "L" => "L", "K" => "K", "P" => "P", "A" => "A")

    return labelmap[oldlabel]
end


function add_units(labels...)
    unitlabels = []
    for name in labels
        if occursin("kᵣ", name)
            push!(unitlabels,name * " (s⁻¹)")
        elseif occursin("kcat", name)
            push!(unitlabels,name * " (s⁻¹)")
        elseif occursin("k", name)
            push!(unitlabels,name * " (μM⁻¹s⁻¹)")
        else
            push!(unitlabels,name * " (μM)")
        end
    end
    return unitlabels
end



"""
    plot_3fixed_makie(dfarray::Vector{DataFrame}, colorvar = :average_period)
Uses GLMakie to plot three side by side 3D scatter plots for each dataframe in array, each with a different value of DF.
"""
function plot_3fixed_makie(dfarray::Vector{DataFrame}, colorvar = :average_period)

    xlabel, ylabel, zlabel = get_new_label.(names(dfarray[1])[1:3])
    dfvals = [df.DF[1] for df in dfarray]

    xlabel_units, ylabel_units, zlabel_units = add_units(xlabel, ylabel, zlabel)

    # Determine global color limits
    clims = (minimum([minimum(filter(!isnan,df[:, colorvar])) for df in dfarray]),
    maximum([maximum(filter(!isnan, df[:, colorvar])) for df in dfarray]))

    # Determine global size limits
    minsize, maxsize = (minimum([minimum(filter(!iszero,df[:, :num_oscillatory_points])) for df in dfarray]),
    maximum([maximum(filter(!iszero, df[:, :num_oscillatory_points])) for df in dfarray]))

    # Color map
    # cmap = gm.Reverse(:berlin)
    cmap = :redsblues


    #* make the figure
    fig = with_theme(theme_black()) do
        fig = Figure(resolution = (2000, 800))

        #* make 3 axes, one for each DF value
        axs = [Axis3(fig[1,i]; aspect = :data, perspectiveness=0.5, title="($xlabel vs. $ylabel vs $zlabel) at DF = $(dfvals[i])", xlabel = "log $xlabel_units", ylabel = "log $ylabel_units", zlabel = "log $zlabel_units", xlabelfont = :bold, ylabelfont = :bold, zlabelfont=:bold, 
                                        titlesize=25, xlabelsize = 22, ylabelsize=22, zlabelsize=22,  
                                        xticks = unique(log10.(dfarray[i][:,1])), yticks = unique(log10.(dfarray[i][:,2])), zticks = unique(log10.(dfarray[i][:,3])),        
                                        xtickformat = values -> [string(round(10^x; digits=3)) for x in values], #* converts the log10 values back to the original values
                                        ytickformat = values -> [string(round(10^x; digits=3)) for x in values],
                                        ztickformat = values -> [string(round(10^x; digits=3)) for x in values]) for i in 1:3]

        


        for (i,ax) in enumerate(axs)
            df = deepcopy(dfarray[i])

            # Log transform the data
            xlog = log10.(df[:,1]) 
            ylog = log10.(df[:,2]) 
            zlog = log10.(df[:,3])

            # Replace any outliers with NaN
            if occursin("amplitude", string(colorvar))
                outlier_indices = findall(x -> x > 300., df[:, colorvar])
                df[outlier_indices, colorvar] .= NaN
            end



            #* make the scatter plot
            # Identify NaN indices 
            nan_indices = findall(isnan, df[:, colorvar])

            # Non-NaN indices just all indices that are not in nan_indices
            nonan_indices = setdiff(1:size(df, 1), nan_indices)

            # Number of oscillatory points controls size
            nonan_numpoints = df[:, :num_oscillatory_points][nonan_indices]
            
            # nonan_periods = df[:, :average_period][nonan_indices]
            # nonan_amplitudes = df[:, :average_amplitude][nonan_indices]

            # Normalize sizes for non-NaN values
            sizes = fill(0.1, size(df, 1))
            sizes[nonan_indices] = 0.1 .+ (((nonan_numpoints .- minsize) ./ (maxsize - minsize)) ./ 3.9) 


            # Scatter plot for non-NaN values
            pl = meshscatter!(ax, xlog[nonan_indices], ylog[nonan_indices], zlog[nonan_indices]; markersize=sizes[nonan_indices], ssao=true, color=df[:, colorvar][nonan_indices], colormap=cmap, transparency=false,
                                diffuse = Vec3f(0.7), specular = Vec3f(0.3), shininess = 100f0, ambient = Vec3f(0.1), shading=true, alpha=0.9, overdraw=false)#, backlight=1.0f0)

            # Scatter with NaN values, high transparency
            meshscatter!(ax, xlog[nan_indices], ylog[nan_indices], zlog[nan_indices]; markersize=0.05, ssao=true, color=:white, transparency=true,
                                        diffuse = Vec3f(0.6), specular = Vec3f(0.2), shininess = 100f0, ambient = Vec3f(0.1), shading=true, alpha=0.4)

            #* 2D projection heatmaps
            # Use actual scatter plot data limits for defining heatmap position
            zmin = minimum(zlog)
            xmax = maximum(xlog)
            ymax = maximum(ylog)

            # Projection on the x-y plane (z = minimum of zlog)
            xyplane, xzplane, yzplane = generate_z_matrices(df, colorvar)

            # Define x, y, and z positions for the heatmap. They should match the actual scatter plot data points.
            x_positions = sort(unique(xlog))
            y_positions = sort(unique(ylog))
            z_positions = sort(unique(zlog))

            # Plot the heatmap for the x-y plane at z = zmin
            heatmap_args = (colormap=cmap, nan_color=:black, alpha=0.8)
            gm.heatmap!(ax, x_positions, y_positions, xyplane; transformation = (:xy, zmin + (zmin/4.)), heatmap_args...)

            # Plot the heatmap for the x-z plane at y = ymin
            gm.heatmap!(ax, x_positions, z_positions, xzplane; transformation = (:xz, ymax + (ymax/4.)), heatmap_args...)

            # Plot the heatmap for the y-z plane at x = xmin
            gm.heatmap!(ax, y_positions, z_positions, yzplane; transformation = (:yz, xmax + (xmax/4.)), heatmap_args...)

            gm.hidespines!(ax)
            # gm.hidedecorations!(ax; label=false, ticklabels=false, ticks=false)

            # Colorbar and labels
            # Colorbar(fig[2, i], pl, label=colorlabel, vertical = false)
        end

        colorlabel = string(colorvar)
        if occursin("period", colorlabel)
            colorlabel *= " (s)"
        elseif occursin("amplitude", colorlabel)
            colorlabel *= " (μM)"
        end

        Colorbar(fig[2, 1:3]; label=colorlabel, vertical = false, limits=clims, colormap= cmap)
        colgap!(fig.layout, 6)
        rowgap!(fig.layout, 1)
        fig
    end
    fig
end


fig = plot_3fixed_makie(dfarray, :period_range)
save("test_glmakie.png", fig)


df = dfarray[3]
findall(isnan, df[:, :average_period])
findall(x -> x <= 1, df[:, :num_oscillatory_points])


#< Loop through all the summary dataframes for each fixed combination and plot them
for dir in readdir("/home/local/WIN/jfisch27/Desktop/Julia/OscillatorPaper/FixedConstraintAnalysis/ROCKFISH_DATA/3Fixed/PopSize_15000", join=true)
    for subdir in readdir(dir, join=true)
        if basename(subdir) == "SummaryResults" && length(readdir(subdir)) == 3
            # println("Plotting $subdir")
            dfarray = read_csvs_in_directory(subdir)

            for df in dfarray
                df.period_range = df.maximum_period .- df.minimum_period
            end

            try
                colorvar = :period_range
                fig = plot_3fixed_makie(dfarray, colorvar)
                pathstring = "./SUMMARY_SCATTERPLOTS_HEATMAPS_" * string(colorvar)
                mkpath(pathstring)
                save("$pathstring/$(basename(dir))_summary_plot.png", fig)
            catch e
                println(e)
            end
        end
    end
end
#> END ###





gm.set_theme!()




#< CLUSTERING #####

rawdf_array = read_csvs_in_directory("OscillatorPaper/FixedConstraintAnalysis/ROCKFISH_DATA/3Fixed/PopSize_15000/L_K_A/RawData/DF=1000.0")
combined_rawdf = vcat(rawdf_array...)

# log10.(combined_rawdf)

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
    # df = log10.(ogdf[!,Not(fixed_cols)])
    # @show df
    dfmat = log10.(df_to_matrix(df, fixed_cols))
    
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
    ax = Axis3(fig[1,1]; aspect=:data, perspectiveness = 0.5, title="PCA (Log transformed data)")
    pc_variances = principalvars(pca_model)
    pc_percentages = pc_variances / sum(pc_variances) * 100
    ax.xlabel = "PC1 ($(round(pc_percentages[1], digits=2))%)"
    ax.ylabel = "PC2 ($(round(pc_percentages[2], digits=2))%)"
    ax.zlabel = "PC3 ($(round(pc_percentages[3], digits=2))%)"
    # ax = Axis3(fig[1,1]; aspect=:data, perspectiveness = 0.5, xlabel="PC1", ylabel="PC2", zlabel="PC3", title="PCA Scatter Plot with K-means Clusters (Log transformed data)")
    sp = gm.scatter!(ax, x_coords, y_coords, z_coords, color = df.per, alpha=0.5)
    # gm.scatter!(ax, centroid_x, centroid_y, centroid_z, marker=:x, label="Centroids", color=:red, markersize = 50, overdraw=true)

    Colorbar(fig[1, 2]; label="Period (s)", limits=(minimum(df.per), maximum(df.per)))

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
    # fixed_cols = fixed_cols[[1,2,3,5]]
    

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









using GLMakie
function ODE_movie(sol)
    points = Observable(Point3f[]) # Signal that can be used to update plots efficiently
    colors = Observable(Int[])

    set_theme!(theme_black())

    fig, ax, l = lines(points, color = colors,
        colormap = :inferno, transparency = true, 
        axis = (; type = Axis3, protrusions = (0, 0, 0, 0), 
                viewmode = :fit, limits = (1, 4, 0, 1, 0, 2)))
    
    Amem = sol[6,:] .+ sol[9,:] .+ sol[10,:] .+ sol[11,:] .+ sol[12,:] .+ sol[15,:] .+ sol[16,:]
    # @show typeof(Amem)

    record(fig, "odetest.mp4", 1:120) do frame
        for i in eachindex(ogsol.t)
            # update arrays inplace
            push!(points[], Point3f(sol.u[i][1], sol.u[i][2], sol.u[i][6] + sol.u[i][9] + sol.u[i][10] + sol.u[i][11] + sol.u[i][12] + sol.u[i][15] + sol.u[i][16]))
            push!(colors[], frame)
        end
        ax.azimuth[] = 1.7pi + 0.3 * sin(2pi * frame / 120) # set the view angle of the axis
        notify(points); notify(colors) # tell points and colors that their value has been updated
        l.colorrange = (0, frame) # update plot attribute directly
    end
end

ogprobjac = make_ODE_problem()
ogsol = solve(ogprobjac, Rodas5P(), saveat=0.1)
ODE_movie(ogsol)








using Plots


"""Converts a rate constant from nm^3/us) to 1/(uM*s)"""
function convert_to_macrorate(microrate::Float64)
    return microrate * 0.602214076
end
psym = [:ka1 => convert_to_macrorate(0.027267670545107203)#5.453534109021441e-05 #ka1, 1
    :kb1 => 0.0643048008980449 #kb1, 2
    :kcat1 => 286.6382253193995 #kcat1, 3
    :ka2 => convert_to_macrorate(1.6605390671738467) #0.0033210781343476934 #ka2, 4
    :kb2 => 0.39569337786534897 #kb2, 5
    :ka3 => convert_to_macrorate(0.04115484270564225) #ka3, 6
    :kb3 => 0.5393197910059361 #kb3, 7
    :ka4 => convert_to_macrorate(0.05448532670350475) #ka4, 8
    :kb4 => 0.2897657637531564 #kb4, 9
    :ka7 => convert_to_macrorate(0.19013582091039463) #ka7, 10
    :kb7 => 0.0028126177315505618 #kb7, 11
    :kcat7 => 1.2733781341040291 #kcat7, 12
    :DF => 25500.] #DF, 13
p = [x[2] for x in psym]


# usym = [:L => 634.6315289074139, :K => 47.42150049582334, :P => 239.66312964177104,  :A => 838.7724702072363, :Lp => 0.0, :LpA => 0.0, :LK => 0.0, #:Lp => 790.5014385747756,
#         :LpP => 0.0, :LpAK => 0.0, :LpAP => 0.0, :LpAKL => 0.0, :LpAPLp => 0.0, :AK => 0.0, :AP => 0.0, :AKL => 0.0, :APLp => 0.0]
usym = [:L => 2.0, :K => 0.05, :P => 0.55,  :A => 2.0, :Lp => 0.0, :LpA => 0.0, :LK => 0.0, #:Lp => 790.5014385747756,
        :LpP => 0.0, :LpAK => 0.0, :LpAP => 0.0, :LpAKL => 0.0, :LpAPLp => 0.0, :AK => 0.0, :AP => 0.0, :AKL => 0.0, :APLp => 0.0]
u0 = [x[2] for x in usym]

newprob = remake(ogprobjac, p=p, u0=u0, tspan= (0.0, 100.0))
newsol = solve(newprob, Rosenbrock23(), saveat=0.1)
plot(newsol)

gaprob = GAProblem()
allconstraints = AllConstraints()
population = generate_population(allconstraints, 20000)
garesults = run_GA(gaprob, population)
gadf = make_ga_dataframe(garesults, allconstraints)

function get_prob(dfrow::DataFrameRow, prob::ODEProblem)
    newp = [param for param in dfrow[Between(:ka1, :DF)]]
    newu0 = [ic for ic in dfrow[Between(:L, :A)]]

    remake(prob, p = newp, u0 = [newu0; zeros(length(prob.u0) - length(newu0))])
end

probarray = []
for row in eachrow(gadf)
    push!(probarray, get_prob(row, ogprobjac))
end
probarray


import CairoMakie as cm
fig = cm.Figure(resolution = (1000, 600))
ax = cm.Axis(fig[1,1], title="ODE Solution")
cm.lines!(ax, newsol[1,:], newsol[3,:], color=:blue, label="L")
fig

function make_tuneplot(prob)
    fig = cm.Figure(resolution = (1000, 600))
    ax = cm.Axis(fig[1,1], title="ODE Solution")


    params = copy(prob.p)
    for i in 1:6
        params[end] = params[end] * 1.1
        newprob = remake(prob, p=params)
        newsol = solve(newprob, Rosenbrock23(), saveat=0.1)
        cm.lines!(ax, newsol[1,:], newsol[3,:])
    end
    fig
end
    
make_tuneplot(probarray[3])


testsol = solve(probarray[1500], Rosenbrock23(), saveat=0.1)
plot(testsol)