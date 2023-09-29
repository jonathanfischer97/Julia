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




# using GLMakie; GLMakie.activate!()
import GLMakie as gm

function plot_3fixed_makie(df::DataFrame)

    xlog = log10.(df[:,1]) 
    ylog = log10.(df[:,2]) 
    zlog = log10.(df[:,3])

    xname, yname, zname = names(df)[1:3]
    dfval = df.DF[1]


    #* make the plot
    fig = gm.Figure(resolution = (1000, 600))

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
    
    ax = gm.Axis3(fig[1:3,1:3]; aspect=:data, perspectiveness=0.5, title="$xname vs. $yname vs $zname at DF = $dfval", xlabel = xname, ylabel = yname, zlabel = zname)

    # Scatter plot for non-NaN values
    hm = gm.meshscatter!(ax, xlog, ylog, zlog; markersize=sizes, ssao=true, color=df.average_period, colormap=:thermal, transparency=true, nan_color=:gray,
                        diffuse = gm.Vec3f(0.5, 0.5, 0.5), specular = gm.Vec3f(0.3, 0.3, 0.3), shininess = 100f0, ambient = gm.Vec3f(0.1), shading=true, alpha=0.7)

    # Colorbar and labels
    gm.Colorbar(fig[2, 4], hm, label="Period (s)", height=gm.Relative(2.0))
    gm.colgap!(fig.layout, 6)

    fig
end


plot_3fixed_makie(dfarray[3])


function plot_3fixed_makie(dfarray::Vector{DataFrame})

    xname, yname, zname = names(dfarray[1])[1:3]
    dfvals = [df.DF[1] for df in dfarray]


    #* make the figure
    fig = gm.Figure(resolution = (1200, 800))

    #* make 3 axes, one for each DF value
    axs = [gm.Axis3(fig[1,i]; aspect = :data, perspectiveness=0.5, title="$xname vs. $yname vs $zname at DF = $(dfvals[i])", xlabel = xname, ylabel = yname, zlabel = zname,            
                                    xtickformat = values -> [string(round(10^x; digits=2)) for x in values], #* converts the log10 values back to the original values
                                    ytickformat = values -> [string(round(10^x; digits=2)) for x in values],
                                    ztickformat = values -> [string(round(10^x; digits=2)) for x in values],) for i in 1:3]

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
        pl = gm.meshscatter!(ax, xlog[nonan_indices], ylog[nonan_indices], zlog[nonan_indices]; markersize=sizes[nonan_indices], ssao=true, color=df.average_period[nonan_indices], colormap=:thermal, transparency=true, nan_color=:gray,
                            diffuse = gm.Vec3f(0.5, 0.5, 0.5), specular = gm.Vec3f(0.3, 0.3, 0.3), shininess = 100f0, ambient = gm.Vec3f(0.1), shading=true, alpha=0.7)

        # Scatter with NaN values, high transparency
        gm.meshscatter!(ax, xlog[nan_indices], ylog[nan_indices], zlog[nan_indices]; markersize=0.1, ssao=true, colormap=:thermal, transparency=true,
                            diffuse = gm.Vec3f(0.5, 0.5, 0.5), specular = gm.Vec3f(0.3, 0.3, 0.3), shininess = 100f0, ambient = gm.Vec3f(0.1), shading=true, alpha=0.1)

        # Define shared attributes
        # shared_sizes = sizes[nonan_indices]
        # shared_colors = df.average_period[nonan_indices]

        #* Adding 2D projections
        # Projection on the x-y plane (z = minimum of zlog)
        # gm.scatter!(ax, xlog[nonan_indices], ylog[nonan_indices], fill(minimum(zlog), length(nonan_indices)); 
        #             color=shared_colors, marker=:circle, alpha=0.5)

        # # Projection on the x-z plane (y = minimum of ylog)
        # gm.scatter!(ax, xlog[nonan_indices], fill(minimum(ylog), length(nonan_indices)), zlog[nonan_indices]; 
        #             markersize=shared_sizes, color=shared_colors, marker=:circle, alpha=0.5)

        # # Projection on the y-z plane (x = minimum of xlog)
        # gm.scatter!(ax, fill(minimum(xlog), length(nonan_indices)), ylog[nonan_indices], zlog[nonan_indices]; 
        #             markersize=shared_sizes, color=shared_colors, marker=:circle, alpha=0.5)

        # Colorbar and labels
        gm.Colorbar(fig[i, i+1], pl, label="Period (s)", height=gm.Relative(1.0))
        gm.colgap!(fig.layout, 4)

        
    end

    # Colorbar and labels
    # gm.Colorbar(fig[1, 4], label="Period (s)", height=gm.Relative(1.0), colormap=:thermal)
    gm.colgap!(fig.layout, 6)

    fig
end


plot_3fixed_makie(dfarray)















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
    fig = gm.Figure(;resolution = (1000, 600))
    ax = gm.Axis3(fig[1,1]; aspect=:data, perspectiveness = 0.5)
    gm.scatter!(ax, x_coords, y_coords, z_coords, label="", xlabel="PC1", ylabel="PC2", title="PCA Scatter Plot with K-means Clusters")
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


