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



#* Clustering test 
dfarray = read_csvs_in_directory("OscillatorPaper/FigureGenerationScripts/ROCKFISH_DATA/4Fixed/PopSize_10000/4FixedICRawSets/DF=1000.0")



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

fig, ax, plt = meshscatter(x_coords, y_coords, z_coords, zcolor=result.assignments, colorbar=true, label="", xlabel="PC1", ylabel="PC2", title="PCA Scatter Plot with K-means Clusters")

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



#* Make figure for Maggie from Ezra's model 

"""
Converts row in oscillatory solution dataframe to an ODESolution
- `df` Dataframe of oscillatory solutions
- `row` Row of dataframe to solve
"""
function entryToSol(df, row)
    ka7Range = 10.0 .^((2:5) ./ 5)
    Km1exp = 15.0
    ka2exp = 0.7*10^-3
    kb2exp = 2.0 * 10^-3
    ka3exp = 0.118
    kb3exp = 0.0609
    Kd4exp =29.0
    Km7exp = 39.0
    kcat7exp = 85.3


    currow = df[row,:]
    u0 = zeros(16)
    u0[1] = currow[:L]
    u0[2] = currow[:K]
    u0[3] = currow[:P]
    u0[4] = currow[:A]
    ka1est = currow[:ka1]
    kb1est = currow[:kb1]
    ka4est = currow[:ka4]
    ka7est = currow[:ka7]
    dfest = currow[:df]
    psym = [:ka1 => ka1est, :kb1 => kb1est, :kcat1 => Km1exp * ka1est - kb1est, :ka2 => ka2exp, :kb2 => kb2exp,
                            :ka3 => ka3exp, :kb3 => kb3exp, :ka4 => ka4est, :kb4 => Kd4exp * ka4est, :ka7 => ka7est, 
                            :kb7 => Km7exp * ka7est - kcat7exp, :kcat7exp => 85.3, :y => dfest]
    p = [x[2] for x in psym]
    # return solve(remake(prob, u0=u0, tspan=(0.0, tspan), p=p), Rodas4(), abstol=1e-8, reltol=1e-12, saveat=0.1, save_idxs=1, maxiters=200 * tspan, verbose=false, callback=TerminateSteadyState(1e-8, 1e-12))
    return u0, p 
end


ezdf = CSV.read("/home/local/WIN/jfisch27/Desktop/Julia/OscillatorPaper/EzraSolutions.csv", DataFrame)


oprob = make_ODE_problem()

slow_newu, slow_newp = entryToSol(ezdf, 5)
slowprob = remake(oprob, u0=slow_newu, p=slow_newp, tspan=(0.,2000.))
slowsol = solve(slowprob, Rodas5P(), saveat=0.1)

fast_newu, fast_newp = entryToSol(ezdf, 3)
fastprob = remake(oprob, u0=fast_newu, p=fast_newp, tspan=(0.,2000.))
fastsol = solve(fastprob, Rodas5P(), saveat=0.1)

using Plots

function apply_default_settings(p)
    plot!(p, lw=4, size=(1000, 600), dpi=200,
          bottom_margin=12px, left_margin=16px, top_margin=10px, right_margin=8px)
    return p
end

function plotsol_maggie(sol::ODESolution; title = "")

    #* Sum up all A in solution 
    Asol = sol[4,:] .+ sol[13,:] .+ sol[14, :]

    #* Sum up all A on membrane
    Amem = OscTools.calculate_Amem(sol)
    
    p = plot(sol, idxs = [1,5,2,3], title = title, xlabel = "Time (s)", ylabel = "Concentration (µM)",
            color = [:blue :orange :purple :gold], label = ["PIP" "PIP2" "PIP5K" "Synaptojanin"], alpha = 0.5, lw=3)

    # p = plot(xlabel = "Time (s)", ylabel = "Concentration (µM)")

    plot!(p, sol.t, Asol, label="AP2 in solution", ls = :dash, alpha=1.0, color=:gray, lw=3)
    plot!(p, sol.t, Amem, label = "AP2 on membrane", ls = :dash, alpha=1.0, color=:black, lw=3)

    return p |> apply_default_settings
end


plotsol_maggie(fastsol; title= "Fast")

import CairoMakie as cm 

publication_theme() = cm.Theme(
    fontsize=16,
    font="CMU Serif",
    cm.Axis=(
        xlabelsize=20,xlabelpadding=-5,
        xgridstyle=:dash, ygridstyle=:dash,
        xtickalign=1, ytickalign=1,
        yticksize=10, xticksize=10,
        xlabel="x", ylabel="y"
        ),
    cm.Legend=(framecolor=(:black, 0.5), bgcolor=(:white, 0.5)),
    cm.Colorbar=(ticksize=16, tickalign=1, spinewidth=0.5),
)

function maggie_plot(slowsol, fastsol)
    fig = cm.Figure(;resolution = (1000, 600))

    ax = cm.Axis(fig[1,1]; title= "Membrane binding of AP2",  xlabel = "Time (s)", ylabel = "Concentration (µM)",
    titlesize = 25, xlabelsize = 22, ylabelsize = 22, xgridstyle=:dash, ygridstyle=:dash, xtickalign=1, ytickalign=1)

    cm.lines!(ax, slowsol.t, OscTools.calculate_Amem(slowsol), color=:blue, linewidth=3, label = "Period: 9.5 s")

    cm.lines!(ax, fastsol.t, OscTools.calculate_Amem(fastsol), color=:tomato, linewidth=3, label = "Period: 91 s")

    cm.axislegend(ax)

    fig
end

maggie_plot(slowsol, fastsol)


f, ax, l1 = cm.lines(slowsol.t, OscTools.calculate_Amem(slowsol);
                                                    figure = (; resolution = (1000, 600)),

                                                    axis = (; title= "Membrane binding of AP2", xlabel = "Time (s)", ylabel = "Concentration (µM)", xlabelsize = 20, ylabelsize = 20, ticklabelsize = 20, legendfontsize = 20),

                                                    color=:blue, linewidth=3, label = "Period: 9.5 s")
cm.lines!(ax, fastsol.t, OscTools.calculate_Amem(fastsol), color=:tomato, linewidth=3, label = "Period: 91 s")
cm.axislegend(ax)
f