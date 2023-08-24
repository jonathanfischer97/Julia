begin 
    using CSV 
    using DataFrames 
    using DataFrameMacros
    using DifferentialEquations
    using Catalyst
    using Evolutionary, FFTW
    using Random
    using Distributions
    using ProgressMeter
    using CSV

    using LinearAlgebra

    using ColorSchemes#, PlotMeasures

    using Combinatorics

    # default(lw = 2, size = (1000, 600), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)


    include("../../UTILITIES/EvolutionaryOverloads.jl")

    # import the Catalyst model "fullrn"
    include("../../UTILITIES/ReactionNetwork.jl")

    # import the cost function and other evaluation functions
    include("../../UTILITIES/EvaluationFunctions.jl")
    # using .EvaluationFunctions

    # import the genetic algorithm and associated functions
    include("../../UTILITIES/GA_functions.jl")

    include("../../UTILITIES/TestBenchPlotUtils.jl")


    # using GLMakie; GLMakie.activate!(ssao=true)
    using Plots, Plots.PlotMeasures; default(lw = 2, size = (1000, 600), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)
    # using CairoMakie; CairoMakie.activate!()

    # const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))
end


pnames = ["ka1","kb1","kcat1"]

# load CSV into DataFrame 
df = CSV.read("OscillatorPaper/FigureGenerationScripts/3FixedResultsCSVs/fixed_triplet_results-$(pnames[1])$(pnames[2])$(pnames[3]).csv", DataFrame)
excessL_df = CSV.read("OscillatorPaper/FigureGenerationScripts/3FixedResultsCSVsExcessL/fixed_triplet_results-$(pnames[1])$(pnames[2])$(pnames[3]).csv", DataFrame)


#* read in test data for fixed ka1, kcat1, kcat7 and initial conditions freely optimized
testfig_df = CSV.read("OscillatorPaper/FigureGenerationScripts/test.csv", DataFrame)

fullrn = make_fullrn()

#* Get the solution to a row of the DataFrame
function get_solution_from_row(dfrow::DataFrameRow, rn)
    #* get the parameters
    p = [dfrow.ka1, dfrow.kb1, dfrow.kcat1, dfrow.ka2, dfrow.kb2, dfrow.ka3, dfrow.kb3, dfrow.ka4, dfrow.kb4, dfrow.ka7, dfrow.kb7, dfrow.kcat7, dfrow.DF]
    #* get the initial conditions
    u0 = [[dfrow.L, dfrow.K, dfrow.P, dfrow.A]; zeros(12)]
    #* get the solution
    sol = solve(ODEProblem(rn, u0, (0.0,2000.0), p), Rosenbrock23(), saveat = 0.1)
    return sol
end

sol = get_solution_from_row(testfig_df[1,:], fullrn)


fig = Figure(resolution=(1600, 1000))




ezra_df = CSV.read("OscillatorPaper/FigureGenerationScripts/AllExpOsc .csv", DataFrame)
Acol = log10.(Vector(ezra_df.A))
Kcol = log10.(Vector(ezra_df.K))
Pcol = log10.(Vector(ezra_df.P))
Lcol = log10.(Vector(ezra_df.L))

function create_3d_scatter_with_shadows_ezra(x,y,z,colors; angle=30)
    # Convert to log scale and add epsilon to avoid log of zero
    # pnames = names(df)[1:3]

    # Identify non-NaN indices
    # nonan_indices = findall(!isnan, df[:, :average_period])
    # nonan_amplitudes = df[:, :average_amplitude][nonan_indices]
    # nonan_periods = df[:, :average_period][nonan_indices]

    # Normalize sizes for non-NaN values
    sizes = fill(5., length(x))
    # sizes[nonan_indices] = (nonan_amplitudes .- minimum(nonan_amplitudes)) ./ (maximum(nonan_amplitudes) - minimum(nonan_amplitudes)) * 10 .+ 5

    # Normalize periods for non-NaN values
    # norm_periods_nonan = (nonan_periods .- minimum(nonan_periods)) ./ (maximum(nonan_periods) - minimum(nonan_periods))

    # Create a color gradient for non-NaN values
    color_gradient = cgrad(:viridis, colors, rev=true)

    # Main plot for non-NaN values
    p = scatter3d(x,y,z, color=color_gradient, legend=false, alpha=1.0, markerstrokewidth=0)

    # Plot NaN values in gray
    # nan_indices = findall(isnan, df[:, :average_period])
    # scatter3d!(p, x[nan_indices], y[nan_indices], z[nan_indices], markersize=sizes[nan_indices], color=:gray, legend=false, alpha=1.0, markerstrokewidth=0)

    # Logic for handling NaN values
    # for (xi, yi, zi, si, ci) in zip(x, y, z, sizes, colors)
    #     # color_proj = isnan(ci) ? :gray : colormap[ci]
    #     scatter3d!(p,[xi], [yi], [minimum(z)], markersize=[si], color=color_proj, alpha=0.4, markerstrokewidth=0)
    #     scatter3d!(p,[xi], [yi], [minimum(y)], [zi], markersize=[si], color=color_proj, alpha=0.4, markerstrokewidth=0)
    #     scatter3d!(p,[minimum(x)], [yi], [zi], markersize=[si], color=color_proj, alpha=0.4, markerstrokewidth=0)
    # end

    # Project the points onto each plane and connect them with lines
    for (xi, yi, zi, si, ci) in zip(x, y, z, sizes, color_gradient)
        color_proj = color_gradient
        scatter3d!(p,[xi], [yi], [minimum(z)], markersize=[si], color=color_proj, alpha=0.4, markerstrokewidth=0)
        scatter3d!(p,[xi], [minimum(y)], [zi], markersize=[si], color=color_proj, alpha=0.4, markerstrokewidth=0)
        scatter3d!(p,[minimum(x)], [yi], [zi], markersize=[si], color=color_proj, alpha=0.4, markerstrokewidth=0)
        plot3d!(p,[xi, xi], [yi, yi], [zi, minimum(z)], color=color_proj, lw=0.5, linealpha=0.5, linestyle=:dash)
        plot3d!(p,[xi, xi], [yi, minimum(y)], [zi, zi], color=color_proj, lw=0.5, linealpha=0.5, linestyle=:dash)
        plot3d!(p,[xi, minimum(x)], [yi, yi], [zi, zi], color=color_proj, lw=0.5, linealpha=0.5, linestyle=:dash)
    end

    # Set axis labels
    xlabel!(p,"log10(A)")
    ylabel!(p,"log10(K)")
    zlabel!(p,"log10(P)")

    # Set view angle
    plot!(p,camera=(30, angle))

    # Set title
    title!(p,"Log-Scaled Parameters with Projections & Shadows (Angle: $angle) - Size & Color: Average Amplitude & Period")
    display(p)
end

# Sample call to the modified function
create_3d_scatter_with_shadows_ezra(Acol, Kcol, Pcol, Lcol, angle=30)

ezfig, ax, pl = meshscatter(Acol, Kcol, Pcol; color=Lcol, colormap=:thermal, markersize=0.1, ssao=true, transparency=false)

Colorbar(ezfig, label="L", height=Relative(2.0))
ezfig
xlabel!(ax3, "log10(kb3)")
ylabel!(ax3, "log10(kb4)")
zlabel!(ax3, "log10(DF)")



# # Convert to log scale and normalize
# log_kcat1 = log10.(df[:,pnames[1]] .+ 1e-9) # Add epsilon to avoid log of zero
# log_kcat7 = log10.(df[:,pnames[2]] .+ 1e-9)
# log_DF = log10.(df[:,pnames[3]] .+ 1e-9)

# nonan_amplitudes = df.average_amplitude[.!isnan.(df.average_amplitude)]
# nonan_periods = df.average_period[.!isnan.(df.average_period)]

# # Normalize sizes and colors
# sizes = (nonan_amplitudes .- minimum(nonan_amplitudes)) ./ (maximum(nonan_amplitudes) - minimum(nonan_amplitudes)) * 10 .+ 5
# colors = nonan_periods
# norm_colors = (colors .- minimum(colors)) ./ (maximum(colors) - minimum(colors))

# numpoints = df.num_oscillatory_points
# sizes = (numpoints .- minimum(numpoints)) ./ (maximum(numpoints) - minimum(numpoints))

#< GLMakie plots
# using Makie

function scatters_in_3D(df; fig = Figure(resolution=(1600, 1000)),  angle=30)
    pnames = names(df)[1:3]
    xlog = log10.(df[:, pnames[1]])
    ylog = log10.(df[:, pnames[2]])
    zlog = log10.(df[:, pnames[3]])

    x = df[:, pnames[1]]
    y = df[:, pnames[2]]
    z = df[:, pnames[3]]

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
    
    ax = Axis3(fig[1:2,3]; aspect=:data, perspectiveness=0.5, title="3 Fixed Parameter Oscillatory Regions Excess L", xlabel = pnames[1], ylabel = pnames[2], zlabel = pnames[3])

    # Scatter plot for non-NaN values
    hm = meshscatter!(ax, xlog, ylog, zlog; markersize=sizes, ssao=true, color=df.average_period, colormap=:thermal, transparency=false, nan_color=:gray,
                        diffuse = Vec3f(0.5, 0.5, 0.5), specular = Vec3f(0.3, 0.3, 0.3), shininess = 100f0, ambient = Vec3f(0.1), shading=true)

    # meshscatter!(ax, xlog, ylog; markersize=sizes, marker= Rect3f(Vec3f(0.,0.,0.1)), ssao=true, color=df.average_period, colormap=:thermal, transparency=false, nan_color=:gray,
    #                     diffuse = Vec3f(0.0), specular = Vec3f(0.0), shininess = 0, ambient = Vec3f(0.0))

    # Scatter plot for NaN values in gray
    # nan_indices = findall(isnan, df[:, :average_period])
    # meshscatter!(ax3, x[nan_indices], y[nan_indices], z[nan_indices]; markersize=sizes[nan_indices], color=:gray)

    # Colorbar and labels
    # Colorbar(fig[2, 4], hm, label="Period (s)", height=Relative(2.0))
    colgap!(fig.layout, 5)
    # xlabel!(ax3, "log10(kb3)")
    # ylabel!(ax3, "log10(kb4)")
    # zlabel!(ax3, "log10(DF)")

    # Display plot
    fig
end

fig = scatters_in_3D(df)

fig2 = scatters_in_3D(excessL_df; fig=fig)


GLMakie.save("3FixedParameterOscillatoryRegionsExcessLCompare.png", fig2)



function scatters_in_3D(x, y, z, sizes, norm_colors)
    n = length(x)
    aspect=(1)
    perspectiveness=0.5
    # the figure
    fig = Figure(; resolution=(800, 800))
    # ax1 = Axis3(fig[1, 1]; aspect, perspectiveness)
    # ax2 = Axis3(fig[1, 2]; aspect, perspectiveness)
    ax3 = Axis3(fig[1,2]; aspect=:data, perspectiveness)
    # scatter!(ax1, x, y, z; markersize=15)
    # meshscatter!(ax2, x, y, z; markersize=0.25)
    hm = meshscatter!(ax3, x, y, z; markersize=0.25, ssao = true,
        marker=Rect3f(Vec3f(0), Vec3f(1)), color=norm_colors,
        colormap=:plasma, transparency=false)
    Colorbar(fig[1, 4], hm, label="values", height=Relative(0.5))
    colgap!(fig.layout, 5)
    fig
end

scatters_in_3D(log_kcat1, log_kcat7, log_DF, sizes, colors)

using Plots

function create_3d_scatter_with_shadows(df; angle=30)
    # Convert to log scale and add epsilon to avoid log of zero
    pnames = names(df)[1:3]
    x = log10.(df[:, pnames[1]])
    y = log10.(df[:, pnames[2]])
    z = log10.(df[:, pnames[3]])

    # Identify non-NaN indices
    nonan_indices = findall(!isnan, df[:, :average_period])
    nonan_amplitudes = df[:, :average_amplitude][nonan_indices]
    nonan_periods = df[:, :average_period][nonan_indices]

    # Normalize sizes for non-NaN values
    sizes = fill(5., size(df, 1))
    sizes[nonan_indices] = (nonan_amplitudes .- minimum(nonan_amplitudes)) ./ (maximum(nonan_amplitudes) - minimum(nonan_amplitudes)) * 10 .+ 5

    # Normalize periods for non-NaN values
    norm_periods_nonan = (nonan_periods .- minimum(nonan_periods)) ./ (maximum(nonan_periods) - minimum(nonan_periods))

    # Create a color gradient for non-NaN values
    color_gradient = cgrad(:viridis, norm_periods_nonan, rev=true)

    # Main plot for non-NaN values
    p = scatter3d(x[nonan_indices], y[nonan_indices], z[nonan_indices], markersize=sizes[nonan_indices], color=color_gradient[norm_periods_nonan], legend=false, alpha=1.0, markerstrokewidth=0)

    # Plot NaN values in gray
    nan_indices = findall(isnan, df[:, :average_period])
    scatter3d!(p, x[nan_indices], y[nan_indices], z[nan_indices], markersize=sizes[nan_indices], color=:gray, legend=false, alpha=1.0, markerstrokewidth=0)

    # Logic for handling NaN values
    # for (xi, yi, zi, si, ci) in zip(x, y, z, sizes, colors)
    #     # color_proj = isnan(ci) ? :gray : colormap[ci]
    #     scatter3d!(p,[xi], [yi], [minimum(z)], markersize=[si], color=color_proj, alpha=0.4, markerstrokewidth=0)
    #     scatter3d!(p,[xi], [yi], [minimum(y)], [zi], markersize=[si], color=color_proj, alpha=0.4, markerstrokewidth=0)
    #     scatter3d!(p,[minimum(x)], [yi], [zi], markersize=[si], color=color_proj, alpha=0.4, markerstrokewidth=0)
    # end

    # Project the points onto each plane and connect them with lines
    for (xi, yi, zi, si, ci) in zip(x, y, z, sizes, color_gradient)
        color_proj = color_gradient
        scatter3d!(p,[xi], [yi], [minimum(z)], markersize=[si], color=color_proj, alpha=0.4, markerstrokewidth=0)
        scatter3d!(p,[xi], [minimum(y)], [zi], markersize=[si], color=color_proj, alpha=0.4, markerstrokewidth=0)
        scatter3d!(p,[minimum(x)], [yi], [zi], markersize=[si], color=color_proj, alpha=0.4, markerstrokewidth=0)
        plot3d!(p,[xi, xi], [yi, yi], [zi, minimum(z)], color=color_proj, lw=0.5, linealpha=0.5, linestyle=:dash)
        plot3d!(p,[xi, xi], [yi, minimum(y)], [zi, zi], color=color_proj, lw=0.5, linealpha=0.5, linestyle=:dash)
        plot3d!(p,[xi, minimum(x)], [yi, yi], [zi, zi], color=color_proj, lw=0.5, linealpha=0.5, linestyle=:dash)
    end

    # Set axis labels
    xlabel!(p,"log10($(pnames[1]))")
    ylabel!(p,"log10($(pnames[2]))")
    zlabel!(p,"log10($(pnames[3]))")

    # Set view angle
    plot!(p,camera=(30, angle))

    # Set title
    title!(p,"Log-Scaled Parameters with Projections & Shadows (Angle: $angle) - Size & Color: Average Amplitude & Period")
    display(p)
end

# Sample call to the modified function
create_3d_scatter_with_shadows(df)


color_gradient = cgrad(:viridis, reverse(nonan_periods); rev=true)

color_gradient[[1,2]]

# Plotting function
#* Plotting function
function create_3d_scatter_with_shadows(angle, x, y, z, sizes, norm_colors)
    # Colormap
    colormap = cgrad(:viridis, norm_colors, rev=true)
    p = scatter3d(x, y, z, markersize=sizes, color=colormap[norm_colors], legend=false, alpha=1.0, markerstrokewidth=0)

    # Project the points onto each plane and connect them with lines
    for (xi, yi, zi, si, ci) in zip(x, y, z, sizes, norm_colors)
        color_proj = colormap[ci]
        scatter3d!(p,[xi], [yi], [minimum(z)], markersize=[si], color=color_proj, alpha=0.4, markerstrokewidth=0)
        scatter3d!(p,[xi], [minimum(y)], [zi], markersize=[si], color=color_proj, alpha=0.4, markerstrokewidth=0)
        scatter3d!(p,[minimum(x)], [yi], [zi], markersize=[si], color=color_proj, alpha=0.4, markerstrokewidth=0)
        plot3d!(p,[xi, xi], [yi, yi], [zi, minimum(z)], color=color_proj, lw=0.5, linealpha=0.5, linestyle=:dash)
        plot3d!(p,[xi, xi], [yi, minimum(y)], [zi, zi], color=color_proj, lw=0.5, linealpha=0.5, linestyle=:dash)
        plot3d!(p,[xi, minimum(x)], [yi, yi], [zi, zi], color=color_proj, lw=0.5, linealpha=0.5, linestyle=:dash)
    end

    # Set axis labels
    xlabel!(p,"log10(kcat1)")
    ylabel!(p,"log10(kcat7)")
    zlabel!(p,"log10(DF)")

    # Set view angle
    plot!(p,camera=(30, angle))

    # Set title
    title!(p,"Log-Scaled Parameters with Projections & Shadows (Angle: $angle) - Size & Color: Average Amplitude & Period")

    display(p)
end

# Example usage
create_3d_scatter_with_shadows(30, log_kcat1, log_kcat7, log_DF, sizes, norm_colors)






# Plot from different angles with log scaling and projections
angles = [0, 30, 60, 90, 120, 150]
for angle in angles
    create_3d_scatter_with_shadows(angle, log_kcat1, log_kcat7, log_DF, sizes, colors)
end



testplot = scatter3d(df.kcat1, df.kcat7, df.DF, color=cgrad(:viridis, colors, rev=true)[colors], xscale = :log10, yscale = :log10)