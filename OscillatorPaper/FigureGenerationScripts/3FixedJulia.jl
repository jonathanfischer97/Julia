begin 
    using CSV 
    using DataFrames 
    using DataFrameMacros
    using GLMakie; GLMakie.activate!()
    # using Plots 
    # using CairoMakie
end



# load CSV into DataFrame 
df = CSV.read("OscillatorPaper/FigureGenerationScripts/3FixedResultsCSVs/fixed_triplet_results-kcat1kcat7DF.csv", DataFrame)

testdf = testdf[(testdf.kcat1 .!= 0.0) .& (testdf.kcat7 .!= 0.0) .& (testdf.DF .!= 0.0), :]
df = testdf[log10.(testdf.kcat1) .> -10, :]




plot(testdf.kcat1, testdf.kcat7, testdf.DF, st = :scatter, xlabel = "kcat1", ylabel = "kcat7", zlabel="DF", title = "Average Period vs. Average Amplitude", legend = false)




# Convert to log scale and normalize
log_kcat1 = log10.(df.kcat1 .+ 1e-9) # Add epsilon to avoid log of zero
log_kcat7 = log10.(df.kcat7 .+ 1e-9)
log_DF = log10.(df.DF .+ 1e-9)

# Normalize sizes and colors
sizes = (df.average_amplitude .- minimum(df.average_amplitude)) ./ (maximum(df.average_amplitude) - minimum(df.average_amplitude)) * 10 .+ 5
colors = df.average_period
norm_colors = (colors .- minimum(colors)) ./ (maximum(colors) - minimum(colors))



#< GLMakie plots

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