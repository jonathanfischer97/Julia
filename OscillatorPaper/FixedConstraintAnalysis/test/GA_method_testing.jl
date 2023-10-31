begin 
    using Plots
    using Catalyst
    using OrdinaryDiffEq, ModelingToolkit
    using Statistics

    using Random
    using Distributions
    using DataFrames
    using CSV
    using StaticArrays
    using BenchmarkTools, ProgressMeter


    using LinearAlgebra
    using MultivariateStats
    using Clustering

    
    using ColorSchemes, Plots.PlotMeasures
    # default(lw = 2, size = (1000, 600), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)


    using OscTools 


    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))

    using FFTW

    numthreads = Threads.nthreads()
    println("Threads detected: $numthreads")
    numcores = min(1,numthreads÷2)
    BLAS.set_num_threads(numcores)
    FFTW.set_num_threads(numcores)
end


#* First testing different group sizes for tournament selection

ga_df = test4fixedGA(10000, Symbol[], Float64[]; num_tournament_groups=5)

"""Tournament selection with unique winners using BitArray"""
function unique_tournament_bitarray(groupSize::Int; select=argmax)
    @assert groupSize > 0 "Group size must be positive"
    function tournamentN(fitness::AbstractVecOrMat{<:Real}, N::Int;
                         rng::AbstractRNG=Random.GLOBAL_RNG)
        sFitness = size(fitness)
        d, nFitness = length(sFitness) == 1 ? (1, sFitness[1]) : sFitness
        selected_flags = falses(nFitness)  # BitArray
        selection = Vector{Int}(undef, N)
        
        count = 1
        while count <= N
            tour = randperm(rng, nFitness)
            j = 1
            while (j+groupSize) <= nFitness && count <= N
                idxs = tour[j:j+groupSize-1]
                idxs = filter(x -> !selected_flags[x], idxs)  # Remove already selected
                
                if isempty(idxs)
                    j += groupSize
                    continue
                end
                
                selected = d == 1 ? view(fitness, idxs) : view(fitness, :, idxs)
                winner = select(selected)
                winner_idx = idxs[winner]
                
                if !selected_flags[winner_idx]
                    selection[count] = winner_idx
                    selected_flags[winner_idx] = true
                    count += 1
                end
                
                j += groupSize
            end
        end
        
        # If needed, fill the rest of the selection with random individuals
        while count <= N
            rand_idx = rand(1:nFitness)
            if !selected_flags[rand_idx]
                selection[count] = rand_idx
                count += 1
            end
        end
        
        return selection
    end
    return tournamentN
end

"""Tournament selection with unique winners using BitArray, modified to select N_selected individuals"""
function unique_tournament_bitarray(groupSize::Int; select=argmax)
    @assert groupSize > 0 "Group size must be positive"
    function tournamentN(fitness::AbstractVecOrMat{<:Real}, N_selected::Int;
                         rng::AbstractRNG=Random.GLOBAL_RNG)
        sFitness = size(fitness)
        d, nFitness = length(sFitness) == 1 ? (1, sFitness[1]) : sFitness
        selected_flags = falses(nFitness)  # BitArray
        selection = Vector{Int}(undef, N_selected)
        
        count = 1
        while count <= N_selected
            tour = randperm(rng, nFitness)
            j = 1
            while (j+groupSize) <= nFitness && count <= N_selected
                idxs = tour[j:j+groupSize-1]
                idxs = filter(x -> !selected_flags[x], idxs)  # Remove already selected
                
                if isempty(idxs)
                    j += groupSize
                    continue
                end
                
                selected = d == 1 ? view(fitness, idxs) : view(fitness, :, idxs)
                winner = select(selected)
                winner_idx = idxs[winner]
                
                if !selected_flags[winner_idx]
                    selection[count] = winner_idx
                    selected_flags[winner_idx] = true
                    count += 1
                end
                
                j += groupSize
            end
        end
        
        return selection
    end
    return tournamentN
end

tourfunc = unique_tournament_bitarray(10)
fitvals = zeros(100)
fitvals[[4,9,50]] .= rand(3)

selected_idxs = tourfunc(fitvals, 90)

oldtourfunc = OscTools.tournament(10)
selected_idxs = oldtourfunc(fitvals, 90)



"""Tournament selection with unique winners using BitArray"""
function unique_tournament_bitarray(groupSize::Int; select=argmax)
    @assert groupSize > 0 "Group size must be positive"

    function tournamentN(fitness::AbstractVecOrMat{<:Real}, N::Int, n_newInds::Int;
                         rng::AbstractRNG=Random.GLOBAL_RNG)
        # Initialize dimensions and BitArray to track selected individuals
        sFitness = size(fitness)
        d, nFitness = length(sFitness) == 1 ? (1, sFitness[1]) : sFitness
        selected_flags = falses(nFitness)

        # Initialize output array for selected individual indices
        N_selected = N - n_newInds  # Number of individuals to be selected
        selection = zeros(Int, N_selected)

        count = 1  # Counter for unique individuals

        # Loop to fill up selection array with unique winners
        for _ in N_selected
            # Randomly permute population indices for the tournament
            tour = randperm(rng, nFitness)

            # Loop through the permuted indices in chunks of 'groupSize'
            for k in 1:groupSize:(nFitness - groupSize + 1)
                idxs = filter(x -> !selected_flags[tour[x]], k:(k + groupSize - 1))
                idxs = tour[idxs]

                # Skip if all are already selected
                if isempty(idxs)
                    continue
                end

                # Get the winner and update BitArray and selection vector
                selected = d == 1 ? view(fitness, idxs) : view(fitness, :, idxs)
                @info "Length of idxs: ", length(idxs)
                @info "Output of select(selected): ", select(selected)
                winner_idx = idxs[select(selected)]
                selection[count] = winner_idx
                selected_flags[winner_idx] = true
                count += 1
                
            end
        end
        return selection  # Return only the required unique individuals
    end
    return tournamentN
end




ga_df = test4fixedGA(100000, Symbol[], Float64[]; num_tournament_groups=10, selection_method=OscTools.unique_tournament_bitarray, n_newInds=0.1)
ga_df = test4fixedGA(50000, [:DF], [1000.]; num_tournament_groups=10, selection_method=OscTools.tournament, n_newInds=0.95)

allconstraints = AllConstraints()
testpop = generate_population(allconstraints, 1000)








"""
    coverage_metric(population::AbstractMatrix{Float64})

Calculate the coverage of the population in the parameter space.
Each parameter's range is calculated, and the sum of these ranges is returned.
"""
function get_coverage(population::AbstractMatrix{Float64}, constraints::CT) where CT <: ConstraintSet
    # Log Transformation and Normalization as before
    # log_population = log.(population)
    log_population = population
    min_vals = minima(constraints)
    max_vals = maxima(constraints)
    normalized_population = (log_population .- min_vals) ./ (max_vals - min_vals)

    # Calculate the range of each parameter
    param_ranges = maximum(normalized_population, dims=2) - minimum(normalized_population, dims=2)
    # @info "Parameter ranges: ", param_ranges

    # Sum the ranges to get the overall coverage metric
    return sum(param_ranges)
end

function get_coverage(df::DataFrame, constraints::CT, exclude_cols::Vector{Symbol} = [:gen, :fit, :per, :amp, :relamp]) where CT <: ConstraintSet
    dfmat = df_to_matrix(df, exclude_cols)
    return get_coverage(dfmat, constraints)
end


get_coverage(ga_df, allconstraints)


# Function to calculate entropy
function entropy(values::Vector{Float64})
    counts = Dict{Float64, Int}()
    for v in values
        counts[v] = get(counts, v, 0) + 1
    end
    probs = [count/length(values) for count in values(counts)]
    H = -sum(p * log2(p) for p in probs)
    return H
end



log_ga_df = log.(ga_df[:, Not([:gen, :fit, :per, :amp, :relamp])])

ga_df2 = copy(ga_df)

using DataFrames
transform!(ga_df2, Not([:gen, :fit, :per, :amp, :relamp]) .=> ByRow(log), renamecols=false)



"""
    calculate_generation_metrics(gen_df::DataFrame)

Calculate various metrics for a single generation of a Genetic Algorithm.
"""
function calculate_generation_metrics(gen_df)
    excluded_cols = [:gen, :fit, :per, :amp, :relamp, :DF]

    # Log transform parameter columns
    log_gen_df = copy(gen_df)
    transform!(log_gen_df, Not(excluded_cols) .=> ByRow(log), renamecols=false)


    # Max Pairwise Distance Diversity of Parameters
    max_pairwise_distance = getmax_pairwise_diversity(log_gen_df, excluded_cols)

    # Get the max, min, and mean of periods
    period_max = maximum(log_gen_df.per)
    period_min = minimum(log_gen_df.per)
    period_mean = mean(log_gen_df.per)

    # Get the max, min, and mean of amplitudes
    amplitude_max = maximum(log_gen_df.relamp)
    amplitude_min = minimum(log_gen_df.relamp)
    amplitude_mean = mean(log_gen_df.relamp)

    # Get the mean fitness value
    mean_fitness_val = mean(log_gen_df.fit)

    # Get the number of clusters
    cluster_number = get_optimal_clusters(log_gen_df, 20, excluded_cols)

    # Get the max and mean cluster distances
    result = get_kmeans(log_gen_df, cluster_number, excluded_cols)
    cluster_distance_matrix = get_cluster_distances(result)
    max_cluster_distance = maximum(cluster_distance_matrix)
    mean_cluster_distance = mean(cluster_distance_matrix)
    

    return max_pairwise_distance, period_max, period_min, period_mean, amplitude_max, amplitude_min, amplitude_mean, mean_fitness_val, cluster_number, max_cluster_distance, mean_cluster_distance
end


"""
    calculate_ga_metrics(df::DataFrame)

Calculate various metrics to assess the performance and diversity of a Genetic Algorithm.
"""
function calculate_ga_metrics(df::DataFrame)
    # Initialize arrays to store metrics
    generations = unique(df.gen)
    num_points = Vector{Int}(undef, length(generations))
    max_pairwise_distances = Vector{Float64}(undef, length(generations))
    period_maxes = similar(max_pairwise_distances)
    period_mins = similar(max_pairwise_distances)
    period_means = similar(max_pairwise_distances)
    amplitude_maxes = similar(max_pairwise_distances)
    amplitude_mins = similar(max_pairwise_distances)
    amplitude_means = similar(max_pairwise_distances)
    mean_fitness_vals = similar(max_pairwise_distances)
    number_of_clusters = similar(num_points)
    max_cluster_distances = similar(max_pairwise_distances)
    mean_cluster_distances = similar(max_pairwise_distances)


    # Calculate metrics for each generation
    for (i,gen) in enumerate(generations)
        gen_df = df[df.gen .== gen, :]

        max_pairwise_distance, period_max, period_min, mean_period, amplitude_max, amplitude_min, mean_amplitude, mean_fitness_val, cluster_number, max_cluster_distance, mean_cluster_distance = calculate_generation_metrics(gen_df)

        num_points[i] = nrow(gen_df)
        max_pairwise_distances[i] = max_pairwise_distance
        period_maxes[i] = period_max
        period_mins[i] = period_min
        period_means[i] = mean_period
        amplitude_maxes[i] = amplitude_max
        amplitude_mins[i] = amplitude_min
        amplitude_means[i] = mean_amplitude
        mean_fitness_vals[i] = mean_fitness_val
        number_of_clusters[i] = cluster_number
        max_cluster_distances[i] = max_cluster_distance
        mean_cluster_distances[i] = mean_cluster_distance
    end
    return (generations = generations, num_points=num_points, max_pairwise_distances = max_pairwise_distances, period_maxes = period_maxes, period_mins = period_mins, period_means = period_means, amplitude_maxes = amplitude_maxes, amplitude_mins = amplitude_mins, amplitude_means = amplitude_means, mean_fitness_vals = mean_fitness_vals, number_of_clusters = number_of_clusters, max_cluster_distances = max_cluster_distances, mean_cluster_distances = mean_cluster_distances)
end

"""
    calculate_cumulative_ga_metrics(df::DataFrame)

Calculate various cumulative metrics to assess the performance and diversity of a Genetic Algorithm.
"""
function calculate_cumulative_ga_metrics(df::DataFrame)
    # Initialize arrays to store metrics
    generations = unique(df.gen)
    cum_num_points = Vector{Int}(undef, length(generations))
    cum_maxpairwise_distances = Vector{Float64}(undef, length(generations))
    cum_period_maxes = similar(cum_maxpairwise_distances)
    cum_period_mins = similar(cum_maxpairwise_distances)
    cum_period_means = similar(cum_maxpairwise_distances)
    cum_amplitude_maxes = similar(cum_maxpairwise_distances)
    cum_amplitude_mins = similar(cum_maxpairwise_distances)
    cum_amplitude_means = similar(cum_maxpairwise_distances)
    cum_mean_fitness_vals = similar(cum_maxpairwise_distances)
    cum_number_of_cluster = similar(cum_num_points)
    cum_max_cluster_distances = similar(cum_maxpairwise_distances)
    cum_mean_cluster_distances = similar(cum_maxpairwise_distances)

    # Calculate cumulative metrics for each generation
    for (i, gen) in enumerate(generations)
        gen_df = df[df.gen .<= gen, :]

        cum_max_pairwise_distance, cum_period_max, cum_period_min, cum_period_mean, cum_amplitude_max, cum_amplitude_min, cum_amplitude_mean, cum_mean_fitness_val, cum_cluster_number, cum_max_cluster_distance, cum_mean_cluster_distance = calculate_generation_metrics(gen_df)

        cum_num_points[i] = nrow(gen_df)
        cum_maxpairwise_distances[i] = cum_max_pairwise_distance
        cum_period_maxes[i] = cum_period_max
        cum_period_mins[i] = cum_period_min
        cum_period_means[i] = cum_period_mean
        cum_amplitude_maxes[i] = cum_amplitude_max
        cum_amplitude_mins[i] = cum_amplitude_min
        cum_amplitude_means[i] = cum_amplitude_mean
        cum_mean_fitness_vals[i] = cum_mean_fitness_val
        cum_number_of_cluster[i] = cum_cluster_number
        cum_max_cluster_distances[i] = cum_max_cluster_distance
        cum_mean_cluster_distances[i] = cum_mean_cluster_distance
    end
    return (generations = generations, cum_num_points=cum_num_points, cum_maxpairwise_distances = cum_maxpairwise_distances, cum_period_maxes = cum_period_maxes, cum_period_mins = cum_period_mins, cum_period_means = cum_period_means, cum_amplitude_maxes = cum_amplitude_maxes, cum_amplitude_mins = cum_amplitude_mins, cum_amplitude_means = cum_amplitude_means, cum_mean_fitness_vals = cum_mean_fitness_vals, cum_number_of_cluster = cum_number_of_cluster, cum_max_cluster_distances = cum_max_cluster_distances, cum_mean_cluster_distances = cum_mean_cluster_distances)
end


"""
    plot_ga_metrics(df::DataFrame)

Plot various metrics to assess the performance and diversity of a Genetic Algorithm.
"""
function plot_ga_metrics(df::DataFrame; cum = true)
    metrics = cum ? calculate_cumulative_ga_metrics(df) : calculate_ga_metrics(df)

    generations, max_pairwise_distances, period_maxes, period_mins, period_means, amplitude_maxes, amplitude_mins, amplitude_means, mean_fitness_vals, number_of_clusters, max_cluster_distances, mean_cluster_distances = metrics

    # Create the plots
    p1 = plot(generations, max_pairwise_distances, label="", xlabel="Generation", ylabel="Max Pairwise Distance", title="Parameter Diversity Over Generations", color=:blue)
    p2 = plot(generations, period_maxes, label="Max", xlabel="Generation", ylabel="Seconds", title="Period Over Generations", color=:green)
    plot!(p2, generations, period_mins, label="Min", color=:green, linestyle=:dash)
    plot!(p2, generations, period_means, label="Mean", color=:green, linestyle=:dot)
    p3 = plot(generations, amplitude_maxes, label="Max", xlabel="Generation", ylabel="% of Initial AP2", title="Relative Amplitude Over Generations", color=:orange)
    plot!(p3, generations, amplitude_mins, label="Min", color=:orange, linestyle=:dash)
    plot!(p3, generations, amplitude_means, label="Mean", color=:orange, linestyle=:dot)
    p4 = plot(generations, mean_fitness_vals, label="", xlabel="Generation", ylabel="Fitness", title="Mean Fitness Over Generations", color=:purple)
    p5 = plot(generations, number_of_clusters, label="", xlabel="Generation", ylabel="N Optimal Clusters", title="Optimal Clusters N Over Generations", color=:brown)
    p6 = plot(generations, max_cluster_distances, label="Max", xlabel="Generation", ylabel="Max Cluster Distances", title="Max & Mean Cluster Distances Over Generations", color=:crimson)
    plot!(p6, generations, mean_cluster_distances, label="Mean", color=:crimson, linestyle=:dash)

    # Combine the plots into a single figure
    suptitle = cum ? "Cumulative " : "Non-cumulative"

    plot(p1, p2, p3, p4, p5, p6, suptitle=suptitle, layout=(3, 2), legend=:bottomright, size = (1200, 800), lw=4, dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)
end


# Generate the plot
plot_ga_metrics(ga_df; cum=true)
savefig("pop50000_20groups_0.95newinds_metrics_CUM.png")

for n_groups in 5:5:20
    ga_df = test4fixedGA(50000, [:DF], [1000.]; num_tournament_groups=n_groups, selection_method=OscTools.tournament, n_newInds=0.90)
    plot_ga_metrics(ga_df; cum=true)
    savefig("pop50000_$(n_groups)groups_0.90newinds_DF1000_metrics_CUM.png")
end




popsizes = []
dfs = []
for popsize in [1000, 10000, 50000, 100000]
    ga_df = test4fixedGA(popsize, [:DF], [10000.]; num_tournament_groups=10, selection_method=OscTools.tournament, n_newInds=0.90)
    push!(popsizes, popsize)
    push!(dfs, ga_df)
end

group_sizes = []
dfs = []
for group_size in [5, 10, 20, 50]
    ga_df = test4fixedGA(50000, [:DF], [10000.]; num_tournament_groups=group_size, selection_method=OscTools.tournament, n_newInds=0.90)
    push!(group_sizes, group_size)
    push!(dfs, ga_df)
end
group_size_dfs = dfs

n_newInds_array = []
new_inds_dfs = []
for n_newInds in [0.5, 0.75, 0.9, 0.95]
    ga_df = test4fixedGA(50000, [:DF], [10000.]; num_tournament_groups=10, selection_method=OscTools.tournament, n_newInds=n_newInds)
    push!(n_newInds_array, n_newInds)
    push!(dfs, ga_df)
end

groupsize_fig = plot_comparative_ga_metrics(group_size_dfs, group_sizes, "Number of Tournament Groups")
newInds_fig = plot_comparative_ga_metrics(new_inds_dfs, n_newInds_array, "Fractional Size of Scouting Party")

cm.save("n_tourneygroups_comparison.svg", groupsize_fig)
cm.save("n_newInds_comparison.svg", newInds_fig)



# """
#     plot_comparative_ga_metrics(dfarray::Vector{DataFrame}, hyperparam_vals::Vector{<:Real}, hyperparam_name; cum = true)

# Plot various metrics between GA runs with different hyperparameters (ex. popsize, group_size).
# """
# function plot_comparative_ga_metrics(dfarray::Vector, hyperparam_vals::Vector, hyperparam_name; cum = true)

#     #Initialize subplots 
#     num_points_plot = plot(title = "Number of Points", xlabel = "Generation", ylabel = "Number of Points", label = "")
#     pairwise_plot = plot(title = "Max Pairwise Distance", xlabel = "Generation", ylabel = "Max Pairwise Distance", label = "")
#     period_plot = plot(title = "Period", xlabel = "Generation", ylabel = "Seconds", label = "")
#     amplitude_plot = plot(title = "Relative Amplitude", xlabel = "Generation", ylabel = "% of Initial AP2", label = "")
#     fitness_plot = plot(title = "Mean Fitness", xlabel = "Generation", ylabel = "Fitness", label = "")
#     cluster_plot = plot(title = "Optimal Clusters N", xlabel = "Generation", ylabel = "N Optimal Clusters", label = "")
#     clusterdist_plot = plot(title = "Max & Mean Cluster Distances", xlabel = "Generation", ylabel = "Max Cluster Distances", label = "")

#     # Calculate metrics for each GA run
#     for (i, df) in enumerate(dfarray)
#         metrics = cum ? calculate_cumulative_ga_metrics(df) : calculate_ga_metrics(df)

#         generations, num_points, max_pairwise_distances, period_maxes, period_mins, period_means, amplitude_maxes, amplitude_mins, amplitude_means, mean_fitness_vals, number_of_clusters, max_cluster_distances, mean_cluster_distances = metrics

#         # Add the metrics to the subplots
#         plot!(num_points_plot, generations, num_points, color = i, legend=false)
#         plot!(pairwise_plot, generations, max_pairwise_distances, label = "", color = i, legend=false)
#         plot!(period_plot, generations, period_maxes, label = "$(hyperparam_vals[i])", color = i, legend=false)
#         plot!(amplitude_plot, generations, amplitude_maxes, label = "", color = i, legend=false)
#         plot!(fitness_plot, generations, mean_fitness_vals, label = "", color = i, legend=false)
#         plot!(cluster_plot, generations, number_of_clusters, label = "", color = i, legend=false)
#         plot!(clusterdist_plot, generations, max_cluster_distances, label = "", color = i, legend=false)

#     end

#     # Combine the subplots into a single figure
#     suptitle = cum ? "Cumulative " : "Non-cumulative"

#     plot(pairwise_plot, period_plot, amplitude_plot, fitness_plot, cluster_plot, clusterdist_plot, suptitle=suptitle, layout=(3, 2), legend=:outertopright, legendtitle=hyperparam_name, size = (1200, 800), lw=4, dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)
# end






import CairoMakie as cm
cm.activate!()
"""
    plot_comparative_ga_metrics(dfarray::Vector{DataFrame}, hyperparam_vals::Vector{<:Real}, hyperparam_name; cum = true)

Plot various metrics between GA runs with different hyperparameters (ex. popsize, group_size).
"""
function plot_comparative_ga_metrics(dfarray::Vector, hyperparam_vals::Vector, hyperparam_name; cum = true)
    fig = cm.Figure(resolution = (1400, 1200); fontsize=20)

    xticks = collect(1:6)

    num_ax = cm.Axis(fig[1, 1]; title= "Number of Points", xlabel = "Generation", ylabel = "Number of Points", xticks=xticks)
    pairwise_ax = cm.Axis(fig[1, 2]; title="Max Pairwise Distance", xlabel = "Generation", ylabel = "Max Pairwise Distance",  xticks=xticks)
    period_ax = cm.Axis(fig[2, 1]; title="Max Period", xlabel = "Generation", ylabel = "Seconds",  xticks=xticks)
    amplitude_ax = cm.Axis(fig[2, 2]; title="Max Relative Amplitude", xlabel = "Generation", ylabel = "% of Initial AP2",  xticks=xticks)
    fitness_ax = cm.Axis(fig[3, 1]; title="Mean Fitness", xlabel = "Generation", ylabel = "Fitness",  xticks=xticks)
    cluster_ax = cm.Axis(fig[3, 2]; title="Optimal Number of Clusters", xlabel = "Generation", ylabel = "N Optimal Clusters",  xticks=xticks)
    clusterdist_ax = cm.Axis(fig[4, 1]; title="Max Inter-Cluster Distance", xlabel = "Generation", ylabel = "Max Cluster Distances",  xticks=xticks)

    cm.linkxaxes!(num_ax, pairwise_ax, period_ax, amplitude_ax, fitness_ax, cluster_ax, clusterdist_ax)

    colorpal = palette([:blue, :green, :orange, :red], length(hyperparam_vals))

    # Calculate metrics for each GA run
    for (i, df) in enumerate(dfarray)
        metrics = cum ? calculate_cumulative_ga_metrics(df) : calculate_ga_metrics(df)

        generations, num_points, max_pairwise_distances, period_maxes, period_mins, period_means, amplitude_maxes, amplitude_mins, amplitude_means, mean_fitness_vals, number_of_clusters, max_cluster_distances, mean_cluster_distances = metrics

        # Add the metrics to the subplots
        cm.lines!(num_ax, generations, num_points, label = "$(hyperparam_vals[i])", color=colorpal[i])
        cm.scatter!(num_ax, generations, num_points, color=colorpal[i])
        cm.lines!(pairwise_ax, generations, max_pairwise_distances, color=colorpal[i])
        cm.scatter!(pairwise_ax, generations, max_pairwise_distances, color=colorpal[i])
        cm.lines!(period_ax, generations, period_maxes, color=colorpal[i])
        cm.scatter!(period_ax, generations, period_maxes, color=colorpal[i])
        cm.lines!(amplitude_ax, generations, amplitude_maxes, color=colorpal[i])
        cm.scatter!(amplitude_ax, generations, amplitude_maxes, color=colorpal[i])
        cm.lines!(fitness_ax, generations, mean_fitness_vals, color=colorpal[i])
        cm.scatter!(fitness_ax, generations, mean_fitness_vals, color=colorpal[i])
        cm.lines!(cluster_ax, generations, number_of_clusters, color=colorpal[i])
        cm.scatter!(cluster_ax, generations, number_of_clusters, color=colorpal[i])
        cm.lines!(clusterdist_ax, generations, max_cluster_distances, color=colorpal[i])
        cm.scatter!(clusterdist_ax, generations, max_cluster_distances, color=colorpal[i])

    end

    # Add legend
    cm.Legend(fig[5,1:end], num_ax, hyperparam_name,  orientation = :horizontal, tellwidth = false, tellheight = true)

    fig
end

fig = plot_comparative_ga_metrics(dfs, popsizes, "Population Size")

cm.save("popsize_comparison.svg", fig)















ogprobjac = make_ODE_problem()
plot_everything(ga_df, ogprobjac; jump =15 , pathdir = "test4fixedGA_10000_5groups_0.0newinds")

rowprob = get_row_prob(ga_df[end,:], ogprobjac)
rowsol = solve_odeprob(rowprob, [6, 9, 10, 11, 12, 15, 16])
Amem = map(sum, rowsol.u)[1:end]
plot(Amem, label = "")
import OscTools: findextrema
max_idxs, max_vals, min_idxs, min_vals = findextrema(Amem; min_height = 0.5)
scatter!(max_idxs, max_vals, label = "", color = :red, markersize = 5)
scatter!(min_idxs, min_vals, label = "", color = :red, markersize = 5)

plotboth(ga_df[end, :], ogprobjac)



#< testing peak filtering
sol = solve_odeprob(ogprobjac, [6, 9, 10, 11, 12, 15, 16])
Amem = map(sum, sol.u)[1:100]
using Plots
import OscTools: findextrema
plot(Amem, label = "")
max_idxs, max_vals, min_idxs, min_vals = findextrema(Amem; min_height = 2.0)
scatter!(max_idxs, max_vals, label = "", color = :red, markersize = 5)
scatter!(min_idxs, min_vals, label = "", color = :red, markersize = 5)












CSV.write("test4fixedGA.csv", ga_df)
ga_df = test4fixedGA(50000, Symbol[], Float64[]; num_tournament_groups=10, selection_method=unique_tournament_bitarray, mutationRate=0.95)

testnestarray = [[1,2,3], [4,5,6], [7,8,9]]
copynestarray = copy(testnestarray)
copynestarray[1] = [10,11,12]


constraints = AllConstraints()
boxconstraints = OscTools.BoxConstraints([constraint.min for constraint in constraints if !constraint.isfixed], [constraint.max for constraint in constraints if !constraint.isfixed])

propertynames(boxconstraints.bounds)

for minval in boxconstraints.bounds.bx[1:2:end]
    println(minval)
end

testpop = generate_population(constraints, 100000)

import OscTools: BoxConstraints

"""
    generate_new_individuals!(offspring::Vector{Vector{Float64}}, constraints::CT) where CT <: BoxConstraints

Generates `n_newInds` individuals to fill out the `offspring` array through log-uniform sampling.
"""
function generate_new_individuals!(new_inds, constraints::CT) where CT <: BoxConstraints

    rand_vals = Vector{Float64}(undef, length(new_inds))
    
    # Populate the array
    i = 1
    for minidx in 1:2:length(constraints.bounds.bx)
        min_val, max_val = log10(constraints.bounds.bx[minidx]), log10(constraints.bounds.bx[minidx+1])
        rand_vals .= exp10.(rand(Uniform(min_val, max_val), length(new_inds)))
        
        for j in eachindex(new_inds)
            new_inds[j][i] = rand_vals[j]
        end
        i += 1
    end
    return new_inds
end

@benchmark generate_new_individuals!($testpop, $boxconstraints)





testpop = generate_population(constraints, 10)
new_inds = @view testpop[6:end]

generate_new_individuals!(new_inds, boxconstraints)
generate_new_individuals!(testpop, boxconstraints)

for i in eachindex(testpop)
    testpop[i][1] = rand()
end

testpop

testvec = zeros(10)
for i in 1:10
    testvec[i] = rand()
end
testvec