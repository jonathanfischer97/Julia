begin 
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

ga_df = test4fixedGA(10000, Symbol[], Float64[]; num_tournament_groups=4)

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
ga_df = test4fixedGA(50000, Symbol[], Float64[]; num_tournament_groups=10, selection_method=OscTools.tournament, n_newInds=0.95)

allconstraints = AllConstraints()
testpop = generate_population(allconstraints, 1000)




diversities = Float64[]
coverages = Float64[]
for n_ind_fraction in 0.1:0.1:0.9
    ga_df = test4fixedGA(10000, Symbol[], Float64[]; num_tournament_groups=5, selection_method=OscTools.tournament, n_newInds=n_ind_fraction)
    println("Diversity metric: ", diversity_metric(ga_df))
    println("Coverage metric: ", coverage_metric(ga_df))
    push!(diversities, diversity_metric(ga_df))
    push!(coverages, coverage_metric(ga_df))
end

#Plot diversity metric vs. fraction of new individuals
using Plots
plot(0.1:0.1:0.9, diversities, label = "", xlabel = "Fraction of new individuals", ylabel = "Diversity metric", legend = :bottomright, size = (1200, 800), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)

savefig("diversity_metric_vs_fraction_new_individuals.png")




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


"""
    calculate_generation_metrics(gen_df::DataFrame)

Calculate various metrics for a single generation of a Genetic Algorithm.
"""
function calculate_generation_metrics(gen_df)
    # Max Pairwise Distance Diversity of Parameters
    max_pairwise_distance = getmax_pairwise_diversity(gen_df, [:gen, :fit, :per, :amp, :relamp])

    # Coverage
    # coverage = get_coverage(gen_df, constraints)
    # push!(coverages, coverage)

    # Get the range of periods and amplitudes
    period_range = maximum(gen_df.per) - minimum(gen_df.per)
    amplitude_range = maximum(gen_df.amp) - minimum(gen_df.amp)

    # Get the mean periods and amplitudes
    mean_period = mean(gen_df.per)
    mean_amplitude = mean(gen_df.amp)

    # Get the mean fitness value
    mean_fitness_val = mean(gen_df.fit)

    # Get the number of clusters
    cluster_number = get_optimal_clusters(gen_df, 10, [:gen, :fit, :per, :amp, :relamp])

    # Get the max and mean cluster distances
    result = get_kmeans(gen_df, cluster_number, [:gen, :fit, :per, :amp, :relamp])
    cluster_distance_matrix = get_cluster_distances(result)
    max_cluster_distance = maximum(cluster_distance_matrix)
    mean_cluster_distance = mean(cluster_distance_matrix)

    return max_pairwise_distance, period_range, amplitude_range, mean_period, mean_amplitude, mean_fitness_val, cluster_number, max_cluster_distance, mean_cluster_distance
end




"""
    calculate_ga_metrics(df::DataFrame)

Calculate various metrics to assess the performance and diversity of a Genetic Algorithm.
"""
function calculate_ga_metrics(df::DataFrame)
    # Initialize arrays to store metrics
    generations = unique(df.gen)
    max_pairwise_distances = Vector{Float64}(undef, length(generations))
    # coverages = []
    period_ranges = similar(max_pairwise_distances)
    mean_periods = similar(max_pairwise_distances)
    amplitude_ranges = similar(max_pairwise_distances)
    mean_amplitudes = similar(max_pairwise_distances)
    mean_fitness_vals = similar(max_pairwise_distances)
    number_of_clusters = Vector{Int}(undef, length(generations))
    max_cluster_distances = similar(max_pairwise_distances)
    mean_cluster_distances = similar(max_pairwise_distances)


    # Calculate metrics for each generation
    for (i,gen) in enumerate(generations)
        gen_df = df[df.gen .== gen, :]

        max_pairwise_distance, period_range, amplitude_range, mean_period, mean_amplitude, mean_fitness_val, cluster_number, max_cluster_distance, mean_cluster_distance = calculate_generation_metrics(gen_df)

        max_pairwise_distances[i] = max_pairwise_distance
        period_ranges[i] = period_range
        amplitude_ranges[i] = amplitude_range
        mean_periods[i] = mean_period
        mean_amplitudes[i] = mean_amplitude
        mean_fitness_vals[i] = mean_fitness_val
        number_of_clusters[i] = cluster_number
        max_cluster_distances[i] = max_cluster_distance
        mean_cluster_distances[i] = mean_cluster_distance
    end
    return (generations = generations, max_pairwise_distances = max_pairwise_distances, period_ranges = period_ranges, mean_periods = mean_periods, amplitude_ranges = amplitude_ranges, mean_amplitudes = mean_amplitudes, mean_fitness_vals = mean_fitness_vals, number_of_clusters = number_of_clusters, max_cluster_distances = max_cluster_distances, mean_cluster_distances = mean_cluster_distances)
end

"""
    calculate_cumulative_ga_metrics(df::DataFrame)

Calculate various cumulative metrics to assess the performance and diversity of a Genetic Algorithm.
"""
function calculate_cumulative_ga_metrics(df::DataFrame)
    # Initialize arrays to store metrics
    generations = unique(df.gen)
    cum_maxpairwise_distances = Vector{Float64}(undef, length(generations))
    # coverages = []
    cum_period_ranges = similar(cum_maxpairwise_distances)
    cum_mean_periods = similar(cum_maxpairwise_distances)
    cum_amplitude_ranges = similar(cum_maxpairwise_distances)
    cum_mean_amplitudes = similar(cum_maxpairwise_distances)
    cum_mean_fitness_vals = similar(cum_maxpairwise_distances)
    cum_number_of_cluster = Vector{Int}(undef, length(generations))
    cum_max_cluster_distances = similar(cum_maxpairwise_distances)
    cum_mean_cluster_distances = similar(cum_maxpairwise_distances)

    # Calculate cumulative metrics for each generation
    for (i, gen) in enumerate(generations)
        gen_df = df[df.gen .<= gen, :]

        cum_max_pairwise_distance, cum_period_range, cum_amplitude_range, cum_mean_period, cum_mean_amplitude, cum_mean_fitness_val, cum_cluster_number, cum_max_cluster_distance, cum_mean_cluster_distance = calculate_generation_metrics(gen_df)

        cum_maxpairwise_distances[i] = cum_max_pairwise_distance
        cum_period_ranges[i] = cum_period_range
        cum_amplitude_ranges[i] = cum_amplitude_range
        cum_mean_periods[i] = cum_mean_period
        cum_mean_amplitudes[i] = cum_mean_amplitude
        cum_mean_fitness_vals[i] = cum_mean_fitness_val
        cum_number_of_cluster[i] = cum_cluster_number
        cum_max_cluster_distances[i] = cum_max_cluster_distance
        cum_mean_cluster_distances[i] = cum_mean_cluster_distance
    end
    return (generations = generations, cum_maxpairwise_distances = cum_maxpairwise_distances, cum_period_ranges = cum_period_ranges, cum_mean_periods = cum_mean_periods, cum_amplitude_ranges = cum_amplitude_ranges, cum_mean_amplitudes = cum_mean_amplitudes, cum_mean_fitness_vals = cum_mean_fitness_vals, cum_number_of_cluster = cum_number_of_cluster, cum_max_cluster_distances = cum_max_cluster_distances, cum_mean_cluster_distances = cum_mean_cluster_distances)
end


using Plots
"""
    plot_ga_metrics(df::DataFrame)

Plot various metrics to assess the performance and diversity of a Genetic Algorithm.
"""
function plot_ga_metrics(df::DataFrame; cum = true)
    metrics = cum ? calculate_cumulative_ga_metrics(df) : calculate_ga_metrics(df)
    generations, max_pairwise_distances, period_ranges, mean_periods, amplitude_ranges, mean_amplitudes, mean_fitness_vals, number_of_clusters, max_cluster_distances, mean_cluster_distances = metrics

    # Create the plots
    p1 = plot(generations, max_pairwise_distances, label="", xlabel="Generation", ylabel="Max Pairwise Distance", title="Parameter Diversity Over Generations", color=:blue)
    # p2 = plot(generations, coverages, label="", xlabel="Generation", ylabel="Coverage", title="Coverage Over Generations", color=:red)
    p2 = plot(generations, period_ranges, label="Range", xlabel="Generation", ylabel="Range", title="Period Range and Mean Over Generations", color=:green)
    plot!(p2, generations, mean_periods, label="Mean", color=:green, linestyle=:dash)
    p3 = plot(generations, amplitude_ranges, label="Range", xlabel="Generation", ylabel="Range", title="Amplitude Range and Mean Over Generations", color=:orange)
    plot!(p3, generations, mean_amplitudes, label="Mean", color=:orange, linestyle=:dash)
    p4 = plot(generations, mean_fitness_vals, label="", xlabel="Generation", ylabel="Fitness", title="Mean Fitness Over Generations", color=:purple)
    p5 = plot(generations, number_of_clusters, label="", xlabel="Generation", ylabel="Optimal Number of Kmeans Clusters", title="Number of Clusters Over Generations", color=:brown)
    p6 = plot(generations, max_cluster_distances, label="Max", xlabel="Generation", ylabel="Max Cluster Distances", title="Max & Mean Cluster Distances Over Generations", color=:crimson)
    plot!(p6, generations, mean_cluster_distances, label="Mean", color=:crimson, linestyle=:dash)

    # Combine the plots into a single figure
    plot(p1, p2, p3, p4, p5, p6, layout=(3, 2), legend=:bottomright, size = (1200, 800), lw=3, dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)
end


# Generate the plot
plot_ga_metrics(ga_df)
savefig("test4fixedGA_50000_5groups_0.95newinds_metrics.png")


gen_df = ga_df[ga_df.gen .== 3, :]

"""
    plot_cumulative_ga_metrics(df::DataFrame)

Plot various cumulative metrics to assess the performance and diversity of a Genetic Algorithm.
"""
function plot_cumulative_ga_metrics(df::DataFrame)
    metrics = calculate_ga_metrics(df)
    cumulative_metrics = [cumsum(x) for x in metrics]
    generations, pairwise_distances, period_ranges, mean_periods, amplitude_ranges, mean_amplitudes, mean_fitness_vals, number_of_cluster, max_cluster_distances, mean_cluster_distances = cumulative_metrics

    # Create the plots
    p1 = plot(generations, pairwise_distances, label="", xlabel="Generation", ylabel="Max Pairwise Distance", title="Cumulative Parameter Diversity Over Generations", color=:blue)
    # p2 = plot(generations, coverages, label="", xlabel="Generation", ylabel="Coverage", title="Coverage Over Generations", color=:red)
    p2 = plot(generations, period_ranges, label="Range", xlabel="Generation", ylabel="Range", title="Cumulative Period Range and Mean Over Generations", color=:green)
    plot!(p2, generations, mean_periods, label="Mean", color=:green, linestyle=:dash)
    p3 = plot(generations, amplitude_ranges, label="Range", xlabel="Generation", ylabel="Range", title="Cumalative Amplitude Range and Mean Over Generations", color=:orange)
    plot!(p3, generations, mean_amplitudes, label="Mean", color=:orange, linestyle=:dash)
    p4 = plot(generations, mean_fitness_vals, label="", xlabel="Generation", ylabel="Fitness", title="Cumulative Mean Fitness Over Generations", color=:purple)
    p5 = plot(generations, number_of_cluster, label="", xlabel="Generation", ylabel="Optimal Number of Kmeans Clusters", title="Number of Clusters Over Generations", color=:brown)
    p6 = plot(generations, max_cluster_distances, label="Max", xlabel="Generation", ylabel="Max Cluster Distances", title="Max & Mean Cluster Distances Over Generations", color=:crimson)
    plot!(p6, generations, mean_cluster_distances, label="Mean", color=:crimson, linestyle=:dash)

    # Combine the plots into a single figure
    plot(p1, p2, p3, p4, p5, p6, layout=(3, 2), legend=:bottomright, size = (1200, 1200), lw=3, dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)
end


groupdf = groupby(ga_df, :gen)


groupdf[3]

cluster_number = get_optimal_clusters(ga_df, 10, [:gen, :fit, :per, :amp, :relamp])
result = get_kmeans(ga_df, cluster_number, [:gen, :fit, :per, :amp, :relamp])

get_cluster_distances(result)


cumsum([1,2,3,4,5])


ogprobjac = make_ODE_problem()
plot_everything(ga_df, ogprobjac; jump =15 , pathdir = "test4fixedGA_10000_10groups_0.5newinds")

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

boxconstraints.bounds.valx

testpop = generate_population(constraints, 1000)

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

generate_new_individuals!(testpop, boxconstraints)

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