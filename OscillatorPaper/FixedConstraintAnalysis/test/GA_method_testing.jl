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
ga_df = test4fixedGA(10000, Symbol[], Float64[]; num_tournament_groups=5, selection_method=OscTools.tournament, n_newInds=0.90)

allconstraints = AllConstraints()
testpop = generate_population(allconstraints, 1000)
hcat(testpop...)
stack(testpop)

dfmat = df_to_matrix(ga_df, [:gen, :fit, :per, :amp, :relamp])
typeof(dfmat)
typeof(dfmat) <: AbstractMatrix{Float64}

function diversity_metric(population::Vector{Vector{Float64}})
    pop_matrix = stack(population)

    diversity_metric(pop_matrix)
end


function diversity_metric(population::AbstractMatrix{Float64})
    # Step 1: Log Transformation 
    log_population = log.population

    # Step 2: Normalization
    min_vals = minimum(log_population, dims=2)
    max_vals = maximum(log_population, dims=2)
    normalized_population = (log_population .- min_vals) ./ (max_vals - min_vals)

    # Step 3 & 4: Compute Average Pairwise Distances
    num_params, n = size(normalized_population)
    @info "Number of individuals: ", n
    distances = [norm(normalized_population[:, i] - normalized_population[:, j]) for i in 1:n for j in (i+1):n]
    
    return mean(distances)
end



function diversity_metric(df::DataFrame)
    dfmat = df_to_matrix(df, [:gen, :fit, :per, :amp, :relamp])
    return diversity_metric(dfmat)
end


diversity_metric(ga_df)


diversities = Float64[]
for n_ind_fraction in 0.1:0.1:0.9
    ga_df = test4fixedGA(10000, Symbol[], Float64[]; num_tournament_groups=5, selection_method=OscTools.tournament, n_newInds=n_ind_fraction)
    println("Diversity metric: ", diversity_metric(ga_df))
    push!(diversities, diversity_metric(ga_df))
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
function coverage_metric(population::AbstractMatrix{Float64})
    # Log Transformation and Normalization as before
    log_population = log.(population .+ 1e-9)
    min_vals = minimum(log_population, dims=2)
    max_vals = maximum(log_population, dims=2)
    normalized_population = (log_population .- min_vals) ./ (max_vals - min_vals)

    # Calculate the range of each parameter
    param_ranges = maximum(normalized_population, dims=2) - minimum(normalized_population, dims=2)

    # Sum the ranges to get the overall coverage metric
    return sum(param_ranges)
end







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