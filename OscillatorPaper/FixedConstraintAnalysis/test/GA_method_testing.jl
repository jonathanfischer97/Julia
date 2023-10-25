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
    function tournamentN(fitness::AbstractVecOrMat{<:Real}, N::Int, N_selected::Int;
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




ga_df = test4fixedGA(10000, Symbol[], Float64[]; num_tournament_groups=10, selection_method=OscTools.unique_tournament_bitarray, n_newInds=0.1)
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