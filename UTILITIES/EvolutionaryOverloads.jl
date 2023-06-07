
#! OVERRIDES FOR Evolutionary.jl ##
"""Custom GA state type that captures additional data from the objective function in the extradata field\n
    - `T` is the type of the fitness value\n
    - `IT` is the type of the individual\n
    - `TT` is the type of the additional data from the objective function\n"""
mutable struct CustomGAState <: Evolutionary.AbstractOptimizerState  
    const N::Int  
    eliteSize::Int  
    fitness::Float64  
    fitpop::Vector{Float64}  
    extradata::Vector{Vector{Float64}}
    fittest::Vector{Float64}  
end  
Evolutionary.value(s::CustomGAState) = s.fitness #return the fitness of the fittest individual
Evolutionary.minimizer(s::CustomGAState) = s.fittest #return the fittest individual

"""Trace override function"""
function Evolutionary.trace!(record::Dict{String,Any}, objfun, state, population::Vector{Vector{Float64}}, method::GA, options) 
    oscillatory_population_idxs = findall(fit -> fit < -0.1, state.fitpop) #find the indices of the oscillatory individuals
    record["staterecord"] = [(ind=population[i], fit=state.fitpop[i], per=state.extradata[i][1], amp=state.extradata[i][2]) for i in oscillatory_population_idxs]
    record["num_oscillatory"] = length(oscillatory_population_idxs)
end

function Evolutionary.trace!(record::Dict{String,Any}, objfun, state, population, method::GA, options)
    @info "WRONG trace override function called"
    record["staterecord"] = [(ind=population[i], fit=state.fitpop[i]) for i in eachindex(population)]
end

"""Show override function to prevent printing large arrays"""
function Evolutionary.show(io::IO, t::Evolutionary.OptimizationTraceRecord)
    print(io, lpad("$(t.iteration)",6))
    print(io, "   ")
    print(io, lpad("$(t.value)",14))
    for (key, value) in t.metadata
        if !isa(value, AbstractArray)
            print(io, "\n * $key: $value")
        end
    end
    return
end


"""Modified value! function from Evolutionary.jl to allow for multiple outputs from the objective function to be stored"""
function Evolutionary.NLSolversBase.value!(obj::EvolutionaryObjective{TC,TF,TX,Val{:thread}},
                                F::AbstractVector, E::AbstractVector, xs::AbstractVector{TX}) where {TC,TF<:AbstractVector,TX}
    n = length(xs)
    # @info "Evaluating $(n) individuals in parallel"
    Threads.@threads for i in 1:n
        F[i], E[i]... = Evolutionary.value(obj, xs[i])  # get the vector
    end
    # F, E
end

function Evolutionary.NLSolversBase.value!(obj::EvolutionaryObjective{TC,TF,TX,Val{:serial}},
                                F::AbstractVector, E::AbstractVector, xs::AbstractVector{TX}) where {TC,TF<:AbstractVector,TX}
    n = length(xs)
    # @info "Evaluating $(n) individuals in serial"
    for i in 1:n
        F[i], E[i]... = Evolutionary.value(obj, xs[i])  # get the vector
    end
    # F, E
end

"""Initialization of my custom GA algorithm state that captures additional data from the objective function\n
    - `method` is the GA method\n
    - `options` is the options dictionary\n
    - `objfun` is the objective function\n
    - `population` is the initial population, specifically a Vector for dispatch\n
    - `extradata` is the additional data from the objective function\n
    - `fittest` is the fittest individual\n"""
function Evolutionary.initial_state(method::GA, options, objfun, population) 
    # @show T = typeof(value(objfun))

    N = length(first(population))
    fitness = zeros(Float64, method.populationSize)
    extradata = Vector{Vector{Float64}}(undef, method.populationSize)


    # setup state values
    eliteSize = isa(method.ɛ, Int) ? method.ɛ : round(Int, method.ɛ * method.populationSize)

    # Evaluate population fitness, extradata (period and amplitude)
    Threads.@threads for i in eachindex(population)
        fitness[i], extradata[i]... = Evolutionary.value(objfun, population[i])
    end
    minfit, fitidx = findmin(fitness)

    # setup initial state
    return CustomGAState(N, eliteSize, minfit, fitness, extradata, copy(population[fitidx]))
end

"""Update state function that captures additional data from the objective function"""
function Evolutionary.update_state!(objfun, constraints, state::CustomGAState, parents::AbstractVector{IT}, method::GA, options, itr) where {IT}
    populationSize = method.populationSize
    evaltype = options.parallelization
    rng = options.rng
    offspring = similar(parents)

    # select offspring
    selected = method.selection(state.fitpop, populationSize, rng=rng)

    # perform mating
    offspringSize = populationSize - state.eliteSize
    Evolutionary.recombine!(offspring, parents, selected, method, offspringSize, rng=rng)

    # Elitism (copy population individuals before they pass to the offspring & get mutated)
    fitidxs = sortperm(state.fitpop)
    for i in 1:state.eliteSize
        subs = offspringSize+i
        offspring[subs] = copy(parents[fitidxs[i]])
    end
    # perform mutation
    Evolutionary.mutate!(offspring, method, constraints, rng=rng)

    # calculate fitness and extradata of the population
    # @info "Evaluating offspring in update_state!"
    Evolutionary.evaluate!(objfun, state.fitpop, state.extradata, offspring, constraints)
    # @info "Finished evaluating offspring in update_state!"

    # select the best individual
    minfit, fitidx = findmin(state.fitpop)
    state.fittest = offspring[fitidx]
    state.fitness = state.fitpop[fitidx]
    
    # replace population
    parents .= offspring

    return false
end

"""Modified evaluate! function from Evolutionary.jl to allow for multiple outputs from the objective function to be stored"""
function Evolutionary.evaluate!(objfun, fitness, extradata, population, constraints)
    # calculate fitness of the population
    Evolutionary.value!(objfun, fitness, extradata, population)
    # apply penalty to fitness
    Evolutionary.penalty!(fitness, constraints, population)
end