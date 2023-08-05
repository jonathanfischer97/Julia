
import Evolutionary
#! OVERRIDES FOR Evolutionary.jl ##
"""Custom GA state type that captures additional data from the objective function in the extradata field\n
    - `T` is the type of the fitness value\n
    - `IT` is the type of the individual\n
    - `TT` is the type of the additional data from the objective function\n"""
mutable struct CustomGAState <: Evolutionary.AbstractOptimizerState  
    N::Int  
    eliteSize::Int  
    fitness::Float64  
    fitpop::Vector{Float64}  
    fittest::Vector{Float64}  
    periods::Vector{Float64}
    amplitudes::Vector{Float64}
end  
Evolutionary.value(s::CustomGAState) = s.fitness #return the fitness of the fittest individual
Evolutionary.minimizer(s::CustomGAState) = s.fittest #return the fittest individual

"""Trace override function"""
function Evolutionary.trace!(record::Dict{String,Any}, objfun, state, population::Vector{Vector{Float64}}, method::GA, options) 
    oscillatory_population_idxs = findall(fit -> fit != 0.0, state.fitpop) #find the indices of the oscillatory individuals
    record["staterecord"] = [(ind=population[i], fit=state.fitpop[i], per=state.periods[i], amp=state.amplitudes[i]) for i in oscillatory_population_idxs]
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


"""
    EvolutionaryObjective(f, x[, F])

Constructor for an objective function object around the function `f` with initial parameter `x`, and objective value `F`.
"""
function Evolutionary.EvolutionaryObjective(f::TC, x::AbstractArray, F::Vector{Float64};
                               eval::Symbol = :serial) where {TC}
    @info "Using custom EvolutionaryObjective constructor"
    defval = Evolutionary.default_values(x)
    # convert function into the in-place one
    TF = typeof(F)
    fn = (Fv,xv) -> (Fv .= f(xv))
    TN = typeof(fn)
    EvolutionaryObjective{TN,TF,typeof(x),Val{eval}}(fn, F, defval, 0)
end

"""Modified value! function from Evolutionary.jl to allow for multiple outputs from the objective function to be stored"""
function Evolutionary.value!(obj::EvolutionaryObjective{TC,TF,TX,Val{:thread}},
                                F::AbstractVector, xs::AbstractVector{TX},  P::AbstractVector, A::AbstractVector) where {TC,TF<:AbstractVector,TX}
    n = length(xs)
    # @info "Evaluating $(n) individuals in parallel"
    Threads.@threads for i in 1:n
        F[i], P[i], A[i] = Evolutionary.value(obj, xs[i])  # get the vector
    end
    # F, E
end

function Evolutionary.value!(obj::EvolutionaryObjective{TC,TF,TX,Val{:serial}},
                                F::AbstractVector, xs::AbstractVector{TX}, P::AbstractVector, A::AbstractVector) where {TC,TF<:AbstractVector,TX}
    n = length(xs)
    # @info "Evaluating $(n) individuals in serial"
    for i in 1:n
        F[i], P[i], A[i] = Evolutionary.value(obj, xs[i])  # get the vector
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
    # extradata = Vector{Vector{Float64}}(undef, method.populationSize)
    periods = zeros(Float64, method.populationSize)
    amplitudes = zeros(Float64, method.populationSize)


    # setup state values
    eliteSize = isa(method.ɛ, Int) ? method.ɛ : round(Int, method.ɛ * method.populationSize)

    # Evaluate population fitness, extradata (period and amplitude)
    Threads.@threads for i in eachindex(population)
        fitness[i], periods[i], amplitudes[i] = Evolutionary.value(objfun, population[i])
    end
    minfit, fitidx = findmin(fitness)

    # setup initial state
    return CustomGAState(N, eliteSize, minfit, fitness, copy(population[fitidx]), periods, amplitudes,)
end

"""Modified evaluate! function from Evolutionary.jl to allow for multiple outputs from the objective function to be stored"""
function Evolutionary.evaluate!(objfun, fitness, population::Vector{Vector{Float64}}, periods, amplitudes, constraints)
    # calculate fitness of the population
    Evolutionary.value!(objfun, fitness, population, periods, amplitudes)
    # apply penalty to fitness
    Evolutionary.penalty!(fitness, constraints, population)
end

"""Update state function that captures additional data from the objective function"""
function Evolutionary.update_state!(objfun, constraints, state::CustomGAState, parents::Vector{Vector{Float64}}, method::GA, options, itr)
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
    Evolutionary.evaluate!(objfun, state.fitpop, offspring, state.periods, state.amplitudes, constraints)

    # select the best individual
    minfit, fitidx = findmin(state.fitpop)
    state.fittest = offspring[fitidx]
    state.fitness = state.fitpop[fitidx]
    
    # replace population
    parents .= offspring

    return false
end

