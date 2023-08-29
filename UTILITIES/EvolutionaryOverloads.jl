
import Evolutionary

#< OVERRIDES FOR Evolutionary.jl ##
"""Custom GA state type that captures additional data from the objective function in the extradata field\n
    - `T` is the type of the fitness value\n
    - `IT` is the type of the individual\n
    - `TT` is the type of the additional data from the objective function\n"""
mutable struct CustomGAState <: Evolutionary.AbstractOptimizerState  
    N::Int  #* number of elements in an individual
    eliteSize::Int  #* number of individuals that are copied to the next generation
    fittestValue::Float64  #* fitness of the fittest individual
    fitvals::Vector{Float64}  #* fitness values of the population
    fittestInd::Vector{Float64}  #* fittest individual
    periods::Vector{Float64} #* periods of the individuals
    amplitudes::Vector{Float64} #* amplitudes of the individuals
end  
Evolutionary.value(s::CustomGAState) = s.fittestValue #return the fitness of the fittest individual
Evolutionary.minimizer(s::CustomGAState) = s.fittestInd #return the fittest individual


"""Trace override function"""
function Evolutionary.trace!(record::Dict{String,Any}, objfun, state, population::Vector{Vector{Float64}}, method::GA, options) 
    oscillatory_population_idxs = findall(fit -> fit > 0.0, state.fitvals) #find the indices of the oscillatory individuals

    record["population"] = deepcopy(population[oscillatory_population_idxs])
    record["fitvals"] = copy(state.fitvals[oscillatory_population_idxs])
    record["periods"] = copy(state.periods[oscillatory_population_idxs])
    record["amplitudes"] = copy(state.amplitudes[oscillatory_population_idxs])

    # record["population"] = deepcopy(population)
    # record["fitvals"] = deepcopy(state.fitvals)
    # record["periods"] = deepcopy(state.periods)
    # record["amplitudes"] = deepcopy(state.amplitudes)

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
    print(io, "\n * num_oscillatory: $(length(t.metadata["fitvals"]))")
    return
end


"""
    EvolutionaryObjective(f, x[, F])

Constructor for an objective function object around the function `f` with initial parameter `x`, and objective value `F`.
"""
function Evolutionary.EvolutionaryObjective(f::TC, x::AbstractArray, F::Vector{Float64};
                               eval::Symbol = :serial) where {TC}
    # @info "Using custom EvolutionaryObjective constructor"
    defval = Evolutionary.default_values(x)

    #* convert function into the in-place one
    TF = typeof(F)
    fn = (Fv,xv) -> (Fv .= f(xv))
    TN = typeof(fn)
    EvolutionaryObjective{TN,TF,typeof(x),Val{eval}}(fn, F, defval, 0)
end

"""Override of the multiobjective check"""
Evolutionary.ismultiobjective(obj) = false

"""Modified value! function from Evolutionary.jl to allow for multiple outputs from the objective function to be stored"""
function Evolutionary.value!(obj::EvolutionaryObjective{TC,TF,TX,Val{:thread}},
                                F::AbstractVector, xs::AbstractVector{TX},  P::AbstractVector, A::AbstractVector) where {TC,TF<:AbstractVector,TX}
    n = length(xs)
    Threads.@threads for i in 1:n
        F[i], P[i], A[i] = Evolutionary.value(obj, xs[i])  #* evaluate the fitness, period, and amplitude for each individual
    end
    F, P, A
end

function Evolutionary.value!(obj::EvolutionaryObjective{TC,Vector{Float64},Vector{Float64},Val{:thread}},
                                F::Vector{Float64}, xs::Vector{Vector{Float64}},  P::Vector{Float64}, A::Vector{Float64}) where {TC}
    # @info "Using custom value! function"
    n = length(xs)
    Threads.@threads for i in 1:n
        F[i], P[i], A[i] = Evolutionary.value(obj, xs[i])  #* evaluate the fitness, period, and amplitude for each individual
    end
    F, P, A
end

"""Same value! function but with serial eval"""
function Evolutionary.value!(obj::EvolutionaryObjective{TC,TF,TX,Val{:serial}},
                                F::AbstractVector, xs::AbstractVector{TX}, P::AbstractVector, A::AbstractVector) where {TC,TF<:AbstractVector,TX}
    n = length(xs)
    for i in 1:n
        F[i], P[i], A[i] = Evolutionary.value(obj, xs[i])  #* evaluate the fitness, period, and amplitude for each individual
        # println("Ind: $(xs[i]) fit: $(F[i]) per: $(P[i]) amp: $(A[i])")
    end
    F, P, A
end

"""Initialization of my custom GA algorithm state that captures additional data from the objective function\n
    - `method` is the GA method\n
    - `options` is the options dictionary\n
    - `objfun` is the objective function\n
    - `population` is the initial population, specifically a Vector for dispatch\n
    - `extradata` is the additional data from the objective function\n
    - `fittest` is the fittest individual\n"""
function Evolutionary.initial_state(method::GA, options, objfun, population) #TODO something wrong with the fitness assignment in the first gen

    N = length(first(population))
    fitvals = zeros(Float64, method.populationSize)
    
    periods = zeros(Float64, method.populationSize)
    amplitudes = zeros(Float64, method.populationSize)
    # @info "Initializing GA state"

    #* setup state values
    eliteSize = isa(method.ɛ, Int) ? method.ɛ : round(Int, method.ɛ * method.populationSize)

    #* Evaluate population fitness, period and amplitude
    Evolutionary.value!(objfun, fitvals, population, periods, amplitudes)

    minfit, fitidx = findmin(fitvals)

    #* setup initial state
    return CustomGAState(N, eliteSize, minfit, fitvals, copy(population[fitidx]), periods, amplitudes)
end

"""Modified evaluate! function from Evolutionary.jl to allow for multiple outputs from the objective function to be stored"""
function Evolutionary.evaluate!(objfun, fitvals, population::Vector{Vector{Float64}}, periods, amplitudes, constraints)

    #* calculate fitness of the population
    Evolutionary.value!(objfun, fitvals, population, periods, amplitudes)

    #* apply penalty to fitness
    Evolutionary.penalty!(fitvals, constraints, population)
end

"""Update state function that captures additional data from the objective function"""
function Evolutionary.update_state!(objfun, constraints, state::CustomGAState, parents::Vector{Vector{Float64}}, method::GA, options, itr)
    populationSize = method.populationSize
    evaltype = options.parallelization
    rng = options.rng
    offspring = similar(parents)

    #* select offspring
    selected = method.selection(state.fitvals, populationSize, rng=rng)

    #* perform mating
    offspringSize = populationSize - state.eliteSize
    Evolutionary.recombine!(offspring, parents, selected, method, offspringSize, rng=rng)

    #* Elitism (copy population individuals before they pass to the offspring & get mutated)
    fitidxs = sortperm(state.fitvals)
    for i in 1:state.eliteSize
        subs = offspringSize+i
        offspring[subs] = copy(parents[fitidxs[i]])
    end

    #* perform mutation
    Evolutionary.mutate!(offspring, method, constraints, rng=rng)

    #* calculate fitness and extradata of the population
    Evolutionary.evaluate!(objfun, state.fitvals, offspring, state.periods, state.amplitudes, constraints)

    #* select the best individual
    minfit, fitidx = findmin(state.fitvals)
    state.fittestInd = offspring[fitidx]
    state.fittestValue = state.fitvals[fitidx]
    
    #* replace population
    parents .= offspring

    return false
end




