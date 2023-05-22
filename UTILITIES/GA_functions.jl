"""Function that generates population of log-uniform sampled random parameter values"""
function generate_population(param_values::OrderedDict, n::Int)
    params = keys(param_values)
    n_params = length(params)
    population = Matrix{Float64}(undef, n_params, n)
    for i in 1:n
        for (ind, param) in enumerate(params)
            min_val = log(param_values[param]["min"])
            max_val = log(param_values[param]["max"])
            min_val < max_val ? population[ind, i] = exp(rand(Uniform(min_val, max_val))) : population[ind, i] = exp(min_val)
        end
    end
    return [collect(population[:, i]) for i in 1:n]
end


#! OVERRIDES FOR Evolutionary.jl ##
"""Trace override function"""
function Evolutionary.trace!(record::Dict{String,Any}, objfun, state, population, method::GA, options)
    record["populationmap"] = [(ind=population[i], fit=state.fitpop[i]) for i in eachindex(population)]
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



"""Returns the function factory for the cost function, referencing the ODE problem and tracker with closure"""
function make_fitness_function(evalfunc!::Function, prob::ODEProblem, peramp_tracker::PeriodAmplitudes)
    function fitness_function_factory(params::Vector{Float64})
        #? Returns a cost function method that takes in just a vector of parameters and references the ODE problem and tracker
        return evalfunc!(peramp_tracker::PeriodAmplitudes, params::Vector{Float64}, prob::ODEProblem)
    end
    return fitness_function_factory
end



"""Runs the genetic algorithm, returning the result, array of all evaluated point vectors with their fitness, and the tracker values"""
function run_GA(param_values::OrderedDict, prob::ODEProblem, evalfunc!::Function, peramp_tracker::PeriodAmplitudes; population_size = 10000,
                    abstol=1e-2, reltol=1.00, successive_f_tol = 5, iterations=10)
    #? Run the genetic algorithm
    #? param_values: Dict of parameter names and their min/max values
    #? prob: ODEProblem to solve
    #? evalfunc!: Function that evaluates the fitness of an individual, and updates the period and amplitude in the tracker
    #? Returns the result, array of all evaluated point vectors with their fitness, and the tracker values

    #? Set up the GA
    keys = keys(param_values)

    pop = generate_population(param_values, population_size)

    constraints = BoxConstraints([param_values[p]["min"] for p in keys(param_values)], [param_values[p]["max"] for p in keys(param_values)])
    
    opts = Evolutionary.Options(abstol=abstol, reltol=reltol, successive_f_tol = successive_f_tol, iterations=iterations, 
                                store_trace = true, show_trace=true, show_every=1, parallelization=:thread)

    common_range = 0.5; valrange = fill(common_range, length(param_values))
    mthd = GA(populationSize = population_size, selection = tournament(Int(population_size/10)),
            crossover = TPX, crossoverRate = 0.5,
            mutation  = BGA(valrange, 2), mutationRate = 0.7)

    fitness_function = make_fitness_function(evalfunc!, prob, peramp_tracker)

    #? Run optimization
    result = Evolutionary.optimize(fitness_function, constraints, mthd, pop, opts)

    #? Get the population map, AKA all point vectors with their fitness values in a single array
    allpoints = [gen.metadata["populationmap"] for gen in result.trace]
    allpoints = reduce(vcat, allpoints)

    return result, allpoints, peramp_tracker.peramps
end



