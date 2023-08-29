#< CONSTRAINT RANGE STRUCTS

abstract type ConstraintType end

"""
    ConstraintRange

Struct for defining parameter ranges. Each instance contains a name, and a range defined by a minimum and maximum value.

# Fields
- `name::String`: Name of the parameter.
- `symbol::Symbol`: Symbol of the parameter.
- `min::Float64`: Minimum value for the parameter.
- `max::Float64`: Maximum value for the parameter.
- `nominal::Float64`: Nominal value for the parameter.
"""
struct ConstraintRange
    name::String
    min::Float64
    max::Float64
    # nominal::Union{Nothing, Float64}
end

"""
    ParameterConstraints

Struct encapsulating parameter constraints. Each field represents a different parameter, holding a `ConstraintRange` object that defines the valid range for that parameter.
"""
mutable struct ParameterConstraints <: ConstraintType
    ranges::Vector{ConstraintRange}
end

"""
    InitialConditionConstraints

Struct encapsulating initial condition constraints. Each field represents a different initial condition, holding a `ConstraintRange` object that defines the valid range for that initial condition.
"""
mutable struct InitialConditionConstraints <: ConstraintType 
    ranges::Vector{ConstraintRange} 
end

"""
    AllConstraints

Struct encapsulating all constraints. Each field represents a different parameter or initial condition, holding a `ConstraintRange` object that defines the valid range for that parameter or initial condition.
"""
mutable struct AllConstraints <: ConstraintType
    ranges::Vector{ConstraintRange}
end
#> END 



#< CONSTRAINT RANGE CONSTRUCTORS
"""
    define_parameter_constraints(; kwargs...)

Define parameter constraints. Each keyword argument represents a different parameter, where the value is a tuple defining the valid range for that parameter.

# Example
```julia
constraints = define_parameter_constraints(
    karange = (-3.0, 1.0), 
    kbrange = (-3.0, 3.0), 
    kcatrange = (-3.0, 3.0), 
    dfrange = (1.0, 5.0)
)
```
"""
function define_parameter_constraints(; karange = (1e-3, 1e1), kbrange = (1e-3, 1e3), kcatrange = (1e-3, 1e3), dfrange = (1e3, 1e5))#, nominalvals = repeat([Nothing],13))
    #* Define parameter constraint ranges
    ka_min, ka_max = karange  # uM^-1s^-1, log scale
    kb_min, kb_max = kbrange  # s^-1, log scale
    kcat_min, kcat_max = kcatrange # s^-1, log scale
    df_min, df_max = dfrange # for DF, log scale

    return ParameterConstraints(
        [
        ConstraintRange("ka1", ka_min, ka_max),
        ConstraintRange("kb1", kb_min, kb_max),
        ConstraintRange("kcat1", kcat_min, kcat_max),
        ConstraintRange("ka2", ka_min, ka_max),
        ConstraintRange("kb2", kb_min, kb_max),
        ConstraintRange("ka3", ka_min, ka_max),
        ConstraintRange("kb3", kb_min, kb_max),
        ConstraintRange("ka4", ka_min, ka_max),
        ConstraintRange("kb4", kb_min, kb_max),
        ConstraintRange("ka7", ka_min, ka_max),
        ConstraintRange("kb7", kb_min, kb_max),
        ConstraintRange("kcat7", kcat_min, kcat_max),
        ConstraintRange("DF", df_min, df_max)
        ]
    )
end

# define_parameter_constraints(prob::ODEProblem; karange = (1e-3, 1e1), kbrange = (1e-3, 1e3), kcatrange = (1e-3, 1e3), dfrange = (1e3, 1e5)) = 
#                         define_parameter_constraints(; karange=karange, kbrange=kbrange, kcatrange=kcatrange, dfrange=dfrange, nominalvals = prob.p)


"""
    define_initialcondition_constraints(; kwargs...)

Define initial condition constraints. Each keyword argument represents a different initial condition, where the value is a tuple defining the valid range for that initial condition.

# Example
```julia
constraints = define_initialcondition_constraints(
    lipidrange = (0.1, 10.0), 
    kinaserange = (0.1, 10.0), 
    phosphataserange = (0.1, 10.0), 
    ap2range = (0.1, 10.0)
)
```
"""
function define_initialcondition_constraints(;Lrange = (0.1, 10.0), Krange = (0.1, 5.0), Prange = (0.1, 5.0), Arange = (0.1, 10.0))#, nominalvals = repeat([Nothing],4))
    # Define parameter constraint ranges
    lipid_min, lipid_max = Lrange  # uM
    kinase_min, kinase_max = Krange  # uM
    phosphatase_min, phosphatase_max = Prange # uM
    ap2_min, ap2_max = Arange # uM

    return InitialConditionConstraints(
        [
        ConstraintRange("L", lipid_min, lipid_max),
        ConstraintRange("K", kinase_min, kinase_max),
        ConstraintRange("P", phosphatase_min, phosphatase_max),
        ConstraintRange("A", ap2_min, ap2_max)
        ]
    )
end

# define_initialcondition_constraints(prob::ODEProblem; Lrange = (0.1, 10.0), Krange = (0.1, 5.0), Prange = (0.1, 5.0), Arange = (0.1, 10.0)) = 
#                                     define_initialcondition_constraints(;Lrange=Lrange, Krange=Krange, Prange= Prange, Arange=Arange, nominalvals = prob.u0[1:4])



function AllConstraints(paramconstraints::ParameterConstraints, icconstraints::InitialConditionConstraints) 
    AllConstraints([paramconstraints.ranges; icconstraints.ranges])
end
#> END




#< GA PROBLEM TYPE
"""
    GAProblem{T <: ConstraintType}

Struct encapsulating a Genetic Algorithm (GA) optimization problem. It holds the constraints for the problem, the ODE problem to be solved, the fitness function, and any additional keyword arguments.

# Fields
- `constraints::T`: Constraints for the problem. Either `ParameterConstraints` or `InitialConditionConstraints`.
- `ode_problem::ODEProblem`: ODE problem to be solved.
- `fitness_function::Function`: Fitness function, automatically generated with constructor
- `options::NamedTuple`: Additional keyword arguments for `Evolutionary.Options`.
"""
struct GAProblem{T <: ConstraintType}
    constraints::T
    ode_problem::ODEProblem
    eval_function::Function

    function GAProblem(constraints::ParameterConstraints, ode_problem::ODEProblem) 
        new{ParameterConstraints}(constraints, ode_problem, eval_param_fitness)
    end

    function GAProblem(constraints::InitialConditionConstraints, ode_problem::ODEProblem) 
        new{InitialConditionConstraints}(constraints, ode_problem, eval_ic_fitness)
    end

    function GAProblem(constraints::AllConstraints, ode_problem::ODEProblem) 
        new{AllConstraints}(constraints, ode_problem, eval_all_fitness)
    end
end


function Base.show(io::IO, ::MIME"text/plain", prob::GAProblem) #TODO add labels for nominal values
    printstyled(io, "GAProblem with constraints:\n"; bold = true, underline=true, color = :green)
    printstyled(io, prob.constraints, "\n")
    printstyled(io, "\nNominal parameter values:\n"; bold = true, color = :blue)
    printstyled(io, prob.ode_problem.p, "\n")
    printstyled(io, "\nNominal initial conditions:\n"; bold = true, color = :blue)
    printstyled(io, prob.ode_problem.u0, "\n")
end

#> END

#< POPULATION GENERATION METHODS
"""
    generate_population(constraints::ParameterConstraints, n::Int)

Generate a population of `n` individuals for the given parameter constraints. Each individual is sampled from a log-uniform distribution within the valid range for each parameter.

# Example
```julia
constraints = define_parameter_constraints(karange = (-3.0, 1.0), kbrange = (-3.0, 3.0))
population = generate_population(constraints, 100)
```
"""
function generate_population(constraint::ConstraintType, n::Int)
    population = [exp10.(rand(Uniform(log10(conrange.min), log10(conrange.max)), n)) for conrange in constraint.ranges]
    population = transpose(hcat(population...))
    return [population[:, i] for i in 1:n]
end

"""
generate_population(constraints::InitialConditionConstraints, n::Int)

Generate a population of `n` individuals for the given initial condition constraints. Each individual is sampled from a uniform distribution within the valid range for each initial condition.

# Example
```julia
constraints = define_initialcondition_constraints(lipidrange = (0.1, 10.0), kinaserange = (0.1, 10.0))
population = generate_population(constraints, 100)
```
"""
# function generate_population(constraint::InitialConditionConstraints, n::Int)
#     population = [rand(Uniform(conrange.min, conrange.max), n) for conrange in constraint.ranges]
#     population = transpose(hcat(population...))
#     return [population[:, i] for i in 1:n]
# end

# function generate_population(constraint::AllConstraints, n::Int)
#     #* find index where initial conditions start, will be named either L, K, P, or A
#     icstartidx = findfirst(x -> x.name ∈ ["L", "K", "P", "A"], constraint.ranges)
#     parampopulation = [exp10.(rand(Uniform(log10(conrange.min), log10(conrange.max)), n)) for conrange in constraint.ranges[1:icstartidx-1]]
#     icpopulation = [rand(Uniform(conrange.min, conrange.max), n) for conrange in constraint.ranges[icstartidx:end]]
#     parampopulation = transpose(hcat(parampopulation...))
#     icpopulation = transpose(hcat(icpopulation...))
#     return [[parampopulation[:,i]; icpopulation[:,i]] for i in 1:n]
#     # return [population[:, i] for i in 1:n]
# end

"""For calculating volume when optimizing for NERDSS solutions"""
function generate_population(constraint::Vector{ConstraintRange}, n::Int)
    population = [rand(Uniform(conrange.min, conrange.max), n) for conrange in constraint]
    population = transpose(hcat(population...))
    return [population[:, i] for i in 1:n]
end
#> END


#< DEFAULT FITNESS FUNCTION FACTORY
"""Returns the `fitness function(input)` for the cost function, referencing the ODE problem with closure"""
function make_fitness_function(evalfunc::Function, prob::ODEProblem)
    function fitness_function(input::Vector{Float64})
        #? Returns a cost function method that takes in just a vector of parameters/ICs and references the ODE problem 
        return evalfunc(input, prob)
    end
    return fitness_function
end

"""Returns the `fitness function(input)` for the cost function, referencing the ODE problem with closure, captured with let keyword"""
function make_fitness_function_with_let(evalfunc::Function, prob::ODEProblem)
    let evalfunc = evalfunc, prob = prob
        function fitness_function(input::Vector{Float64})
            return evalfunc(input, prob)
        end
        return fitness_function
    end
end
#> END


"""### Callback function that terminates the GA if the number of oscillations exceeds the threshold, and updates the progress bar"""
function ga_callback(trace::Evolutionary.OptimizationTrace, progressbar::Progress, threshold::Int)
    #? Callback function for the GA, updating the progress bar
    num_oscillation = trace[end].metadata["num_oscillatory"]
    if num_oscillation >= threshold 
        finish!(progressbar)
        return true
    else
        next!(progressbar, step = num_oscillation)
        return false
    end
end


#< RUN GENETIC ALGORITHM OPTIMIZATION ##
"""
Runs the genetic algorithm, returning the `result`, and the `record` named tuple
"""
function run_GA(ga_problem::GAProblem, fitnessfunction_factory::Function=make_fitness_function; 
                                            threshold=10000, population_size = 5000, abstol=1e-4, reltol=1e-2, successive_f_tol = 2, iterations=5, parallelization = :thread, show_trace=true)
    blas_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)

    #* Generate the initial population.
    pop = generate_population(ga_problem.constraints, population_size)

    #* Create constraints using the min and max values from param_values.
    boxconstraints = BoxConstraints([constraint.min for constraint in ga_problem.constraints.ranges], [constraint.max for constraint in ga_problem.constraints.ranges])

    # *Create Progress bar and callback function
    # ga_progress = Progress(threshold; desc = "GA Progress")
    # callback_func = (trace) -> ga_callback(trace, ga_progress, threshold)

    #* Define options for the GA.
    opts = Evolutionary.Options(abstol=abstol, reltol=reltol, successive_f_tol = successive_f_tol, iterations=iterations, 
                        store_trace = true, show_trace=show_trace, show_every=1, parallelization=parallelization)#, callback=callback_func)

    #* Define the range of possible values for each parameter when mutated, and the mutation scalar.
    mutation_scalar = 0.5
    mutation_range = fill(mutation_scalar, length(ga_problem.constraints.ranges))

    #* Define the GA method.
    mthd = GA(populationSize = population_size, selection = tournament(Int(population_size/10), select=argmax),
    crossover = TPX, crossoverRate = 1.0, # Two-point crossover event
    mutation  = BGA(mutation_range, 2), mutationRate = 1.0)

    #* Make fitness function. Makes closure of evaluation function and ODE problem
    fitness_function = fitnessfunction_factory(ga_problem.eval_function, ga_problem.ode_problem)

    #* Run the optimization.
    result = Evolutionary.optimize(fitness_function, [0.0,0.0,0.0], boxconstraints, mthd, pop, opts)
    # return result

    #* Get the individual, fitness, period/amplitude, for each oscillatory individual evaluated
    # record::Vector{NamedTuple{(:ind,:fit,:per,:amp),Tuple{Vector{Float64},Float64, Float64, Float64}}} = reduce(vcat,[gen.metadata["staterecord"] for gen in result.trace[2:end]])
    # num_oscillatory = sum([gen.metadata["num_oscillatory"] for gen in result.trace[2:end]])

    BLAS.set_num_threads(blas_threads)
    # return result
    return GAResults(result, length(ga_problem.constraints.ranges))
    # return trace_to_df(result)
end
#> END

"Struct to hold the results of a GA optimization"
struct GAResults 
    trace::Vector{Evolutionary.OptimizationTraceRecord}
    population::Vector{Vector{Float64}}
    fitvals::Vector{Float64}
    periods::Vector{Float64}
    amplitudes::Vector{Float64}
end

function GAResults(result::Evolutionary.EvolutionaryOptimizationResults, indlength::Int = 13)
    numpoints = sum(length, (gen.metadata["fitvals"] for gen in result.trace))
    population = [Vector{Float64}(undef, indlength) for i in 1:numpoints]
    fitvals = Vector{Float64}(undef, numpoints)
    periods = Vector{Float64}(undef, numpoints)
    amplitudes = Vector{Float64}(undef, numpoints)

    startidx = 1
    for gen in result.trace
        endidx = startidx + length(gen.metadata["population"]) - 1
        # @info startidx, endidx
        # push!(population, gen.metadata["population"]...)
        population[startidx:endidx] .= gen.metadata["population"]
        # push!(fitvals, gen.metadata["fitvals"]...)
        fitvals[startidx:endidx] .= gen.metadata["fitvals"]
        # push!(periods, gen.metadata["periods"]...)
        periods[startidx:endidx] .= gen.metadata["periods"]
        # push!(amplitudes, gen.metadata["amplitudes"]...)
        amplitudes[startidx:endidx] .= gen.metadata["amplitudes"]
        startidx = endidx + 1
    end
    return GAResults(result.trace, population, fitvals, periods, amplitudes)
end

function make_df(results::GAResults)
    return DataFrame(ind = results.population, fit = results.fitvals, per = results.periods, amp = results.amplitudes)
end



# function trace_to_df(results)
#     #* make a dataframe from the trace
#     df = DataFrame(ind = [], fit = [], per = [], amp = [])
#     for gen in results.trace
#         push!(df.ind, gen.metadata["population"]...)
#         push!(df.fit, gen.metadata["fitvals"]...)
#         push!(df.per, gen.metadata["periods"]...)
#         push!(df.amp, gen.metadata["amplitudes"]...)
#     end
#     return df
# end


#< MISCELLANEOUS FUNCTIONS ##
"""Defines logspace function for sampling parameters"""
logrange(start, stop, length::Int) = exp10.(collect(range(start=log10(start), stop=log10(stop), length=length)))



"""Splits ind column into separate columns for each parameter, adds initial conditions for writing DataFrame to CSV"""
function split_dataframe!(df, prob)
    paramsymbols = [:ka1,:kb1,:kcat1,:ka2,:kb2,:ka3,:kb3,:ka4,:kb4,:ka7,:kb7,:kcat7,:DF]
    
    if length(df.ind[1]) == 13

        #* Split ind column into separate columns for each parameter
        for (i,param) in enumerate(paramsymbols)
            df[!, param] .= [x[i] for x in df.ind]
        end

        select!(df, Not(:ind))

        #* Add initial conditions
        df.L .= prob.u0[1]
        df.K .= prob.u0[2]
        df.P .= prob.u0[3]
        df.A .= prob.u0[4]
    elseif length(df.ind[1]) == 4
        concsymbols = [:L,:K,:P,:A]
        for (i,conc) in enumerate(concsymbols)
            df[!, conc] .= [x[i] for x in df.ind]
        end
        select!(df, Not(:ind))

        #* Add parameters
        for (i,param) in paramsymbols
            df[!, param] .= prob.p[i]
        end
    else 
        concsymbols = [:L,:K,:P,:A]
        allsymbols = vcat(paramsymbols, concsymbols)

        for (i,conc) in enumerate(allsymbols)
            df[!, conc] .= [x[i] for x in df.ind]
        end
        select!(df, Not(:ind))
    end
end


"""Find the indices of the inputs in a `NAME` array"""
function find_indices(combination::Vector{String}, NAMES::Vector{String})
    p1idx = findfirst(isequal(combination[1]), NAMES)
    p2idx = findfirst(isequal(combination[2]), NAMES)
    return p1idx, p2idx
end

function find_indices(combination::Vector{Symbol}, NAMES::Vector{String})
    str1 = string(combination[1])
    str2 = string(combination[2])
    p1idx = findfirst(isequal(str1), NAMES)
    p2idx = findfirst(isequal(str2), NAMES)
    return p1idx, p2idx
end

function find_indices(combination::Vector{ConstraintRange}, constraints::Vector{ConstraintRange})::Tuple{Int,Int}
    p1idx = findfirst(x -> x.name == combination[1].name, constraints)
    p2idx = findfirst(x -> x.name == combination[2].name, constraints)
    return p1idx, p2idx
end

"""Triplet version for 3FixedParamCSVMaker"""
function find_indices(param1::String, param2::String, param3::String, constraints::Vector{ConstraintRange})::Tuple{Int,Int,Int}
    p1idx = findfirst(x -> x.name == param1, constraints)
    p2idx = findfirst(x -> x.name == param2, constraints)
    p3idx = findfirst(x -> x.name == param3, constraints)

    return p1idx, p2idx, p3idx
end
#> END
