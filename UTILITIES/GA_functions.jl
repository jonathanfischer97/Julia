#< CONSTRAINT RANGE STRUCTS

abstract type ConstraintType end

"""
    ConstraintRange

Struct for defining parameter ranges. Each instance contains a name, and a range defined by a minimum and maximum value.

# Fields
- `name::String`: Name of the parameter.
- `min::Float64`: Minimum value for the parameter.
- `max::Float64`: Maximum value for the parameter.
"""
struct ConstraintRange
    name::String
    symbol::Symbol
    min::Float64
    max::Float64
    nominal::Float64
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
function define_parameter_constraints(; karange = (-3.0, 1.0), kbrange = (-3.0, 3.0), kcatrange = (-3.0, 3.0), dfrange = (1.0, 5.0))

    nominalvals = (;ka1 = 0.009433439939827041, kb1 = 2.3550169939427845, kcat1 = 832.7213093872278, ka2 = 12.993995997539924, kb2 = 6.150972501791291,
            ka3 = 1.3481451097940793, kb3 = 0.006201726090609513, ka4 = 0.006277294665474662, kb4 = 0.9250191811994848, ka7 = 57.36471615394549, 
            kb7 = 0.04411989797898752, kcat7 = 42.288085868394326, DF = 3631.050539219606)
    # Define parameter constraint ranges
    ka_min, ka_max = karange  # uM^-1s^-1, log scale
    kb_min, kb_max = kbrange  # s^-1, log scale
    kcat_min, kcat_max = kcatrange # s^-1, log scale
    df_min, df_max = dfrange # for DF, log scale


    return ParameterConstraints(
        [
        ConstraintRange("ka1", :ka1, ka_min, ka_max, nominalvals[1]),
        ConstraintRange("kb1", :kb1, kb_min, kb_max, nominalvals[2]),
        ConstraintRange("kcat1", :kcat1, kcat_min, kcat_max, nominalvals[3]),
        ConstraintRange("ka2", :ka2, ka_min, ka_max, nominalvals[4]),
        ConstraintRange("kb2", :kb2, kb_min, kb_max, nominalvals[5]),
        ConstraintRange("ka3", :ka3, ka_min, ka_max, nominalvals[6]),
        ConstraintRange("kb3", :kb3, kb_min, kb_max, nominalvals[7]),
        ConstraintRange("ka4", :ka4, ka_min, ka_max, nominalvals[8]),
        ConstraintRange("kb4", :kb4, kb_min, kb_max, nominalvals[9]),
        ConstraintRange("ka7", :ka7, ka_min, ka_max, nominalvals[10]),
        ConstraintRange("kb7", :kb7, kb_min, kb_max, nominalvals[11]),
        ConstraintRange("kcat7", :kcat7, kcat_min, kcat_max, nominalvals[12]),
        ConstraintRange("DF", :DF, df_min, df_max, nominalvals[13])
        ]
    )
end


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
function define_initialcondition_constraints(;lipidrange = (0.1, 10.0), kinaserange = (0.1, 5.0), phosphataserange = (0.1, 5.0), ap2range = (0.1, 10.0))

    nominalvals = (;L = 3.0, K = 0.5, P = 0.3, A = 2.0)
    # Define parameter constraint ranges
    lipid_min, lipid_max = lipidrange  # uM
    kinase_min, kinase_max = kinaserange  # uM
    phosphatase_min, phosphatase_max = phosphataserange # uM
    ap2_min, ap2_max = ap2range # uM

    return InitialConditionConstraints(
        [
        ConstraintRange("PIP+PIP2", :L, lipid_min, lipid_max, nominalvals[1]),
        ConstraintRange("Kinase", :K, kinase_min, kinase_max, nominalvals[2]),
        ConstraintRange("Phosphatase", :P, phosphatase_min, phosphatase_max, nominalvals[3]),
        ConstraintRange("AP2", :A, ap2_min, ap2_max, nominalvals[4])
        ]
    )
end
#> END

# @code_warntype define_initialcondition_constraints()


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
        # evalfunc = T === ParameterConstraints ? eval_param_fitness : eval_ic_fitness
        # fitness_function = fitness_function_maker(evalfunc, ode_problem)
        new{ParameterConstraints}(constraints, ode_problem, eval_param_fitness)
    end

    function GAProblem(constraints::InitialConditionConstraints, ode_problem::ODEProblem) 
        # evalfunc = T === ParameterConstraints ? eval_param_fitness : eval_ic_fitness
        # fitness_function = fitness_function_maker(evalfunc, ode_problem)
        new{InitialConditionConstraints}(constraints, ode_problem, eval_ic_fitness)
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
function generate_population(constraint::ParameterConstraints, n::Int)
    population = [exp10.(rand(Uniform(conrange.min, conrange.max), n)) for conrange in constraint.ranges]
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
function generate_population(constraint::InitialConditionConstraints, n::Int)
    population = [rand(Uniform(conrange.min, conrange.max), n) for conrange in constraint.ranges]
    population = transpose(hcat(population...))
    return [population[:, i] for i in 1:n]
end

"""For calculating volume"""
function generate_population(constraint::Vector{ConstraintRange}, n::Int)
    population = [rand(Uniform(conrange.min, conrange.max), n) for conrange in constraint]
    population = transpose(hcat(population...))
    return [population[:, i] for i in 1:n]
end
#> END


#< DEFAULT FITNESS FUNCTION FACTORY
"""Returns the `fitness function(input)` for the cost function, referencing the ODE problem and tracker with closure"""
function make_fitness_function(evalfunc::Function, prob::ODEProblem)
    function fitness_function(input::Vector{Float64})
        #? Returns a cost function method that takes in just a vector of parameters/ICs and references the ODE problem and tracker
        return evalfunc(input, prob)
    end
    return fitness_function
end
#> END

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
function run_GA(ga_problem::GAProblem, fitnessfunction_factory::Function=make_fitness_function; threshold=10000, population_size = 10000, abstol=1e-12, reltol=1e-10, successive_f_tol = 1, iterations=10, parallelization = :thread)
    # Generate the initial population.
    pop = generate_population(ga_problem.constraints, population_size)
    # @info "Generated initial population"

    # Create constraints using the min and max values from param_values.
    boxconstraints = BoxConstraints([constraint.min for constraint in ga_problem.constraints.ranges], [constraint.max for constraint in ga_problem.constraints.ranges])
    # @info "Created box constraints"

    # Create Progress bar and callback function
    ga_progress = Progress(threshold; desc = "GA Progress")
    callback_func = (trace) -> ga_callback(trace, ga_progress, threshold)

    # Define options for the GA.
    opts = Evolutionary.Options(abstol=abstol, reltol=reltol, successive_f_tol = successive_f_tol, iterations=iterations, 
                        store_trace = true, show_trace=true, show_every=1, parallelization=parallelization, callback=callback_func)

    # Define the range of possible values for each parameter.
    mutation_scalar = 0.5; mutation_range = fill(mutation_scalar, length(ga_problem.constraints.ranges))

    # Define the GA method.
    mthd = GA(populationSize = population_size, selection = tournament(Int(population_size/10)),
    crossover = TPX, crossoverRate = 0.5,
    mutation  = BGA(mutation_range, 2), mutationRate = 0.7)

    # Make fitness function
    # @code_warntype fitnessfunction_factory(ga_problem.eval_function, ga_problem.ode_problem)
    fitness_function = fitnessfunction_factory(ga_problem.eval_function, ga_problem.ode_problem)
    # @info "Created fitness function"

    # Run the optimization.
    # @info "Starting optimization"
    result = Evolutionary.optimize(fitness_function, boxconstraints, mthd, pop, opts)
    # @info "Finished optimization"
    # return result
    # Get the individual, fitness, and extradata of the population
    record::Vector{NamedTuple{(:ind,:fit,:per,:amp),Tuple{Vector,Float64, Float64, Float64}}} = reduce(vcat,[gen.metadata["staterecord"] for gen in result.trace])
    return record
    # return record, result
end
#> END






#< MISCELLANEOUS FUNCTIONS ##
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

function find_indices(combination::Vector{ConstraintRange}, constraints::Vector{ConstraintRange})
    p1idx = findfirst(x -> x.name == combination[1].name, constraints)
    p2idx = findfirst(x -> x.name == combination[2].name, constraints)
    return p1idx, p2idx
end
#> END
