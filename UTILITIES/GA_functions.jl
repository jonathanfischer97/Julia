#< CONSTRAINT RANGE STRUCTS

abstract type ConstraintType end

Base.iterate(C::ConstraintType, state=1) = state > length(fieldnames(typeof(C))) ? nothing : (getfield(C, state), state + 1)
Base.length(C::ConstraintType) = length(fieldnames(typeof(C)))
Base.eltype(::Type{<:ConstraintType}) = ConstraintRange

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
    min::Float64
    max::Float64
end

"""
    ParameterConstraints

Struct encapsulating parameter constraints. Each field represents a different parameter, holding a `ConstraintRange` object that defines the valid range for that parameter.
"""
mutable struct ParameterConstraints <: ConstraintType
    constraintranges::Vector{ConstraintRange}
end

"""
    InitialConditionConstraints

Struct encapsulating initial condition constraints. Each field represents a different initial condition, holding a `ConstraintRange` object that defines the valid range for that initial condition.
"""
mutable struct InitialConditionConstraints <: ConstraintType
    constraintranges::Vector{ConstraintRange}
    L::ConstraintRange #* PIP
    K::ConstraintRange #* PIP5K kinase
    P::ConstraintRange #* Synaptojanin phosphotase
    A::ConstraintRange #* AP2 adaptor
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
    # Define parameter constraint ranges
    ka_min, ka_max = karange  # uM^-1s^-1, log scale
    kb_min, kb_max = kbrange  # s^-1, log scale
    kcat_min, kcat_max = kcatrange # s^-1, log scale
    df_min, df_max = dfrange # for DF, log scale

    return ParameterConstraints(
        [
        :ka1 => ConstraintRange("ka1", ka_min, ka_max),
        :kb1 => ConstraintRange("kb1", kb_min, kb_max),
        :kcat1 => ConstraintRange("kcat1", kcat_min, kcat_max),
        :ka2 => ConstraintRange("ka2",ka_min, ka_max),
        :kb2 => ConstraintRange("kb2",kb_min, kb_max),
        :ka3 => ConstraintRange("ka3",ka_min, ka_max),
        :kb3 => ConstraintRange("kb3", kb_min, kb_max),
        :ka4 => ConstraintRange("ka4", ka_min, ka_max),
        :kb4 => ConstraintRange("kb4", kb_min, kb_max),
        :ka7 => ConstraintRange("ka7", ka_min, ka_max),
        :kb7 => ConstraintRange("kb7", kb_min, kb_max),
        :kcat7 => ConstraintRange("kcat7", kcat_min, kcat_max),
        :DF => ConstraintRange("DF", df_min, df_max)
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
function define_initialcondition_constraints(;lipidrange = (0.1, 10.0), kinaserange = (0.1, 10.0), phosphataserange = (0.1, 10.0), ap2range = (0.1, 10.0))
    # Define parameter constraint ranges
    lipid_min, lipid_max = lipidrange  # uM
    kinase_min, kinase_max = kinaserange  # uM
    phosphatase_min, phosphatase_max = phosphataserange # uM
    ap2_min, ap2_max = ap2range # uM

    return InitialConditionConstraints(
        [
        :L => ConstraintRange("PIP+PIP2", lipid_min, lipid_max),
        :K => ConstraintRange("Kinase", kinase_min, kinase_max),
        :P => ConstraintRange("Phosphatase", phosphatase_min, phosphatase_max),
        :A => ConstraintRange("AP2", ap2_min, ap2_max)
        ]
    )
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

    function GAProblem(constraints::T, ode_problem::ODEProblem) where T <: ConstraintType
        evalfunc = T <: ParameterConstraints ? eval_param_fitness : eval_ic_fitness
        # fitness_function = fitness_function_maker(evalfunc, ode_problem)
        new{T}(constraints, ode_problem, evalfunc)
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
    population = [exp10.(rand(Uniform(param.min, param.max), n)) for param in constraint]
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
    population = [rand(Uniform(param.min, param.max), n) for param in constraint]
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


#< RUN GENETIC ALGORITHM OPTIMIZATION ##
"""
Runs the genetic algorithm, returning the `result`, and the `record` named tuple
"""
function run_GA(ga_problem::GAProblem, fitnessfunction_factory::Function=make_fitness_function; population_size = 10000, abstol=1e-12, reltol=1e-10, successive_f_tol = 5, iterations=10, parallelization = :thread)
    # Generate the initial population.
    pop = generate_population(ga_problem.constraints, population_size)

    # Create constraints using the min and max values from param_values.
    boxconstraints = BoxConstraints([constraint.min for constraint in ga_problem.constraints], [constraint.max for constraint in ga_problem.constraints])

    # Define options for the GA.
    opts = Evolutionary.Options(abstol=abstol, reltol=reltol, successive_f_tol = successive_f_tol, iterations=iterations, 
                        store_trace = true, show_trace=true, show_every=1, parallelization=parallelization)

    # Define the range of possible values for each parameter.
    mutation_scalar = 0.5; mutation_range = fill(mutation_scalar, length(ga_problem.constraints))

    # Define the GA method.
    mthd = GA(populationSize = population_size, selection = tournament(Int(population_size/10)),
    crossover = TPX, crossoverRate = 0.5,
    mutation  = BGA(mutation_range, 2), mutationRate = 0.7)

    # Make fitness function
    # @code_warntype fitnessfunction_factory(ga_problem.eval_function, ga_problem.ode_problem)
    fitness_function = fitnessfunction_factory(ga_problem.eval_function, ga_problem.ode_problem)

    # Run the optimization.
    result = Evolutionary.optimize(fitness_function, boxconstraints, mthd, pop, opts)

    # Get the individual, fitness, and extradata of the population
    record = reduce(vcat,[gen.metadata["staterecord"] for gen in result.trace])
    
    return result, record
end
#> END







#< MISCALAENOUS FUNCTIONS ##
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
#> END
