"""ParamRange struct for defining parameter ranges"""
struct ConstraintRange
    name::String
    min::Float64
    max::Float64
end

abstract type ConstraintType end

"""Wrapper around NamedTuple for parameter constraints"""
struct ParameterConstraints <: ConstraintType
    data::NamedTuple
end

"""Wrapper around NamedTuple for initial condition constraints"""
struct InitialConditionConstraints <: ConstraintType
    data::NamedTuple
end



#! CONSTRAINT RANGE FUNCTIONS
"""Parameter constraint function, returns a named tuple of parameter names and their ranges"""
function define_parameter_constraints(; karange = (-3.0, 1.0), kbrange = (-3.0, 3.0), kcatrange = (-3.0, 3.0), dfrange = (1.0, 5.0))
    # Define parameter constraint ranges
    ka_min, ka_max = karange  # uM^-1s^-1, log scale
    kb_min, kb_max = kbrange  # s^-1, log scale
    kcat_min, kcat_max = kcatrange # s^-1, log scale
    df_min, df_max = dfrange # for DF, log scale

    return ParameterConstraints((
        ka1 = ConstraintRange("ka1", ka_min, ka_max),
        kb1 = ConstraintRange("kb1", kb_min, kb_max),
        kcat1 = ConstraintRange("kcat1", kcat_min, kcat_max),
        ka2 = ConstraintRange("ka2",ka_min, ka_max),
        kb2 = ConstraintRange("kb2",kb_min, kb_max),
        ka3 = ConstraintRange("ka3",ka_min, ka_max),
        kb3 = ConstraintRange("kb3", kb_min, kb_max),
        ka4 = ConstraintRange("ka4", ka_min, ka_max),
        kb4 = ConstraintRange("kb4", kb_min, kb_max),
        ka7 = ConstraintRange("ka7", ka_min, ka_max),
        kb7 = ConstraintRange("kb7", kb_min, kb_max),
        kcat7 = ConstraintRange("kcat7", kcat_min, kcat_max),
        DF = ConstraintRange("DF", df_min, df_max)
    ))
end

"""Initial condition constraint function, returns a named tuple of variable names (L, K, P, A) and their ranges"""
function define_initialcondition_constraints(;lipidrange = (0.1, 10.0), kinaserange = (0.1, 10.0), phosphataserange = (0.1, 10.0), ap2range = (0.1, 10.0))
    # Define parameter constraint ranges
    lipid_min, lipid_max = lipidrange  # uM
    kinase_min, kinase_max = kinaserange  # uM
    phosphatase_min, phosphatase_max = phosphataserange # uM
    ap2_min, ap2_max = ap2range # uM

    return InitialConditionConstraints((
        L = ConstraintRange("PIP/PIP2", lipid_min, lipid_max),
        K = ConstraintRange("Kinase", kinase_min, kinase_max),
        P = ConstraintRange("Phosphatase", phosphatase_min, phosphatase_max),
        A = ConstraintRange("AP2", ap2_min, ap2_max)
    ))
end



#! GENERATE POPULATION FUNCTIONS ##

"""Function that generates population of log-uniform sampled random parameter values"""
function generate_population(param_constraints::ParameterConstraints, n::Int)
    population = [exp10.(rand(Uniform(param.min, param.max), n)) for param in param_constraints.data]
    population = transpose(hcat(population...))

    return [population[:, i] for i in 1:n]
end

"""Function that generates population of uniform sampled random initial conditions"""
function generate_population(initialcondition_constraints::InitialConditionConstraints, n::Int)
    population = [rand(Uniform(param.min, param.max), n) for param in initialcondition_constraints.data]
    population = transpose(hcat(population...))

    return [population[:, i] for i in 1:n]
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



#! MAIN FUNCTIONS ##
"""Returns the fitness function for the cost function, referencing the ODE problem and tracker with closure"""
function make_fitness_function(evalfunc!::Function, prob::ODEProblem, peramp_tracker::PeriodAmplitudes)
    function fitness_function(input::Vector{Float64})
        #? Returns a cost function method that takes in just a vector of parameters/ICs and references the ODE problem and tracker
        return evalfunc!(input, prob, peramp_tracker)
    end
    return fitness_function
end

"""Runs the genetic algorithm, returning the result, array of all evaluated point vectors with their fitness, and the tracker values"""
function run_GA(constraints::ConstraintType, prob::ODEProblem, fitness_function_factory::Function = make_fitness_function; population_size = 10000,
                    abstol=1e-12, reltol=1e-10, successive_f_tol = 5, iterations=10, parallelization = :thread)

    # Generate the initial population.
    pop = generate_population(constraints, population_size)

    # Create constraints using the min and max values from param_values.
    boxconstraints = BoxConstraints([constraintrange.min for constraintrange in constraints.data], [constraintrange.max for constraintrange in constraints.data])

    # Define options for the GA.
    opts = Evolutionary.Options(abstol=abstol, reltol=reltol, successive_f_tol = successive_f_tol, iterations=iterations, 
                        store_trace = true, show_trace=true, show_every=1, parallelization=parallelization)

    # Define the range of possible values for each parameter.
    mutation_scalar = 0.5; mutation_range = fill(mutation_scalar, length(constraints.data))

    # Define the GA method.
    mthd = GA(populationSize = population_size, selection = tournament(Int(population_size/10)),
    crossover = TPX, crossoverRate = 0.5,
    mutation  = BGA(mutation_range, 2), mutationRate = 0.7)

    # Create a tracker to store the period and amplitude of each individual.
    peramp_tracker = PeriodAmplitudes()

    # Define the fitness function by calling the closure (holds prob, tracker, etc.).
    fitness_function = fitness_function_factory(eval_param_fitness, prob, peramp_tracker)

    # Run the optimization.
    result = Evolutionary.optimize(fitness_function, boxconstraints, mthd, pop, opts)

    # Extract all evaluated point vectors with their fitness values.
    allpoints = [gen.metadata["populationmap"] for gen in result.trace]
    allpoints = reduce(vcat, allpoints)

    return result, allpoints, peramp_tracker.peramps
end



#! MISCALAENOUS FUNCTIONS ##
"""Find the indices of the inputs in a NAME array"""
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

