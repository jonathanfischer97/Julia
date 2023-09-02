#< CONSTRAINT RANGE STRUCTS

abstract type ConstraintType end

"""
    iterate(constraint::ConstraintType)

Iterate over the constraint ranges in a `ConstraintType` object.
"""
# Define the length method
Base.length(constraint::ConstraintType) = length(fieldnames(typeof(constraint)))

# Define the getindex method for index-based access
function Base.getindex(constraint::ConstraintType, idx::Int)
    field_name = fieldnames(typeof(constraint))[idx]
    return getfield(constraint, field_name)
end

# To make it iterable, define start, next and done methods
Base.iterate(constraint::ConstraintType, state=1) = state > length(constraint) ? nothing : (getfield(constraint, fieldnames(typeof(constraint))[state]), state + 1)

# Required for the `in` keyword
Base.eltype(::Type{ConstraintType}) = ConstraintRange

# Required for length
activelength(constraint::ConstraintType) = length(filter(x -> getfield(constraint, x).active, fieldnames(typeof(constraint))))

function get_fixed_indices(constraints::ConstraintType)
    inactive_indices = Int[]  # Initialize an empty array to store the indices of inactive elements
    for (idx, constraint) in enumerate(constraints)  # Loop through each element along with its index
        if !constraint.active  # Check if the element is inactive
            push!(inactive_indices, idx)  # Add the index to the array
        end
    end
    return inactive_indices  # Return the array of indices
end



"""
    ConstraintRange

Struct for defining parameter or initial condition ranges. Each instance contains a name, and a range defined by a minimum and maximum value.

# Fields
- `name::String`: Name of the parameter or initial condtion.
- `min::Float64`: Minimum allowed value.
- `max::Float64`: Maximum allowed value.
- `active::Bool`: Whether or not the constraint will be used for optimization.
- `fixed_value::Float64`: Fixed value if inactive.
"""
mutable struct ConstraintRange
    name::String
    min::Float64
    max::Float64
    active::Bool
    fixed_value::Union{Nothing, Float64}
end

function fix_value!(conrange::ConstraintRange, value::Float64)
    conrange.active = false
    conrange.fixed_value = value
end

"""
    ParameterConstraints

Struct encapsulating parameter constraints. Each field represents a different parameter, holding a `ConstraintRange` object that defines the valid range for that parameter.
"""
# mutable struct ParameterConstraints <: ConstraintType
#     ranges::Vector{ConstraintRange}
# end

@kwdef struct ParameterConstraints <: ConstraintType
    ka1::ConstraintRange
    kb1::ConstraintRange
    kcat1::ConstraintRange
    ka2::ConstraintRange
    kb2::ConstraintRange
    ka3::ConstraintRange
    kb3::ConstraintRange
    ka4::ConstraintRange
    kb4::ConstraintRange
    ka7::ConstraintRange
    kb7::ConstraintRange
    kcat7::ConstraintRange
    DF::ConstraintRange
end

"""
    InitialConditionConstraints

Struct encapsulating initial condition constraints. Each field represents a different initial condition, holding a `ConstraintRange` object that defines the valid range for that initial condition.
"""
@kwdef struct InitialConditionConstraints <: ConstraintType 
    L::ConstraintRange
    K::ConstraintRange
    P::ConstraintRange
    A::ConstraintRange
end

"""
    AllConstraints

Struct encapsulating all constraints. Each field represents a different parameter or initial condition, holding a `ConstraintRange` object that defines the valid range for that parameter or initial condition.
"""
@kwdef struct AllConstraints <: ConstraintType
    ka1::ConstraintRange
    kb1::ConstraintRange
    kcat1::ConstraintRange
    ka2::ConstraintRange
    kb2::ConstraintRange
    ka3::ConstraintRange
    kb3::ConstraintRange
    ka4::ConstraintRange
    kb4::ConstraintRange
    ka7::ConstraintRange
    kb7::ConstraintRange
    kcat7::ConstraintRange
    DF::ConstraintRange

    L::ConstraintRange
    K::ConstraintRange
    P::ConstraintRange
    A::ConstraintRange
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
        ka1 = ConstraintRange("ka1", ka_min, ka_max, true, nothing),
        kb1 = ConstraintRange("kb1", kb_min, kb_max, true, nothing),
        kcat1 = ConstraintRange("kcat1", kcat_min, kcat_max, true, nothing),
        ka2 = ConstraintRange("ka2", ka_min, ka_max, true, nothing),
        kb2 = ConstraintRange("kb2", kb_min, kb_max, true, nothing),
        ka3 = ConstraintRange("ka3", ka_min, ka_max, true, nothing),
        kb3 = ConstraintRange("kb3", kb_min, kb_max, true, nothing),
        ka4 = ConstraintRange("ka4", ka_min, ka_max, true, nothing),
        kb4 = ConstraintRange("kb4", kb_min, kb_max, true, nothing),
        ka7 = ConstraintRange("ka7", ka_min, ka_max, true, nothing),
        kb7 = ConstraintRange("kb7", kb_min, kb_max, true, nothing),
        kcat7 = ConstraintRange("kcat7", kcat_min, kcat_max, true, nothing),
        DF = ConstraintRange("DF", df_min, df_max, true, nothing)
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
function define_initialcondition_constraints(;Lrange = (0.1, 10.0), Krange = (0.1, 5.0), Prange = (0.1, 5.0), Arange = (0.1, 10.0))#, nominalvals = repeat([Nothing],4))
    # Define initial condition constraint ranges
    lipid_min, lipid_max = Lrange  # uM
    kinase_min, kinase_max = Krange  # uM
    phosphatase_min, phosphatase_max = Prange # uM
    ap2_min, ap2_max = Arange # uM

    return InitialConditionConstraints(
        L = ConstraintRange("L", lipid_min, lipid_max, true, nothing),
        K = ConstraintRange("K", kinase_min, kinase_max, true, nothing),
        P = ConstraintRange("P", phosphatase_min, phosphatase_max, true, nothing),
        A = ConstraintRange("A", ap2_min, ap2_max, true, nothing)
    )
end





function AllConstraints(paramconstraints::ParameterConstraints, icconstraints::InitialConditionConstraints) 
    return AllConstraints(
        ka1 = paramconstraints.ka1,
        kb1 = paramconstraints.kb1,
        kcat1 = paramconstraints.kcat1,
        ka2 = paramconstraints.ka2,
        kb2 = paramconstraints.kb2,
        ka3 = paramconstraints.ka3,
        kb3 = paramconstraints.kb3,
        ka4 = paramconstraints.ka4,
        kb4 = paramconstraints.kb4,
        ka7 = paramconstraints.ka7,
        kb7 = paramconstraints.kb7,
        kcat7 = paramconstraints.kcat7,
        DF = paramconstraints.DF,

        L = icconstraints.L,
        K = icconstraints.K,
        P = icconstraints.P,
        A = icconstraints.A
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

    printstyled(io, "\nFixed values:\n"; bold = true, color = :red)
    printstyled(io, [(con.name => con.fixed_value) for con in prob.constraints if !con.active], "\n")
end

#> END

#< POPULATION GENERATION METHODS
"""
    generate_population(constraints::ParameterConstraints, n::Int)

Generate a population of `n` individuals for the given generic `constraints <: ConstraintType`. Each individual is sampled from a log-uniform distribution within the valid range for each parameter or initial condition.

# Example
```julia
constraints = define_parameter_constraints(karange = (-3.0, 1.0), kbrange = (-3.0, 3.0))
population = generate_population(constraints, 100)
```
"""
# function generate_population(constraint::ConstraintType, n::Int)
#     population = [exp10.(rand(Uniform(log10(conrange.min), log10(conrange.max)), n)) for conrange in constraint]
#     population = transpose(hcat(population...))
#     return [population[:, i] for i in 1:n]
# end

function generate_population(constraint::ConstraintType, n::Int)
    num_params = length(constraint)
    
    # Preallocate the population array of arrays
    population = [Vector{Float64}(undef, num_params) for _ in 1:n]
    
    # Populate the array
    for (i, conrange) in enumerate(constraint)
        min_val, max_val = log10(conrange.min), log10(conrange.max)
        rand_vals = exp10.(rand(Uniform(min_val, max_val), n))
        
        for j in 1:n
            population[j][i] = rand_vals[j]
        end
    end
    
    return population
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
    population = [rand(Uniform(conrange.min, conrange.max), n) for conrange in constraint if conrange.active]
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

"""Returns the `fitness function(input)` for the cost function, referencing the GA problem with closure"""
    function make_fitness_function(gaprob::GAProblem)
        fixed_idxs = get_fixed_indices(gaprob.constraints)
        
        function fitness_function(input::Vector{Float64})
            # Preallocate a new input array, merging fixed and variable parameters
            merged_input = Vector{Float64}(undef, length(input) + length(fixed_idxs))
            
            var_idx = 1
            for idx in eachindex(merged_input)
                if idx in fixed_idxs
                    merged_input[idx] = gaprob.constraints[idx].fixed_value
                else
                    merged_input[idx] = input[var_idx]
                    var_idx += 1
                end
            end
            
            return evalfunc(merged_input, gaprob.ode_problem)
        end
        return fitness_function
    end
#> END


# """### Callback function that terminates the GA if the number of oscillations exceeds the threshold, and updates the progress bar"""
# function ga_callback(trace::Evolutionary.OptimizationTrace, progressbar::Progress, threshold::Int)
#     #? Callback function for the GA, updating the progress bar
#     num_oscillation = trace[end].metadata["num_oscillatory"]
#     if num_oscillation >= threshold 
#         finish!(progressbar)
#         return true
#     else
#         next!(progressbar, step = num_oscillation)
#         return false
#     end
# end


#< RUN GENETIC ALGORITHM OPTIMIZATION ##
"""
Runs the genetic algorithm, returning the `result`, and the `record` named tuple
"""
function run_GA(ga_problem::GAProblem, fitnessfunction_factory::Function=make_fitness_function; 
                                            population_size = 5000, abstol=1e-4, reltol=1e-2, successive_f_tol = 2, iterations=5, parallelization = :thread, show_trace=true)#, threshold=10000)
    # blas_threads = BLAS.get_num_threads()
    # BLAS.set_num_threads(1)

    #* Generate the initial population.
    pop = generate_population(ga_problem.constraints, population_size)

    #* Create constraints using the min and max values from constraints if they are active for optimization.
    boxconstraints = BoxConstraints([constraint.min for constraint in ga_problem.constraints if constraint.active], [constraint.max for constraint in ga_problem.constraints if constraint.active])

    # *Create Progress bar and callback function
    # ga_progress = Progress(threshold; desc = "GA Progress")
    # callback_func = (trace) -> ga_callback(trace, ga_progress, threshold)

    #* Define options for the GA.
    opts = Evolutionary.Options(abstol=abstol, reltol=reltol, successive_f_tol = successive_f_tol, iterations=iterations, 
                        store_trace = true, show_trace=show_trace, show_every=1, parallelization=parallelization)#, callback=callback_func)

    #* Define the range of possible values for each parameter when mutated, and the mutation scalar.

    #? BGA mutation scheme
    mutation_scalar = 0.5
    mutation_range = fill(mutation_scalar, length(ga_problem.constraints))
    mutation_scheme = BGA(mutation_range, 2)

    #? PM mutation scheme
    # lowerbound = [constraint.min/10 for constraint in ga_problem.constraints.ranges]
    # upperbound = [constraint.max*10 for constraint in ga_problem.constraints.ranges]
    # mutation_scheme = PM(lowerbound, upperbound, 2.)


    #* Define the GA method.
    mthd = GA(populationSize = population_size, selection = tournament(Int(population_size/10), select=argmax),
                crossover = TPX, crossoverRate = 1.0, # Two-point crossover event
                mutation  = mutation_scheme, mutationRate = 1.0)

    #* Make fitness function. Makes closure of evaluation function and ODE problem
    fitness_function = fitnessfunction_factory(ga_problem.eval_function, ga_problem.ode_problem)

    #* Run the optimization.
    result = Evolutionary.optimize(fitness_function, [0.0,0.0,0.0], boxconstraints, mthd, pop, opts)


    # BLAS.set_num_threads(blas_threads)
    # return result
    return GAResults(result, length(ga_problem.constraints))
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

function GAResults(result::Evolutionary.EvolutionaryOptimizationResults, indlength::Int)
    numpoints = sum(length, (gen.metadata["fitvals"] for gen in result.trace))
    population = [Vector{Float64}(undef, indlength) for i in 1:numpoints]
    fitvals = Vector{Float64}(undef, numpoints)
    periods = Vector{Float64}(undef, numpoints)
    amplitudes = Vector{Float64}(undef, numpoints)

    startidx = 1
    for gen in result.trace
        endidx = startidx + length(gen.metadata["population"]) - 1

        population[startidx:endidx] .= gen.metadata["population"]
  
        fitvals[startidx:endidx] .= gen.metadata["fitvals"]
     
        periods[startidx:endidx] .= gen.metadata["periods"]
    
        amplitudes[startidx:endidx] .= gen.metadata["amplitudes"]

        startidx = endidx + 1
    end
    return GAResults(result.trace, population, fitvals, periods, amplitudes)
end



"""Makes a DataFrame from the results of a GA optimization"""
function make_ga_dataframe(results::GAResults, constraints::ConstraintType)
    df = DataFrame(fit = results.fitvals, per = results.periods, amp = results.amplitudes)
    for (i,conrange) in enumerate(constraints)
        if conrange.active
            df[!, conrange.name] .= [x[i] for x in results.population]
        else
            df[!, conrange.name] .= conrange.fixed_value
        end
    end
    return df
end


# function make_ga_dataframe(results::GAResults, prob::ODEProblem)
#     df = DataFrame(ind = results.population, fit = results.fitvals, per = results.periods, amp = results.amplitudes)
#     split_dataframe!(df, prob)
#     return df
# end

# function make_ga_dataframe(results::GAResults, prob::ODEProblem, fixedDF::Float64)
#     df = DataFrame(ind = results.population, fit = results.fitvals, per = results.periods, amp = results.amplitudes)
#     split_dataframe!(df, prob, fixedDF)
#     return df
# end

# function make_ga_dataframe(results::GAResults, prob::ODEProblem, fixedval_idxs::Vector{Int}, fixedvals::Vector{Float64}, fixedDF::Float64)
#     df = DataFrame(fit = results.fitvals, per = results.periods, amp = results.amplitudes)
#     split_dataframe!(df, prob, fixedDF)
#     insertcols!(df, fixedval_idxs, fixedvals)
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

function split_dataframe!(df, prob, fixedDF)
    paramsymbols = [:ka1,:kb1,:kcat1,:ka2,:kb2,:ka3,:kb3,:ka4,:kb4,:ka7,:kb7,:kcat7]
    
    if length(df.ind[1]) == 12

        #* Split ind column into separate columns for each parameter
        for (i,param) in enumerate(paramsymbols)
            df[!, param] .= [x[i] for x in df.ind]
        end

        df[!, :DF] .= fixedDF

        select!(df, Not(:ind))

        #* Add initial conditions
        df.L .= prob.u0[1]
        df.K .= prob.u0[2]
        df.P .= prob.u0[3]
        df.A .= prob.u0[4]
    elseif length(df.ind[1]) == 4
        #* Add parameters
        for (i,param) in paramsymbols
            df[!, param] .= prob.p[i]
        end
        concsymbols = [:L,:K,:P,:A]
        for (i,conc) in enumerate(concsymbols)
            df[!, conc] .= [x[i] for x in df.ind]
        end

        df[!, :DF] .= fixedDF
        select!(df, Not(:ind))
    else 
        concsymbols = [:L,:K,:P,:A]
        allsymbols = vcat(paramsymbols, concsymbols)

        for (i,conc) in enumerate(allsymbols)
            df[!, conc] .= [x[i] for x in df.ind]
        end

        df[!, :DF] .= fixedDF
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
