#< CONSTRAINT RANGE STRUCTS

abstract type ConstraintSet end

"""
    iterate(constraint::ConstraintSet)

Iterate over the constraint ranges in a `ConstraintSet` object.
"""
# Define the length method
Base.length(constraint::ConstraintSet) = length(fieldnames(typeof(constraint)))

# Define the getindex method for index-based access
function Base.getindex(constraint::ConstraintSet, idx::Int)
    field_name = fieldnames(typeof(constraint))[idx]
    return getfield(constraint, field_name)
end

function find_field_index(field_name::Union{Symbol, String}, constraint::ConstraintSet)
    fields = fieldnames(typeof(constraint))
    idx = findfirst(x -> x == Symbol(field_name), fields)
    
    if idx === nothing
        throw(ArgumentError("Field name '$field_name' not found in ConstraintSet."))
    end
    
    return idx
end


# To make it iterable, define start, next and done methods
Base.iterate(constraint::ConstraintSet, state=1) = state > length(constraint) ? nothing : (getfield(constraint, fieldnames(typeof(constraint))[state]), state + 1)

# Required for the `in` keyword
Base.eltype(::Type{ConstraintSet}) = ConstraintRange

# Required for length
activelength(constraint::ConstraintSet) = length(filter(x -> !isfixed(constraint[x]), fieldnames(typeof(constraint))))

function get_fixed_indices(constraints::ConstraintSet)
    inactive_indices = Int[]  # Initialize an empty array to store the indices of inactive elements
    for (idx, constraint) in enumerate(constraints)  # Loop through each element along with its index
        if isfixed(constraint)  # Check if the element is fixed
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
- `fixed_value::Float64`: Fixed value is a float if inactive. 0.0 if it's to be optimized
"""
mutable struct ConstraintRange
    const name::Symbol
    const min::Float64
    const max::Float64
    fixed_value::Float64
end

"Tests whether a ConstraintRange is fixed or active. Will return true if fixed."
function isfixed(conrange::ConstraintRange)
    return !isnan(conrange.fixed_value)
end

# """Uses Setfield.jl to set the active field of a ConstraintRange"""
# function fix_constraint!(conrange::ConstraintRange, value::Float64)
#     @set conrange.fixed_value = value
# end

"""
    ParameterConstraints

Struct encapsulating parameter constraints. Each field represents a different parameter, holding a `ConstraintRange` object that defines the valid range for that parameter.
"""
struct ParameterConstraints <: ConstraintSet
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
struct InitialConditionConstraints <: ConstraintSet 
    L::ConstraintRange
    K::ConstraintRange
    P::ConstraintRange
    A::ConstraintRange
end

"""
    AllConstraints

Struct encapsulating all constraints. Each field represents a different parameter or initial condition, holding a `ConstraintRange` object that defines the valid range for that parameter or initial condition.
"""
struct AllConstraints <: ConstraintSet
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
function ParameterConstraints(; karange = (1e-3, 1e1), kbrange = (1e-3, 1e3), kcatrange = (1e-3, 1e3), dfrange = (1e3, 1e5))#, nominalvals = repeat([Nothing],13))
    #* Define parameter constraint ranges
    ka_min, ka_max = karange  # uM^-1s^-1, log scale
    kb_min, kb_max = kbrange  # s^-1, log scale
    kcat_min, kcat_max = kcatrange # s^-1, log scale
    df_min, df_max = dfrange # for DF, log scale

    return ParameterConstraints(
        ConstraintRange(:ka1, ka_min, ka_max, NaN),
        ConstraintRange(:kb1, kb_min, kb_max, NaN),
        ConstraintRange(:kcat1, kcat_min, kcat_max, NaN),
        ConstraintRange(:ka2, ka_min, ka_max, NaN),
        ConstraintRange(:kb2, kb_min, kb_max, NaN),
        ConstraintRange(:ka3, ka_min, ka_max, NaN),
        ConstraintRange(:kb3, kb_min, kb_max, NaN),
        ConstraintRange(:ka4, ka_min, ka_max, NaN),
        ConstraintRange(:kb4, kb_min, kb_max, NaN),
        ConstraintRange(:ka7, ka_min, ka_max, NaN),
        ConstraintRange(:kb7, kb_min, kb_max, NaN),
        ConstraintRange(:kcat7, kcat_min, kcat_max, NaN),
        ConstraintRange(:DF, df_min, df_max, NaN)
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
function InitialConditionConstraints(;Lrange = (0.1, 10.0), Krange = (0.1, 5.0), Prange = (0.1, 5.0), Arange = (0.1, 10.0))#, nominalvals = repeat([Nothing],4))
    # Define initial condition constraint ranges
    lipid_min, lipid_max = Lrange  # uM
    kinase_min, kinase_max = Krange  # uM
    phosphatase_min, phosphatase_max = Prange # uM
    ap2_min, ap2_max = Arange # uM

    return InitialConditionConstraints(
        ConstraintRange(:L, lipid_min, lipid_max, NaN),
        ConstraintRange(:K, kinase_min, kinase_max, NaN),
        ConstraintRange(:P, phosphatase_min, phosphatase_max, NaN),
        ConstraintRange(:A, ap2_min, ap2_max, NaN)
    )
end



function AllConstraints(paramconstraints::ParameterConstraints, icconstraints::InitialConditionConstraints) 
    return AllConstraints(
        paramconstraints.ka1,
        paramconstraints.kb1,
        paramconstraints.kcat1,
        paramconstraints.ka2,
        paramconstraints.kb2,
        paramconstraints.ka3,
        paramconstraints.kb3,
        paramconstraints.ka4,
        paramconstraints.kb4,
        paramconstraints.ka7,
        paramconstraints.kb7,
        paramconstraints.kcat7,
        paramconstraints.DF,

        icconstraints.L,
        icconstraints.K,
        icconstraints.P,
        icconstraints.A
    )
end

#> END




"""Fitness function constructor called during GAProblem construction that captures the fixed indices and ODE problem"""
# function make_fitness_function(constraints::ConstraintSet, ode_problem::ODEProblem, eval_function)

#     let fixed_idxs = get_fixed_indices(constraints), ode_problem = ode_problem, eval_function = eval_function
        
#         function fitness_function(input::Vector{Float64})
#             #* Preallocate a new input array, merging fixed and variable parameters
#             merged_input = Vector{Float64}(undef, length(input) + length(fixed_idxs))
            
#             var_idx = 1
#             for idx in eachindex(merged_input)
#                 if idx in fixed_idxs #* If the index is fixed in the constraints, assign that value at the same index in the merged input
#                     merged_input[idx] = constraints[idx].fixed_value
#                 else #* Otherwise, assign the value from the input array
#                     merged_input[idx] = input[var_idx]
#                     var_idx += 1
#                 end
#             end
            
#             return eval_function(merged_input, ode_problem)
#         end
#         return fitness_function
#     end
# end

function make_fitness_function(constraints::ConstraintSet, ode_problem::ODEProblem, eval_function)
    fixed_idxs = get_fixed_indices(constraints)
    fixed_values = [constraints[i].fixed_value for i in fixed_idxs]
    n_fixed = length(fixed_idxs)
    n_total = n_fixed + length(ode_problem.p)  # Assuming ode_problem.p contains the initial variable parameters

    function fitness_function(input::Vector{Float64})
        merged_input = Vector{Float64}(undef, n_total)
        merged_input[fixed_idxs] .= fixed_values  # Fill in fixed values

        # Fill in variable values
        j = 1
        for i in 1:n_total
            if !(i in fixed_idxs)
                merged_input[i] = input[j]
                j += 1
            end
        end

        return eval_function(merged_input, ode_problem)
    end

    return fitness_function
end

make_fitness_function(constraints::ParameterConstraints, ode_problem::ODEProblem) = make_fitness_function(constraints, ode_problem, eval_param_fitness)
make_fitness_function(constraints::InitialConditionConstraints, ode_problem::ODEProblem) = make_fitness_function(constraints, ode_problem, eval_ic_fitness)
make_fitness_function(constraints::AllConstraints, ode_problem::ODEProblem) = make_fitness_function(constraints, ode_problem, eval_all_fitness)


#< GA PROBLEM TYPE
"""
    GAProblem{T <: ConstraintSet}

Struct encapsulating a Genetic Algorithm (GA) optimization problem. It holds the constraints for the problem, the ODE problem to be solved, the fitness function, and any additional keyword arguments.

# Fields
- `constraints::T`: Constraints for the problem. Either `ParameterConstraints` or `InitialConditionConstraints` or `AllConstraints`.
- `ode_problem::ODEProblem`: ODE problem to be solved.
- `fitness_function::Function`: Fitness function, automatically generated with constructor
"""
mutable struct GAProblem{CT <: ConstraintSet, OT <: ODEProblem, FT <: Function}
    constraints::CT
    ode_problem::OT
    fitness_function::FT

    function GAProblem(constraints::CT, ode_problem::OT, fitness_function::FT = make_fitness_function(constraints, ode_problem)) where {CT<:ConstraintSet, OT<:ODEProblem, FT<:Function}
        new{CT,OT,FT}(constraints, ode_problem, fitness_function)
    end
end

# """Fixes constraints prior to GAProblem"""
# function set_fixed_constraints!(constraints::ConstraintSet; fixedinputs...)
#     for (name, value) in pairs(fixedinputs)
#         symbolic_name = Symbol(name)
#         if symbolic_name in fieldnames(typeof(constraints))
#             fix_constraint!(getfield(constraints, symbolic_name), value) #! THIS DOESN'T WORK FOR SOME REASON IN SETTING NESTED IMMUTABLE STRUCT PROPERTIES
#         end
#     end
#     return constraints
# end

function set_fixed_constraints!(constraints::ConstraintSet, fixedinputs)
    for (name, value) in fixedinputs
        symbolic_name = Symbol(name)
        if symbolic_name in fieldnames(typeof(constraints))
            conrange = getfield(constraints, symbolic_name)
            conrange.fixed_value = value
        end
    end
    return constraints
end

"""Sets the fixed values for a GAProblem and reconstructs the fitness function"""
function set_fixed_constraints!(ga_problem::GAProblem; fixedinputs...)
    for (name, value) in pairs(fixedinputs)
        symbolic_name = Symbol(name)
        if symbolic_name in fieldnames(typeof(constraints))
            setfield!(constraints[symbolic_name], :fixed_value, value)

        end
    end
    ga_problem.fitness_function = make_fitness_function(ga_problem.constraints, ga_problem.ode_problem)
    return ga_problem
end


function Base.show(io::IO, ::MIME"text/plain", prob::GAProblem) 
    printstyled(io, typeof(prob), ":\n"; bold = true, underline=true, color = :green)
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

Generate a population of `n` individuals for the given generic `constraints <: ConstraintSet`. Each individual is sampled from a log-uniform distribution within the valid range for each parameter or initial condition.

# Example
```julia
constraints = define_parameter_constraints(karange = (-3.0, 1.0), kbrange = (-3.0, 3.0))
population = generate_population(constraints, 100)
```
"""
# function generate_population(constraint::ConstraintSet, n::Int)
#     population = [exp10.(rand(Uniform(log10(conrange.min), log10(conrange.max)), n)) for conrange in constraint]
#     population = transpose(hcat(population...))
#     return [population[:, i] for i in 1:n]
# end

function generate_empty_population(constraint::ConstraintSet, n::Int)
    num_params = activelength(constraint)
    
    # Preallocate the population array of arrays
    [Vector{Float64}(undef, num_params) for _ in 1:n]
end

function generate_population(constraint::ConstraintSet, n::Int)
    # Preallocate the population array of arrays
    population = generate_empty_population(constraint, n)
    
    generate_population!(population, constraint)
end


function generate_population!(population::Vector{Vector{Float64}}, constraint::ConstraintSet)

    rand_vals = Vector{Float64}(undef, length(population))
    
    # Populate the array
    for (i, conrange) in enumerate(constraint)
        if !isfixed(conrange)
            min_val, max_val = log10(conrange.min), log10(conrange.max)
            rand_vals .= exp10.(rand(Uniform(min_val, max_val), length(population)))
            
            for j in eachindex(population)
                population[j][i] = rand_vals[j]
            end
        end
    end
    return population
end

#! MATRIX VERSIONS
# Generate an empty population as a 2D array
# function generate_empty_population(constraint::ConstraintSet, n::Int)
#     num_params = activelength(constraint)
#     return Matrix{Float64}(undef, num_params, n)
# end

# # Generate a population as a 2D array
# function generate_population(constraint::ConstraintSet, n::Int)
#     population = generate_empty_population(constraint, n)
#     generate_population!(population, constraint)
#     return population
# end

# # Populate the 2D array in-place
# function generate_population!(population::Matrix{Float64}, constraint::ConstraintSet)
#     n = size(population, 2)
    
#     for (i, conrange) in enumerate(constraint)
#         if conrange.active
#             min_val, max_val = log10(conrange.min), log10(conrange.max)
#             population[i, :] .= exp10.(rand(Uniform(min_val, max_val), n))
#         end
#     end
# end
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
        
        return gaprob.eval_function(merged_input, gaprob.ode_problem)
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
function run_GA(ga_problem::GAProblem, population::Vector{Vector{Float64}} = generate_population(ga_problem.constraints, 10000); 
                abstol=1e-4, reltol=1e-2, successive_f_tol = 2, iterations=5, parallelization = :thread, show_trace=true)#, threshold=10000)
    # blas_threads = BLAS.get_num_threads()
    # BLAS.set_num_threads(1)

    population_size = length(population)

    # #* Generate the initial population.
    # pop = generate_population!(population, ga_problem.constraints)

    #* Create constraints using the min and max values from constraints if they are active for optimization.
    boxconstraints = BoxConstraints([constraint.min for constraint in ga_problem.constraints if !isfixed(constraint)], [constraint.max for constraint in ga_problem.constraints if !isfixed(constraint)])

    # *Create Progress bar and callback function
    # ga_progress = Progress(threshold; desc = "GA Progress")
    # callback_func = (trace) -> ga_callback(trace, ga_progress, threshold)

    #* Define options for the GA.
    opts = Evolutionary.Options(abstol=abstol, reltol=reltol, successive_f_tol = successive_f_tol, iterations=iterations, 
                        store_trace = true, show_trace=show_trace, show_every=1, parallelization=parallelization)#, callback=callback_func)

    #* Define the range of possible values for each parameter when mutated, and the mutation scalar.

    #? BGA mutation scheme
    mutation_scalar = 0.5
    mutation_range = fill(mutation_scalar, activelength(ga_problem.constraints))
    mutation_scheme = BGA(mutation_range, 2)

    #? PM mutation scheme
    # lowerbound = [constraint.min/10 for constraint in ga_problem.constraints.ranges]
    # upperbound = [constraint.max*10 for constraint in ga_problem.constraints.ranges]
    # mutation_scheme = PM(lowerbound, upperbound, 2.)


    #* Define the GA method.
    mthd = GA(populationSize = population_size, selection = tournament(Int(population_size/10), select=argmax),
                crossover = TPX, crossoverRate = 1.0, # Two-point crossover event
                mutation  = mutation_scheme, mutationRate = 1.0)



    #* Run the optimization.
    result = Evolutionary.optimize(ga_problem.fitness_function, zeros(3,population_size), boxconstraints, mthd, population, opts)


    # BLAS.set_num_threads(blas_threads)
    # return result
    return GAResults(result, activelength(ga_problem.constraints))
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
function make_ga_dataframe(results::GAResults, constraints::ConstraintSet)
    df = DataFrame(fit = results.fitvals, per = results.periods, amp = results.amplitudes)
    for (i,conrange) in enumerate(constraints)
        if !isfixed(conrange)
            df[!, conrange.name] .= [x[i] for x in results.population]
        else
            df[!, conrange.name] .= conrange.fixed_value
        end
    end
    return df
end




#< MISCELLANEOUS FUNCTIONS ##
"""Defines logspace function for sampling parameters"""
logrange(start, stop, length::Int) = exp10.(collect(range(start=log10(start), stop=log10(stop), length=length)))






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
