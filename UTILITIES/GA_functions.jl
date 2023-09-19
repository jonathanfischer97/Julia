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

"""Gives the active length of a ConstraintSet, i.e. the number of elements that are not fixed"""
activelength(constraints::ConstraintSet) = count(x -> !x.isfixed, constraints)

"""Returns a vector of the numerical indices of the fixed elements in a ConstraintSet"""
function get_fixed_indices(constraints::ConstraintSet)::Vector{Int}
    inactive_indices = Int[]  # Initialize an empty array to store the indices of inactive elements
    for (idx, constraint) in enumerate(constraints)  # Loop through each element along with its index
        if constraint.isfixed  # Check if the element is fixed
            push!(inactive_indices, idx)  # Add the index to the array
        end
    end
    return inactive_indices  # Return the array of indices
end


"""Returns a vector of the constraintranges that are marked as fixed but have yet to be assigned fixed values"""
function get_fixed_constraintranges(constraints::ConstraintSet)::Vector{ConstraintRange}
    fixed_constraintranges = ConstraintRange[]
    for constraintrange in constraints  # Loop through each element along with its index
        if constraintrange.isfixed && isnan(constraintrange.fixed_value)  # Check if the element is fixed but not assigned a value
            push!(fixed_constraintranges, constraintrange)  # Add the index to the array
        end
    end
    return fixed_constraintranges  # Return the array of indices
end


"""
    ConstraintRange

Struct for defining parameter or initial condition ranges. Each instance contains a name, and a range defined by a minimum and maximum value.

# Fields
- `name::String`: Name of the parameter or initial condtion.
- `min::Float64`: Minimum allowed value.
- `max::Float64`: Maximum allowed value.
- `isfixed::Bool`: Whether the parameter or initial condition is fixed. Defaults to false.
- `fixed_value::Float64`: Fixed value is to be used if fixed. Defaults to NaN.
"""
@kwdef mutable struct ConstraintRange
    const name::Symbol
    const min::Float64
    const max::Float64
    isfixed::Bool = false
    fixed_value::Float64 = NaN
end


"""
    ParameterConstraints

Struct encapsulating parameter constraints. Each field represents a different parameter, holding a `ConstraintRange` object that defines the valid range for that parameter.
"""
mutable struct ParameterConstraints <: ConstraintSet
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
mutable struct InitialConditionConstraints <: ConstraintSet 
    L::ConstraintRange
    K::ConstraintRange
    P::ConstraintRange
    A::ConstraintRange
end

"""
    AllConstraints

Struct encapsulating all constraints. Each field represents a different parameter or initial condition, holding a `ConstraintRange` object that defines the valid range for that parameter or initial condition.
"""
mutable struct AllConstraints <: ConstraintSet
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
    ParameterConstraints(; kwargs...)

Define parameter constraints. Each keyword argument represents a different parameter, where the value is a tuple defining the valid range for that parameter.

# Example
```julia
constraints = ParameterConstraints(
    karange = (-3.0, 1.0), 
    kbrange = (-3.0, 3.0), 
    kcatrange = (-3.0, 3.0), 
    dfrange = (1.0, 5.0)
)
```
"""
function ParameterConstraints(; karange = (1e-3, 1e2), kbrange = (1e-3, 1e3), kcatrange = (1e-3, 1e3), dfrange = (1e2, 2e4))#, nominalvals = repeat([Nothing],13))
    #* Define parameter constraint ranges
    ka_min, ka_max = karange  # uM^-1s^-1, log scale
    kb_min, kb_max = kbrange  # s^-1, log scale
    kcat_min, kcat_max = kcatrange # s^-1, log scale
    df_min, df_max = dfrange # for DF, log scale

    return ParameterConstraints(
        ConstraintRange(name = :ka1, min = ka_min, max = ka_max),
        ConstraintRange(name = :kb1, min = kb_min, max = kb_max),
        ConstraintRange(name = :kcat1, min = kcat_min, max = kcat_max),
        ConstraintRange(name = :ka2, min = ka_min, max = ka_max),
        ConstraintRange(name = :kb2, min = kb_min, max = kb_max),
        ConstraintRange(name = :ka3, min = ka_min, max = ka_max),
        ConstraintRange(name = :kb3, min = kb_min, max = kb_max),
        ConstraintRange(name = :ka4, min = ka_min, max = ka_max),
        ConstraintRange(name = :kb4, min = kb_min, max = kb_max),
        ConstraintRange(name = :ka7, min = ka_min, max = ka_max),
        ConstraintRange(name = :kb7, min = kb_min, max = kb_max),
        ConstraintRange(name = :kcat7, min = kcat_min, max = kcat_max),
        ConstraintRange(name = :DF, min = df_min, max = df_max)
    )
end


"""
    InitialConditionConstraints(; kwargs...)

Define initial condition constraints. Each keyword argument represents a different initial condition, where the value is a tuple defining the valid range for that initial condition.

# Example
```julia
constraints = InitialConditionConstraints(
    lipidrange = (0.1, 10.0), 
    kinaserange = (0.1, 10.0), 
    phosphataserange = (0.1, 10.0), 
    ap2range = (0.1, 10.0)
)
```
"""
function InitialConditionConstraints(; Lrange = (1e-1, 1e2), Krange = (1e-2, 1e2), Prange = (1e-2, 1e2), Arange = (1e-1, 1e2))#, nominalvals = repeat([Nothing],4))
    # Define initial condition constraint ranges
    lipid_min, lipid_max = Lrange  # uM
    kinase_min, kinase_max = Krange  # uM
    phosphatase_min, phosphatase_max = Prange # uM
    ap2_min, ap2_max = Arange # uM

    return InitialConditionConstraints(
        ConstraintRange(name = :L, min = lipid_min, max = lipid_max),
        ConstraintRange(name = :K, min = kinase_min, max = kinase_max),
        ConstraintRange(name = :P, min = phosphatase_min, max = phosphatase_max),
        ConstraintRange(name = :A, min = ap2_min, max = ap2_max)
    )
end


function AllConstraints(paramconstraints::ParameterConstraints=ParameterConstraints(), icconstraints::InitialConditionConstraints=InitialConditionConstraints()) 
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
function make_fitness_function(constraints::ConstraintSet, ode_problem::OT, eval_function::FT) where {OT<:ODEProblem, FT<:Function}
    fixed_idxs = get_fixed_indices(constraints)
    fixed_values = [constraints[i].fixed_value for i in fixed_idxs]
    n_fixed = length(fixed_idxs)
    n_total = n_fixed + activelength(constraints) 

    non_fixed_indices = setdiff(1:n_total, fixed_idxs)

    # merged_input = Vector{Float64}(undef, n_total)
    merged_input = zeros(Float64, n_total+12)
    # @info "Merged input length: $(length(merged_input))"

    merged_input[fixed_idxs] .= fixed_values  # Fill in fixed values

    function fitness_function(input::Vector{Float64})
        merged_input[non_fixed_indices] .= input  # Fill in variable values
        # @info "Merged input: $merged_input"
        return eval_function(merged_input, ode_problem)
    end

    return fitness_function
end

make_fitness_function(constraints::ParameterConstraints, ode_problem::ODEProblem) = make_fitness_function(constraints, ode_problem, eval_param_fitness)
make_fitness_function(constraints::InitialConditionConstraints, ode_problem::ODEProblem) = make_fitness_function(constraints, ode_problem, eval_ic_fitness)
make_fitness_function(constraints::AllConstraints, ode_problem::ODEProblem) = make_fitness_function(constraints, ode_problem, eval_all_fitness)

# """Returns in-place function"""
# function make_fitness_function_inplace(constraints::ConstraintSet, ode_problem::OT, eval_function::FT) where {OT<:ODEProblem, FT<:Function}
#     fixed_idxs = get_fixed_indices(constraints)
#     fixed_values = [constraints[i].fixed_value for i in fixed_idxs]
#     n_fixed = length(fixed_idxs)
#     n_total = n_fixed + activelength(constraints)

#     non_fixed_indices = setdiff(1:n_total, fixed_idxs)


#     # Create a ThreadLocal array
#     merged_inputs = [zeros(Float64, n_total+12) for _ in 1:Threads.nthreads()]
#     Fvs = [Vector{Float64}(undef,3) for _ in 1:Threads.nthreads()]

#     # Fill in the fixed values
#     for input in merged_inputs
#         input[fixed_idxs] .= fixed_values  # Fill in fixed values
#     end

#     function fitness_function!(input::Vector{Float64})
#         # Get the merged_input array for the current thread
#         merged_input = merged_inputs[Threads.threadid()]
#         merged_input[non_fixed_indices] .= input  # Fill in variable values
#         Fvs[Threads.threadid()] .= eval_function(merged_input, ode_problem)
#     end

#     return fitness_function!
# end

# make_fitness_function_inplace(constraints::ParameterConstraints, ode_problem::ODEProblem) = make_fitness_function_inplace(constraints, ode_problem, eval_param_fitness)
# make_fitness_function_inplace(constraints::InitialConditionConstraints, ode_problem::ODEProblem) = make_fitness_function_inplace(constraints, ode_problem, eval_ic_fitness)
# make_fitness_function_inplace(constraints::AllConstraints, ode_problem::ODEProblem) = make_fitness_function_inplace(constraints, ode_problem, eval_all_fitness)


"""Multithreaded fitness function, allocated a merged array for each thread"""
function make_fitness_function_threaded(constraints::ConstraintSet, ode_problem::OT, eval_function::FT) where {OT<:ODEProblem, FT<:Function}
    fixed_idxs = get_fixed_indices(constraints)
    fixed_values = [constraints[i].fixed_value for i in fixed_idxs]
    n_fixed = length(fixed_idxs)
    n_total = n_fixed + activelength(constraints) 

    non_fixed_indices = setdiff(1:n_total, fixed_idxs)

    # Create a ThreadLocal array
    merged_inputs = [zeros(Float64, n_total+12) for _ in 1:Threads.nthreads()]

    # Fill in the fixed values
    for input in merged_inputs
        input[fixed_idxs] .= fixed_values  # Fill in fixed values
    end

    function fitness_function(input::Vector{Float64})
        # Get the merged_input array for the current thread
        merged_input = merged_inputs[Threads.threadid()]
        merged_input[non_fixed_indices] .= input  # Fill in variable values

        return eval_function(merged_input, ode_problem)
    end

    return fitness_function
end

make_fitness_function_threaded(constraints::ParameterConstraints, ode_problem::ODEProblem) = make_fitness_function_threaded(constraints, ode_problem, eval_param_fitness)
make_fitness_function_threaded(constraints::InitialConditionConstraints, ode_problem::ODEProblem) = make_fitness_function_threaded(constraints, ode_problem, eval_ic_fitness)
make_fitness_function_threaded(constraints::AllConstraints, ode_problem::ODEProblem) = make_fitness_function_threaded(constraints, ode_problem, eval_all_fitness)


#< GA PROBLEM TYPE
"""
    GAProblem{T <: ConstraintSet}

Struct encapsulating a Genetic Algorithm (GA) optimization problem. It holds the constraints for the problem, the ODE problem to be solved.

# Fields
- `constraints::T`: Constraints for the problem. Either `ParameterConstraints` or `InitialConditionConstraints` or `AllConstraints`.
- `ode_problem::ODEProblem`: ODE problem to be solved.
"""
@kwdef mutable struct GAProblem{CT <: ConstraintSet, OT <: ODEProblem}
    constraints::CT = AllConstraints()
    ode_problem::OT = make_ODE_problem()
end


"""Simply marks the constraints as fixed, without assigning a value"""
function set_fixed_constraints!(constraints::ConstraintSet, fixednames::Vector{Symbol})
    for name in fixednames
        if name in fieldnames(typeof(constraints))
            conrange = getfield(constraints, name)
            conrange.isfixed = true
            conrange.fixed_value = NaN
        end
    end
    return constraints
end

"""Sets the vector of unpacked fixed constraints according to symbol, assigning the given values"""
function set_fixed_values!(fixed_constraintranges::Vector{ConstraintRange}, values...)
    for (conrange, value) in zip(fixed_constraintranges, values)
        conrange.fixed_value = value
    end
    return fixed_constraintranges
end

"""Unsets the fixed constraints according to symbol, resetting both the isfixed and fixed_value fields to default"""
function unset_fixed_constraints!(constraints::ConstraintSet, fixednames::Vector{Symbol})
    for name in fixednames
        if name in fieldnames(typeof(constraints))
            conrange = getfield(constraints, name)
            conrange.isfixed = false
            conrange.fixed_value = NaN
        end
    end
    return constraints
end


function Base.show(io::IO, ::MIME"text/plain", prob::GAProblem) 
    printstyled(io, typeof(prob), ":\n"; bold = true, underline=true, color = :green)
    printstyled(io, prob.constraints, "\n")
    printstyled(io, "\nNominal parameter values:\n"; bold = true, color = :blue)
    printstyled(io, prob.ode_problem.p, "\n")
    printstyled(io, "\nNominal initial conditions:\n"; bold = true, color = :blue)
    printstyled(io, prob.ode_problem.u0, "\n")

    printstyled(io, "\nFixed values:\n"; bold = true, color = :red)
    printstyled(io, [(con.name => con.fixed_value) for con in prob.constraints if con.isfixed], "\n")
end
#> END

#< POPULATION GENERATION METHODS
"""
    generate_population(constraints::ConstraintSet, n::Int)

Generate a population of `n` individuals for the given generic `constraints <: ConstraintSet`. Each individual is sampled from a log-uniform distribution within the valid range for each parameter or initial condition.

# Example
```julia
constraints = ParameterConstraints()
population = generate_population(constraints, 100)
```
"""
function generate_population(constraints::ConstraintSet, n::Int)
    # Preallocate the population array of arrays
    population = generate_empty_population(constraints, n)
    
    generate_population!(population, constraints)
end

"""Generate an empty population array of arrays, where the length of each individual is the number of constraints minus the fixed ones"""
function generate_empty_population(constraints::ConstraintSet, n::Int)
    num_params = activelength(constraints)
    
    # Preallocate the population array of arrays
    [Vector{Float64}(undef, num_params) for _ in 1:n]
end

"""In-place population generation of already allocated empty population"""
function generate_population!(population::Vector{Vector{Float64}}, constraints::ConstraintSet)

    rand_vals = Vector{Float64}(undef, length(population))
    
    # Populate the array
    i = 1
    for conrange in constraints
        if !conrange.isfixed
            min_val, max_val = log10(conrange.min), log10(conrange.max)
            rand_vals .= exp10.(rand(Uniform(min_val, max_val), length(population)))
            
            for j in 1:length(population)
                population[j][i] = rand_vals[j]
            end
            i += 1
        end
    end
    return population
end

function generate_population_stratified!(population::Vector{Vector{Float64}}, constraints::ConstraintSet, n_strata::Int)

    rand_vals = Vector{Float64}(undef, length(population))
    
    # Populate the array
    i = 1
    for conrange in constraints
        if !conrange.isfixed
            min_val, max_val = log10(conrange.min), log10(conrange.max)
            strata_bounds = range(min_val, stop=max_val, length=n_strata+1)
            
            for stratum in 1:n_strata
                lower_bound, upper_bound = strata_bounds[stratum], strata_bounds[stratum+1]
                n_samples_stratum = length(population) ÷ n_strata  # Integer division to get an equal number of samples from each stratum
                n_samples_stratum += (stratum <= length(population) % n_strata)  # Distribute any remaining samples among the first few strata
                
                rand_vals[1:n_samples_stratum] .= exp10.(rand(Uniform(lower_bound, upper_bound), n_samples_stratum))
                
                for j in 1:n_samples_stratum
                    population[(stratum-1)*n_samples_stratum+j][i] = rand_vals[j]
                end
            end
            i += 1
        end
    end
    return population
end


# function generate_population_stratified!(population::Vector{Vector{Float64}}, constraints::ConstraintSet, n_strata::Int)
#     # Pre-allocate rand_vals based on n_strata and population size
#     rand_vals = Vector{Float64}(undef, length(population) ÷ n_strata + 1)
    
#     i = 1
#     for conrange in constraints
#         if !conrange.isfixed
#             min_val, max_val = log10(conrange.min), log10(conrange.max)
#             strata_bounds = range(min_val, stop=max_val, length=n_strata+1)
            
#             for stratum in 1:n_strata
#                 lower_bound, upper_bound = strata_bounds[stratum], strata_bounds[stratum+1]
#                 n_samples_stratum = length(population) ÷ n_strata
#                 n_samples_stratum += (stratum <= length(population) % n_strata)
                
#                 rand_vals[1:n_samples_stratum] .= exp10.(rand(Uniform(lower_bound, upper_bound), n_samples_stratum))
                
#                 for j in 1:n_samples_stratum
#                     population[(stratum-1)*n_samples_stratum+j][i] = rand_vals[j]
#                 end
#             end
#             i += 1
#         end
#     end
#     return population
# end



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
                abstol=1e-4, reltol=1e-2, successive_f_tol = 4, iterations=5, parallelization = :thread, show_trace=true)#, threshold=10000)
    # blas_threads = BLAS.get_num_threads()
    # BLAS.set_num_threads(1)

    population_size = length(population)

    #* Create constraints using the min and max values from constraints if they are active for optimization.
    boxconstraints = BoxConstraints([constraint.min for constraint in ga_problem.constraints if !constraint.isfixed], [constraint.max for constraint in ga_problem.constraints if !constraint.isfixed])

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

    #* Make fitness function
    fitness_function = make_fitness_function_threaded(ga_problem.constraints, ga_problem.ode_problem)

    #* Run the optimization.
    result = Evolutionary.optimize(fitness_function, zeros(3,population_size), boxconstraints, mthd, population, opts)

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
    gen_indices::Vector{Tuple{Int,Int}}
end

function GAResults(result::Evolutionary.EvolutionaryOptimizationResults, indlength::Int)
    numpoints = sum(length, (gen.metadata["fitvals"] for gen in result.trace))
    population = [Vector{Float64}(undef, indlength) for _ in 1:numpoints]
    fitvals = Vector{Float64}(undef, numpoints)
    periods = Vector{Float64}(undef, numpoints)
    amplitudes = Vector{Float64}(undef, numpoints)

    gen_indices = Tuple{Int, Int}[]
    startidx = 1
    for gen in result.trace
        endidx = startidx + length(gen.metadata["population"]) - 1

        push!(gen_indices, (startidx, endidx))

        population[startidx:endidx] .= gen.metadata["population"]
  
        fitvals[startidx:endidx] .= gen.metadata["fitvals"]
     
        periods[startidx:endidx] .= gen.metadata["periods"]
    
        amplitudes[startidx:endidx] .= gen.metadata["amplitudes"]

        startidx = endidx + 1
    end
    return GAResults(result.trace, population, fitvals, periods, amplitudes, gen_indices)
end

function save_to_csv(results::GAResults, constraints::ConstraintSet, filename::String)
    open(filename, "w") do io
        # Write the header with an additional "Generation" column
        write(io, "gen,fit,per,amp")
        for conrange in constraints
            write(io, ",$(conrange.name)")
        end
        write(io, "\n")
        
        # Loop over each generation based on gen_indices
        for (gen, (start_idx, end_idx)) in enumerate(results.gen_indices)
            for i in start_idx:end_idx
                # Write the generation, fitness, period, and amplitude values
                write(io, "$gen,$(results.fitvals[i]),$(results.periods[i]),$(results.amplitudes[i])")
                
                # Write the population and fixed values
                j = 1
                for conrange in constraints
                    if !conrange.isfixed
                        write(io, ",$(results.population[i][j])")
                        j += 1
                    else
                        write(io, ",$(conrange.fixed_value)")
                    end
                end
                
                write(io, "\n")
            end
        end
    end
end


"""Makes a DataFrame from the results of a GA optimization"""
function make_ga_dataframe(results::GAResults, constraints::ConstraintSet)
    df = DataFrame(gen = Vector{Int}(undef, length(results.fitvals)), fit = results.fitvals, per = results.periods, amp = results.amplitudes, relamp = Vector{Float64}(undef, length(results.fitvals)))

    #* Loop over each generation based on gen_indices
    for (gen, (start_idx, end_idx)) in enumerate(results.gen_indices)
        df.gen[start_idx:end_idx] .= gen
    end

    i = 1
    for conrange in constraints
        if !conrange.isfixed
            df[!, conrange.name] .= [x[i] for x in results.population]
            i+=1
        else
            df[!, conrange.name] .= conrange.fixed_value
        end
    end
    #* Calculate the relative amplitude by dividing the amp column by the initial concentration of A
    df.relamp .= df.amp ./ df.A[1]
    return df
end

"""Makes a DataFrame from a raw population, nested vectors"""
function make_pop_dataframe(pop::Vector{Vector{Float64}}, constraints::AllConstraints)
    df = DataFrame()
    i = 1
    for conrange in constraints
        if !conrange.isfixed
            df[!, conrange.name] = [x[i] for x in pop]
            i+=1
        else
            df[!, conrange.name] = conrange.fixed_value
        end
    end
    return df
end


#< MISCELLANEOUS FUNCTIONS ##
"""Defines logspace function for sampling parameters and initial conditions"""
logrange(start, stop, length::Int) = exp10.(collect(range(start=log10(start), stop=log10(stop), length=length)))
#> END
