begin 
    using Plots; #theme(:juno)
    using Catalyst
    using DifferentialEquations
    using Statistics
    using Peaks
    using Evolutionary, FFTW
    using Random
    using Distributions
    using DataFrames
    using CSV
    # using Unitful
    # using Unitful: µM, M, nm, µm, s, μs, Na, L, 𝐍
    using StaticArrays
    using BenchmarkTools, Profile, ProgressMeter
    # using Cthulhu
    # using JET
    # using MultivariateStats, UMAP, TSne, StatsPlots
    # using GlobalSensitivity, QuasiMonteCarlo
    using LinearAlgebra
    # using BifurcationKit, Setfield, ForwardDiff, Parameters; const BK = BifurcationKit
    using ColorSchemes, Plots.PlotMeasures
    # using OrderedCollections
    # using Combinatorics
    # using LazySets, Polyhedra
    default(lw = 2, size = (1000, 600), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)
    # plotlyjs()
    # import CairoMakie as cm 
    # gr()
    # push!(LOAD_PATH, "../../UTILITIES")

    include("../../UTILITIES/EvolutionaryOverloads.jl")

    # import the Catalyst model "fullrn"
    include("../../UTILITIES/ReactionNetwork.jl")

    # import the cost function and other evaluation functions
    include("../../UTILITIES/EvaluationFunctions.jl")
    # using .EvaluationFunctions

    # import the genetic algorithm and associated functions
    include("../../UTILITIES/GA_functions.jl")

    include("../../UTILITIES/ODEProbMaker.jl")

    include("../../UTILITIES/UnitTools.jl")


    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))
end


fullrn = make_fullrn()
ogprob = ODEProblem(fullrn, [], (0.,100.0), [])
# @benchmark solve($ogprob, saveat = 0.1, save_idxs = 1)

# ogsol = solve(ogprob, Rosenbrock23(), save_idxs = 1)
# plot(ogsol)


#* Optimization of parameters to produce data for CSV
param_constraints = define_parameter_constraints(karange = (1e-3, 1e2), kbrange = (1e-3, 1e3), kcatrange = (1e-2, 1e3), dfrange = (1e3, 2e4))


# param_gaproblem = GAProblem(param_constraints,ogprob)

# param_record = run_GA(param_gaproblem)
# CSV.write("OscillatorPaper/FigureGenerationScripts/initialconditions.csv", param_record)




# """Compare the allowed solution space for when each triplet of parameters is fixed"""
# function reachability_analysis(constraints::ConstraintType, prob::ODEProblem) 
#     #* Get all combinations of fixed pairs
#     fixed_triple_combos = combinations(constraints.ranges, 3)

#     combos_length = length(collect(fixed_triple_combos))
#     combos_names = [string(triplet[1].name, "_", triplet[2].name, "_", triplet[3].name) for triplet in fixed_triple_combos]

#     loopprogress = Progress(length(combos_length), desc ="Looping thru fixed pairs: " , color=:blue)

#     #* Make a results DataFrame where fixedpair => (volume, avg_period, avg_amplitude, num_oscillatory_points)
#     results_df = DataFrame(FixedTriplet = combos_names, volume = Vector{Float64}(undef, combos_length), periodstats = Vector{Vector{Float64}}(undef, combos_length), amplitudestats = Vector{Vector{Float64}}(undef, combos_length), num_oscillatory_points = Vector{Int}(undef, combos_length))

#     for (i,fixedpair) in enumerate(fixed_triple_combos)
#         @info "Fixed input pair: $(fixedpair[1].name), $(fixedpair[2].name)"

#         #* Make a copy of the constraints to modify
#         variable_constraints = deepcopy(constraints)

#         #* Remove the fixed parameters from the variable constraints
#         filter!(x -> x.name != fixedpair[1].name && x.name != fixedpair[2].name, variable_constraints.ranges)

#         fixed_ga_problem = GAProblem(variable_constraints, prob)

#         fixedpair_idxs = find_indices(fixedpair, constraints.ranges) # Get the indices of the fixed inputs.

#         #* Create a closure for the fitness function that includes the fixed inputs
#         make_fitness_function_closure(evalfunc,prob) = make_fitness_function_with_fixed_inputs(evalfunc, prob, fixedpair, fixedpair_idxs)

#         #* Run the optimization function to get the evaluated points
#         oscillatory_points_df = run_GA(fixed_ga_problem, make_fitness_function_closure; population_size = 10000, iterations = 8) #TODO: outputting the same number of points for multiple pairs

 
#         #* Calculate the number of oscillatory points
#         num_points = length(oscillatory_points_df.ind)
       

#         # #* Calculate the average period and amplitude
#         # periodstats = quantile(oscillatory_points_df.per, [0.0, 0.25, 0.5, 0.75, 1.0])
        

#         # amplitudestats = quantile(oscillatory_points_df.amp, [0.0, 0.25, 0.5, 0.75, 1.0])
        

#         #* Store the results for the given fixed parameter combination
#         results_df[i, 2:end] .= (volume, oscillatory_points_df.per, oscillatory_points_df.amp, num_points)

#         next!(loopprogress)
#     end
#     return results_df
# end



# Modification to make_fitness_function_with_fixed_inputs function
function make_fitness_function_with_fixed_inputs(evalfunc::Function, prob::ODEProblem, fixed_input_triplet::Vector{Float64}, triplet_idxs::Tuple{Int, Int, Int})
    function fitness_function(input::Vector{Float64})
        # Create a new input vector that includes the fixed inputs.
        new_input = Vector{Float64}(undef, length(input) + length(fixed_input_triplet))

        # Keep track of the number of fixed inputs that have been inserted.
        fixed_inputs_inserted = 0

        for i in eachindex(new_input)
            if i in triplet_idxs
                # If the current index matches the index of a fixed input, insert the fixed input.
                new_input[i] = fixed_input_triplet[fixed_inputs_inserted + 1]
                fixed_inputs_inserted += 1
            else
                # Otherwise, insert the next value from the input vector.
                new_input[i] = input[i - fixed_inputs_inserted]
            end
        end

        return evalfunc(new_input, prob)
    end
    return fitness_function
end

# Modified make_fitness_function_with_fixed_inputs function
function make_fitness_function_with_fixed_inputs(evalfunc::Function, prob::ODEProblem, fixed_input_triplet::Vector{Float64}, triplet_idxs::Tuple{Int, Int, Int})

    idx_map = [i for i in 1:(length(prob.p))]
    j = 1
    for i in eachindex(fixed_input_triplet)
        while j in triplet_idxs
            j += 1
        end
        idx_map[j] = input[i]
        j += 1
    end

    function fitness_function(input::Vector{Float64})
        new_input = Vector{Float64}(undef, length(input) + length(fixed_input_triplet))
        
        for i in eachindex(new_input)
            if i in triplet_idxs
                new_input[i] = fixed_input_triplet[i]
            else
                new_input[i] = input[idx_map[i]]
            end
        end
        return evalfunc(new_input, prob)
    end
    return fitness_function
end



# Modification to fixed_triplet_csv_maker function
function fixed_triplet_csv_maker(param1::String, param2::String, param3::String, constraints::ConstraintType, prob::ODEProblem)
    variable_constraints = deepcopy(constraints)
    fixedtrip = [x for x in variable_constraints.ranges if x.name == param1 || x.name == param2 || x.name == param3]
    filter!(x -> x.name != param1 && x.name != param2 && x.name != param3, variable_constraints.ranges)

    fixed_ga_problem = GAProblem(variable_constraints, prob)
    fixedtrip_idxs = find_indices(param1, param2, param3, constraints.ranges) 

    fixed_values1 = range(fixedtrip[1].min, stop = fixedtrip[1].max, length = 5)
    fixed_values2 = range(fixedtrip[2].min, stop = fixedtrip[2].max, length = 5)
    fixed_values3 = range(fixedtrip[3].min, stop = fixedtrip[3].max, length = 5)

    results_df = DataFrame(fixed_value1 = Float64[], fixed_value2 = Float64[], fixed_value3 = Float64[], average_period = Float64[], average_amplitude = Float64[])

    for val1 in fixed_values1
        for val2 in fixed_values2
            for val3 in fixed_values3
                fixed_values = [val1, val2, val3]
                @info fixed_values
                make_fitness_function_closure(evalfunc,prob) = make_fitness_function_with_fixed_inputs(evalfunc, prob, fixed_values, fixedtrip_idxs)
                oscillatory_points_df = run_GA(fixed_ga_problem, make_fitness_function_closure; population_size = 1000, iterations = 5) 

                average_period::Float64 = mean(oscillatory_points_df.per)
                average_amplitude::Float64 = mean(oscillatory_points_df.amp)

                push!(results_df, (val1, val2, val3, average_period, average_amplitude))
            end
        end
    end
    return results_df
end

results_df = fixed_triplet_csv_maker("kcat1", "kcat7", "DF", param_constraints, ogprob)
results_df = DataFrame(kcat1 = results_df.fixed_value1, kcat7 = results_df.fixed_value2, DF = results_df.fixed_value3, average_period = results_df.average_period, average_amplitude = results_df.average_amplitude)
CSV.write("OscillatorPaper/FigureGenerationScripts/fixed_triplet_results.csv", results_df)

@code_warntype fixed_triplet_csv_maker("ka1", "ka2", "ka3", param_constraints, ogprob)