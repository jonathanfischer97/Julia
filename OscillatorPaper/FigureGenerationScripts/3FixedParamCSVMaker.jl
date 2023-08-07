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
ogprob = ODEProblem(fullrn, [], (0.,500.0), [])
new_u0 = ogprob.u0 .* 10
ogprob = remake(ogprob, u0 = new_u0)
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

    function fitness_function(input::Vector{Float64})
        # Create a new input vector that includes the fixed inputs.
        new_input = Vector{Float64}(undef, length(input) + length(fixed_input_triplet))

        # Pre-calculate the mappings of indices from `input` to `new_input`
        idx_map = [i for i in 1:(length(input) + length(fixed_input_triplet))]
        j = 1
        for i in eachindex(fixed_input_triplet)
            while j in triplet_idxs
                j += 1
            end
            idx_map[j] = i
            j += 1
        end

        # Keep track of the number of fixed inputs that have been inserted.
        fixed_inputs_inserted = 0

        for i in eachindex(new_input)
            if i in triplet_idxs
                # If the current index matches the index of a fixed input, insert the fixed input.
                new_input[i] = fixed_input_triplet[fixed_inputs_inserted + 1]
                fixed_inputs_inserted += 1
            else
                # Otherwise, insert the next value from the input vector.
                new_input[i] = input[idx_map[i]]
            end
        end

        return evalfunc(new_input, prob)
    end
    return fitness_function
end


"""Defines logspace function for sampling parameters"""
logrange(start, stop, length) = exp10.(collect(range(start=log10(start), stop=log10(stop), length=length)))

# logrange(start=1e-2, stop=1e2, length=10)


#* Modification to fixed_triplet_csv_maker function
function fixed_triplet_csv_maker(param1::String, param2::String, param3::String, constraints::ConstraintType, prob::ODEProblem; rangelength = 5)
    variable_constraints = deepcopy(constraints)
    fixedtrip = [x for x in variable_constraints.ranges if x.name == param1 || x.name == param2 || x.name == param3]
    filter!(x -> x.name != param1 && x.name != param2 && x.name != param3, variable_constraints.ranges)

    fixed_ga_problem = GAProblem(variable_constraints, prob)
    fixedtrip_idxs = find_indices(param1, param2, param3, constraints.ranges) 

    # fixed_values1 = range(fixedtrip[1].min, stop = fixedtrip[1].max, length = 5)
    fixed_values1 = logrange(fixedtrip[1].min, fixedtrip[1].max, rangelength)
    # fixed_values2 = range(fixedtrip[2].min, stop = fixedtrip[2].max, length = 5)
    fixed_values2 = logrange(fixedtrip[2].min, fixedtrip[2].max, rangelength)
    # fixed_values3 = range(fixedtrip[3].min, stop = fixedtrip[3].max, length = 5)
    fixed_values3 = logrange(fixedtrip[3].min, fixedtrip[3].max, rangelength)

    results_df = DataFrame(param1 => Vector{Float64}(undef, rangelength^3), param2 => Vector{Float64}(undef, rangelength^3), param3 => Vector{Float64}(undef, rangelength^3), 
                        "average_period" => Vector{Float64}(undef, rangelength^3), "maximum_period"=>Vector{Float64}(undef, rangelength^3), "minimum_period"=>Vector{Float64}(undef, rangelength^3),
                        "average_amplitude" => Vector{Float64}(undef, rangelength^3), "maximum_amplitude"=>Vector{Float64}(undef, rangelength^3), "minimum_amplitude"=>Vector{Float64}(undef, rangelength^3))

    i = 1
    for val1 in fixed_values1
        for val2 in fixed_values2
            for val3 in fixed_values3
                fixed_values = [val1, val2, val3]
                @info fixed_values
                make_fitness_function_closure(evalfunc,prob) = make_fitness_function_with_fixed_inputs(evalfunc, prob, fixed_values, fixedtrip_idxs)
                oscillatory_points_df = run_GA(fixed_ga_problem, make_fitness_function_closure; population_size = 5000, iterations = 5) 

                if isempty(oscillatory_points_df)
                    continue
                else
                    average_period::Float64 = mean(oscillatory_points_df.per)
                    maximum_period::Float64 = maximum(oscillatory_points_df.per)
                    minimum_period::Float64 = minimum(oscillatory_points_df.per)
    
                    average_amplitude::Float64 = mean(oscillatory_points_df.amp)
                    maximum_amplitude::Float64 = maximum(oscillatory_points_df.amp)
                    minimum_amplitude::Float64 = minimum(oscillatory_points_df.amp)
    
                    # push!(results_df, (val1, val2, val3, average_period, maximum_period, minimum_period, average_amplitude, maximum_amplitude, minimum_amplitude))
                    results_df[i, :] = (val1, val2, val3, average_period, maximum_period, minimum_period, average_amplitude, maximum_amplitude, minimum_amplitude)
                    i += 1
                end
            end
        end
    end
    return results_df
end

param_triplet = ["kcat1", "kcat7", "DF"]
param_triplet_symbols = Symbol.(param_triplet)
results_df = fixed_triplet_csv_maker(param_triplet..., param_constraints, ogprob)
results_df = DataFrame(kb1 = results_df.fixed_value1, kb2 = results_df.fixed_value2, kb7 = results_df.fixed_value3, #TODO fix programmatic naming
                            average_period = results_df.average_period, maximum_period = results_df.maximum_period, minimum_period = results_df.minimum_period, 
                            average_amplitude = results_df.average_amplitude, maximum_amplitude = results_df.maximum_amplitude, minimum_amplitude = results_df.minimum_amplitude)
CSV.write("OscillatorPaper/FigureGenerationScripts/fixed_triplet_results-$(param_triplet[1]*param_triplet[2]*param_triplet[3]).csv", results_df)



@code_warntype fixed_triplet_csv_maker("ka1", "ka2", "ka3", param_constraints, ogprob)



#! TESTING GA PERFORMANCE
test_gaproblem = GAProblem(param_constraints, ogprob)


# @time "18 threads, 18 BLAS thread" run_GA(test_gaproblem)

# @time "18 threads, 1 BLAS threads" run_GA(test_gaproblem)
# "18 threads, 1 BLAS threads: 30.316618 seconds (19.58 M allocations: 19.422 GiB, 2.85% gc time)"

# Threads.nthreads()
# @time "36 threads, 1 BLAS thread" run_GA(test_gaproblem)
# "36 threads, 1 BLAS thread: 21.726124 seconds (15.02 M allocations: 19.014 GiB, 2.25% gc time)"


#! TESTING GA FUNCTIONALITY
test_results = run_GA(test_gaproblem; population_size = 5000, iterations = 5)
sort!(test_results, :fit, rev=false)


plotsol(row) = plotsol(row, test_results, ogprob)

plotsol(1)


test_fitness(row) = eval_param_fitness(test_results.ind[row], ogprob)
test_fitness(1)

# reogprob = remake(ogprob, p=test_results.ind[1])
# testsol = solve(reogprob, saveat = 0.1)
testsol = solve(remake(ogprob, p=test_results.ind[1]), saveat=0.01, save_idxs=1)
plot(testsol)
CostFunction(testsol)



#Test the fitness function with scaled amplitudes
testsolx10 = testsol.u .* 10

CostFunction(testsolx10)

function CostFunction(u)
    norm_u = normalize_time_series(u) #* normalize the solution
    #*get the fft of the solution
    fftData = getFrequencies(norm_u)
    fft_peakindexes, fft_peakvals = findmaxima(fftData,1) #* get the indexes of the peaks in the fft
    time_peakindexes, time_peakvals = findmaxima(u,1) #* get the times of the peaks in the fft
    if length(fft_peakindexes) < 3 || length(time_peakindexes) < 3 #* if there are no peaks in either domain, return 0
        return [0.0, 0.0, 0.0]
    end
    std = getSTD(fft_peakindexes, fftData) #* get the average standard deviation of the peaks in frequency domain
    diff = getDif(fft_peakvals) #* get the summed difference between the peaks in frequency domain

    #* Compute the period and amplitude
    # period, amplitude = getPerAmp(sol, time_peakindexes, time_peakvals)

    #* Return cost, period, and amplitude as a vector
    return -std - diff
end


getFrequencies(testsol.u)
getFrequencies(testsolx10)