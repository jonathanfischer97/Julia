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
    # using StaticArrays
    using BenchmarkTools, Profile, ProgressMeter
    # using Cthulhu
    # using JET
    # using MultivariateStats, UMAP, TSne, StatsPlots
    # using GlobalSensitivity, QuasiMonteCarlo
    using LinearAlgebra
    # using BifurcationKit, Setfield, ForwardDiff, Parameters; const BK = BifurcationKit
    # using ColorSchemes, Plots.PlotMeasures
    # using OrderedCollections
    # using Combinatorics
    # using LazySets, Polyhedra
    # default(lw = 2, size = (1000, 600), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)
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

param_constraints = define_parameter_constraints(ogprob; karange = (1e-3, 1e2), kbrange = (1e-2, 1e3), kcatrange = (1e-2, 1e3), dfrange = (1e3, 1e5))

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



"""Defines logspace function for sampling parameters"""
logrange(start, stop, length) = exp10.(collect(range(start=log10(start), stop=log10(stop), length=length)))



#* Modification to fixed_triplet_csv_maker function
function fixed_triplet_csv_maker(param1::String, param2::String, param3::String, constraints::ConstraintType, prob::ODEProblem; rangelength = 4) #TODO add progress bar
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

    results_df = DataFrame(param1 => Vector{Float64}(undef, rangelength^3), param2 => Vector{Float64}(undef, rangelength^3), param3 => Vector{Float64}(undef, rangelength^3), "num_oscillatory_points" => Vector{Int}(undef, rangelength^3), 
                        "average_period" => Vector{Float64}(undef, rangelength^3), "maximum_period"=>Vector{Float64}(undef, rangelength^3), "minimum_period"=>Vector{Float64}(undef, rangelength^3),
                        "average_amplitude" => Vector{Float64}(undef, rangelength^3), "maximum_amplitude"=>Vector{Float64}(undef, rangelength^3), "minimum_amplitude"=>Vector{Float64}(undef, rangelength^3))

    
    #* make folder to hold all the csv files 
    path = mkpath("OscillatorPaper/FigureGenerationScripts/3FixedCSVRawSets/$(param1)_$(param2)_$(param3)")
    i = 1

    #* make progress bar 
    # loopprogress = Progress(rangelength^3, desc ="Looping thru fixed triplets: " , color=:blue)
    for val1 in fixed_values1
        for val2 in fixed_values2
            for val3 in fixed_values3
                fixed_values = [val1, val2, val3]
                # @info fixed_values
                make_fitness_function_closure(evalfunc,prob) = make_fitness_function_with_fixed_inputs(evalfunc, prob, fixed_values, fixedtrip_idxs)
                oscillatory_points_df = run_GA(fixed_ga_problem, make_fitness_function_closure; population_size = 100, iterations = 1) 
                num_oscillatory_points = length(oscillatory_points_df.ind)

                if iszero(num_oscillatory_points)
                    results_df[i, :] = (val1, val2, val3, 0, NaN, NaN, NaN, NaN, NaN, NaN)
                else
                    average_period::Float64 = mean(oscillatory_points_df.per)
                    maximum_period::Float64 = maximum(oscillatory_points_df.per)
                    minimum_period::Float64 = minimum(oscillatory_points_df.per)
    
                    average_amplitude::Float64 = mean(oscillatory_points_df.amp)
                    maximum_amplitude::Float64 = maximum(oscillatory_points_df.amp)
                    minimum_amplitude::Float64 = minimum(oscillatory_points_df.amp)
    
                    results_df[i, :] = (val1, val2, val3, num_oscillatory_points, average_period, maximum_period, minimum_period, average_amplitude, maximum_amplitude, minimum_amplitude)
                
                    #* insert the fixed params into each ind of oscillatory_points_df
                    for ind in oscillatory_points_df.ind
                        for (j,fixedidx) in enumerate(fixedtrip_idxs)
                            if fixedidx <= length(ind)
                                insert!(ind, fixedtrip_idxs[j], fixed_values[j])
                            else
                                push!(ind, fixed_values[j])
                            end
                        end
                    end
                    # return oscillatory_points_df
                    #* split parameter values into separate columns and add initial conditions
                    split_dataframe!(oscillatory_points_df, prob)
                    CSV.write(path*"/$(round(val1; digits = 2))_$(round(val2;digits = 2))_$(round(val3; digits=2)).csv", oscillatory_points_df)
                end
                # next!(loopprogress)
                i += 1
                println(i)
            end
        end
    end
    return results_df
end

param_triplet = ["kcat1", "kcat7", "DF"]

results_df = fixed_triplet_csv_maker(param_triplet..., param_constraints, ogprob)




#! TESTING GA FUNCTIONALITY
test_gaproblem = GAProblem(param_constraints, ogprob)

test_results = run_GA(test_gaproblem; population_size = 5000, iterations = 5)


testind = test_results.trace[2].metadata["staterecord"][1].ind
testsol = solve(remake(ogprob, p=testind), saveat = 0.01, save_idxs = 1)
CostFunction(testsol)

#TODO fix the evaluation of the first generation in initial_state
#* have ruled out the in-place objective function, as well as my replacing the map with a for loop
#* doesn't seem to be an issue with the trace recording either 
#* something wrong at the time of evaluation within initial_state


eval_param_fitness(test_results.ind[2], ogprob)



plotsol(row) = plotsol(row, test_results, ogprob)

plotsol(1)


test_fitness(row) = eval_param_fitness(test_results.ind[row], ogprob)
test_fitness(500)

reogprob = remake(ogprob, p=test_results.ind[1])
testsol = solve(reogprob, saveat = 0.01, save_idxs = 1)
CostFunction(testsol)


plot(testsol)


fitind = [0.005904610505030443, 131.13364934592303, 314.4503575261949, 
        3.3658232279193547, 3.0772188176974846, 1.6497701059703715, 
        0.22809186527536418, 16.53422149313229, 0.8938599354855014, 
        3.34373724825129, 0.6270815692033677, 101.721977264349, 44216.07505733261]

fitprob = remake(ogprob, p=fitind)
fitsol = solve(fitprob, saveat = 0.01, save_idxs = 1)
CostFunction(fitsol)

