begin 
    using Plots; #theme(:juno)
    using Catalyst
    using DifferentialEquations, ModelingToolkit
    using Statistics
    using Peaks
    using FindPeaks1D
    using Evolutionary, FFTW
    using Random
    using Distributions
    using DataFrames, DataFrameMacros
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
    using Combinatorics
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

    include("../../UTILITIES/TestBenchPlotUtils.jl")

    # include("../../UTILITIES/UnitTools.jl")


    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))
end


# debuglogger = ConsoleLogger(stderr, Logging.Debug)
# infologger = ConsoleLogger(stderr, Logging.Info)
# global_logger(infologger)

tspan = (0., 2000.0)
fullrn = make_fullrn()
ogprob = ODEProblem(fullrn, [], tspan, [])



de = modelingtoolkitize(ogprob)

ogprobjac = ODEProblem(de, [], tspan, jac=true)

ogsol = solve(ogprobjac, saveat=0.1, save_idxs = [6, 9, 10, 11, 12, 15, 16])


Array(ogsol)

map(sum,ogsol.u)
# @code_warntype make_fitness_function(eval_param_fitness, ogprobjac; fitidx = 4)

# @btime solve($ogprob, saveat = 0.1, save_idxs = 1)
# @btime solve($ogprobjac, saveat = 0.1, save_idxs = 1)

# newprobjac = remake(ogprobjac, p = ogprobjac.p .* 10)
# @btime solve($newprobjac, saveat = 0.1, save_idxs = 1)



# new_u0 = ogprob.u0 .* 10
# ogprob = remake(ogprob, u0 = new_u0)
# # @benchmark solve($ogprob, saveat = 0.1, save_idxs = 1)

# @benchmark solve($ogprob, saveat=0.1, save_idxs = 1)
# ogsol = solve(ogprob, saveat=0.1, save_idxs = 1)
# peakidxs, props = findpeaks1d(ogsol[1,:]; height = 0.1)
# plot(ogsol)
# testfunc(ogprob) = solve(ogprob, saveat=0.1, save_idxs = 1)

# using Cthulhu
# using ProfileView
# descend_code_warntype(testfunc, (ODEProblem,))
# @code_warntype solve(ogprob, saveat=0.1, save_idxs = 1)
# # plot(ogsol)
# @code_warntype solve_for_fitness_peramp(ogprob)

# @benchmark solve_for_fitness_peramp($ogprob)

# @code_warntype CostFunction(ogsol)
# @benchmark CostFunction($ogsol)
# @benchmark getPerAmp($ogsol)



#* Optimization of parameters to produce data for CSV
param_constraints = define_parameter_constraints(ogprobjac; karange = (1e-3, 1e2), kbrange = (1e-3, 1e3), kcatrange = (1e-3, 1e3), dfrange = (1e2, 2e4))
ic_constraints = define_initialcondition_constraints(ogprobjac; Lrange = (1e-2, 1e2), Krange = (1e-2, 1e2), Prange = (1e-2, 1e2), Arange = (1e-2, 1e2))

allconstraints = AllConstraints(param_constraints, ic_constraints)


#* Modification to make_fitness_function_with_fixed_inputs function
function make_fitness_function_with_fixed_inputs_bothparamsIC(evalfunc::Function, prob::ODEProblem, fixed_input_triplet::Vector{Float64}, triplet_idxs::Tuple{Int, Int, Int}, fixedDF=1000.)
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

        new_param_input = new_input[1:12]
        new_ic_input = new_input[13:end]

        push!(new_param_input, fixedDF)

        return evalfunc(new_param_input, new_ic_input, prob)
    end
    return fitness_function
end




#* Modification to fixed_triplet_csv_maker function
function fixed_triplet_csv_maker(param1::String, param2::String, param3::String, constraints::AllConstraints, prob::ODEProblem; rangelength = 4, fixedDF = 1000.) 
    variable_constraints = deepcopy(constraints)

    #* get the fixed parameters from the variable constraints
    fixedval_constraints = [x for x in variable_constraints.ranges if x.name == param1 || x.name == param2 || x.name == param3]

    #* remove the fixed parameters from the variable constraints so only the variable parameters are used in the GA
    filter!(x -> x.name != param1 && x.name != param2 && x.name != param3, variable_constraints.ranges)

    #* filter out DF because it will be fixed
    filter!(x -> x.name != "DF", variable_constraints.ranges)


    fixed_ga_problem = GAProblem(variable_constraints, prob)

    #* get the indices of the fixed parameters in the input vector so that they can be inserted after GA
    fixedval_idxs = find_indices(param1, param2, param3, constraints.ranges) 

    
    fixed_valrange1 = logrange(fixedval_constraints[1].min, fixedval_constraints[1].max, rangelength)
    
    fixed_valrange2 = logrange(fixedval_constraints[2].min, fixedval_constraints[2].max, rangelength)
    
    fixed_valrange3 = logrange(fixedval_constraints[3].min, fixedval_constraints[3].max, rangelength)

    num_rows = rangelength^3

    vals1 = Vector{Float64}(undef, num_rows)
    vals2 = Vector{Float64}(undef, num_rows)
    vals3 = Vector{Float64}(undef, num_rows)
    num_oscillatory_points_array = Vector{Int}(undef, num_rows)
    average_periods = Vector{Float64}(undef, num_rows)
    maximum_periods = Vector{Float64}(undef, num_rows)
    minimum_periods = Vector{Float64}(undef, num_rows)
    average_amplitudes = Vector{Float64}(undef, num_rows)
    maximum_amplitudes = Vector{Float64}(undef, num_rows)
    minimum_amplitudes = Vector{Float64}(undef, num_rows)

    # results_df = DataFrame(param1 => Vector{Float64}(undef, rangelength^3), param2 => Vector{Float64}(undef, rangelength^3), param3 => Vector{Float64}(undef, rangelength^3), "num_oscillatory_points" => Vector{Int}(undef, rangelength^3), 
    #                     "average_period" => Vector{Float64}(undef, rangelength^3), "maximum_period"=>Vector{Float64}(undef, rangelength^3), "minimum_period"=>Vector{Float64}(undef, rangelength^3),
    #                     "average_amplitude" => Vector{Float64}(undef, rangelength^3), "maximum_amplitude"=>Vector{Float64}(undef, rangelength^3), "minimum_amplitude"=>Vector{Float64}(undef, rangelength^3))

    
    #* make folder to hold all the csv files 
    path = mkpath("OscillatorPaper/FigureGenerationScripts/3FixedParams+ICsRawSets/$(param1)_$(param2)_$(param3)/DF=$(round(fixedDF))")
    i = 1

    #* make progress bar 
    loopprogress = Progress(num_rows, desc ="Looping thru fixed triplets: " , color=:blue)
    for val1 in fixed_valrange1
        for val2 in fixed_valrange2
            for val3 in fixed_valrange3
                fixed_values = [val1, val2, val3]
                @info fixed_values

                #* make fitness function closure with fixed inputs
                make_fitness_function_closure(evalfunc,prob) = make_fitness_function_with_fixed_inputs_bothparamsIC(evalfunc, prob, fixed_values, fixedval_idxs, fixedDF)

                Random.seed!(1234)
                oscillatory_points_results = run_GA(fixed_ga_problem, make_fitness_function_closure; population_size = 100000, iterations = 5) 
                num_oscillatory_points = length(oscillatory_points_results.population)

                #* if there are no oscillatory points, save the results to the results_df and continue
                if iszero(num_oscillatory_points)
                    vals1[i] = val1
                    vals2[i] = val2
                    vals3[i] = val3
                    num_oscillatory_points_array[i] = 0
                    average_periods[i] = NaN
                    maximum_periods[i] = NaN
                    minimum_periods[i] = NaN
                    average_amplitudes[i] = NaN
                    maximum_amplitudes[i] = NaN
                    minimum_amplitudes[i] = NaN
                else
                    average_periods[i]::Float64 = mean(oscillatory_points_results.periods)
                    maximum_periods[i]::Float64 = maximum(oscillatory_points_results.periods; init=0.0)
                    minimum_periods[i]::Float64 = minimum(oscillatory_points_results.periods; init=0.0)

                    average_amplitudes[i]::Float64 = mean(oscillatory_points_results.amplitudes)
                    maximum_amplitudes[i]::Float64 = maximum(oscillatory_points_results.amplitudes; init=0.0)
                    minimum_amplitudes[i]::Float64 = minimum(oscillatory_points_results.amplitudes; init=0.0)
                    
                    #* save the results to the results_df
                    vals1[i] = val1
                    vals2[i] = val2
                    vals3[i] = val3
                    num_oscillatory_points_array[i] = num_oscillatory_points
                
                    #* make dataframe from oscillatory_points_results
                    oscillatory_points_df = make_df(oscillatory_points_results)
                    
                    #* insert the fixed params into each ind of oscillatory_points_df
                    for ind in oscillatory_points_df.ind
                        for (j,fixedidx) in enumerate(fixedval_idxs)
                            if fixedidx <= length(ind)
                                insert!(ind, fixedval_idxs[j], fixed_values[j])
                            else
                                push!(ind, fixed_values[j])
                            end
                        end
                        insert!(ind, 13, fixedDF)
                    end
                    #* split parameter values into separate columns and add initial conditions
                    split_dataframe!(oscillatory_points_df, prob)

                    #* rewrite the L, K, P, A columns with the initial conditions
                    CSV.write(path*"/$(round(val1; digits = 2))_$(round(val2;digits = 2))_$(round(val3; digits=2)).csv", oscillatory_points_df)
                end
                next!(loopprogress)
                i += 1
            end
        end
    end
    return results_df
end

param_triplet = ["kcat1", "kcat7", "DF"]

results_df = fixed_triplet_csv_maker(param_triplet..., allconstraints, ogprobjac)


# CSV.write("OscillatorPaper/FigureGenerationScripts/3FixedResultsCSVs/fixed_triplet_results-$(param_triplet[1]*param_triplet[2]*param_triplet[3]).csv", results_df)


#< loop through combinations of parameters and run the function of each triplet

function run_all_triplets(constraints::AllConstraints, prob::ODEProblem; startinput::String="ka1", rangelength = 4)
    names = ["ka1", "kb1", "kcat1", "ka2", "kb2", "ka3", "kb3", "ka4", "kb4", "ka7", "kb7", "kcat7", "DF", "L", "K", "P", "A"]
    triplets = collect(combinations(names, 3))
    start_idx = findfirst(x -> x[1] == startinput, triplets)

    mainpath = mkpath("OscillatorPaper/FigureGenerationScripts/3FixedParams+ICsFixedDF_CSVs")

    #* loops through these fixed DF values for each triplet
    df_range = [100., 1000., 10000.]

    for triplet in triplets[start_idx:end]
        @info triplet
        tripletpath = mkpath(mainpath*"/$(triplet[1])_$(triplet[2])_$(triplet[3])")
        for df in df_range
            @info df
            results_df = fixed_triplet_csv_maker(triplet..., constraints, prob; rangelength=rangelength, fixedDF=df)
            CSV.write(tripletpath*"/DF=$(df).csv", results_df)
        end
    end
end

# @code_warntype run_all_triplets(allconstraints, ogprobjac)

run_all_triplets(allconstraints, ogprobjac; startinput="ka1", rangelength=4)







#< Bugs to fix
#* 1. The fitness function is not working properly. Optimize, reduce fine tuning, and test it for edge cases against expected results
#! * 2. Fix the "FFTW can't make plan" error 
#! * 3. Logscale projection isn't working, fix it. 
#! 4. Save the optimized parameters, not just the evaluation values 
#* 5. Run GA through debugger to see the sequence of selection, recombination
#! 6. Print out the initial conditions in the CSV 
#* 7. Make test suite for the fitness function. Orthogonal tests will run through all optimized solutions and classify them as correct, false negative, false positive.
#* 8. Save all variables for selected solutions. IDK Maggie mentioned it 
#* Run longer sims 
#* Duplicates in raw 
#* Validate solutions from raw df, then see it written to csv correctly  
#* Idea for FFT: downsampling, or smoothing
#* Look into the recombine, duplicates. More duplicates with each generation 
#* Rerun for twice as long if 2 index check is failed in getPerAmp
#* Fix DF in the hundreds
#* optimize initial conditions too, and fix DF 


