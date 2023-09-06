begin 
    using Plots; #theme(:juno)
    using Catalyst
    using DifferentialEquations, ModelingToolkit
    using Statistics
    using Evolutionary, FFTW
    using Random
    using Distributions
    using DataFrames#, DataFrameMacros
    using CSV
    # using Setfield

    # using StaticArrays
    # using BenchmarkTools, Profile
    using ProgressMeter

    using LinearAlgebra

    # using ColorSchemes, Plots.PlotMeasures

    using Combinatorics

    # default(lw = 2, size = (1000, 600), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)

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


    numthreads = Threads.nthreads()
    @info "Threads detected: $numthreads"
    numcores = numthreads÷2
    BLAS.set_num_threads(numcores)
    FFTW.set_num_threads(numcores)
end





#* Modification to fixed_triplet_csv_maker function
function fixed_triplet_csv_maker(constraints::AllConstraints, ode_prob::ODEProblem, fixed_names; rangelength = 4, fixedDF = 1000., popsize = 20000) 

    fixedval_idxs = (find_field_index(fixed_name, constraints) for fixed_name in fixed_names)
    
    fixed_valrange1, fixed_valrange2, fixed_valrange3 = (logrange(constraints[fixedval_idx].min, constraints[fixedval_idx].max, rangelength) for fixedval_idx in fixedval_idxs)
    

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

    
    #* make folder to hold all the csv files 
    path = mkpath("./ROCKFISH_DATA/3Fixed/3FixedParams+ICsRawSets/$(fixed_names[1])_$(fixed_names[2])_$(fixed_names[3])/DF=$(round(fixedDF))")
    i = 1

    #* make progress bar 
    loopprogress = Progress(num_rows, desc ="Looping thru fixed triplets: " , color=:blue)
    for val1 in fixed_valrange1
        for val2 in fixed_valrange2
            for val3 in fixed_valrange3
                fixed_inputs = [fixed_name => val for (fixed_name, val) in zip(fixed_names, [val1, val2, val3])]
                push!(fixed_inputs, "DF" => fixedDF)

                fixed_constraints = deepcopy(constraints)
                
                set_fixed_constraints!(fixed_constraints; fixed_inputs)

                @info get_fixed_indices(fixed_constraints)

                ga_problem = GAProblem(fixed_constraints, ode_prob)


                Random.seed!(1234)
            
                initial_population = generate_population(fixed_constraints, popsize)
            
                oscillatory_points_results = run_GA(ga_problem, initial_population; iterations = 5)
            
                num_oscillatory_points = length(oscillatory_points_results.fitvals)


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
                    oscillatory_points_df = make_ga_dataframe(oscillatory_points_results, fixed_constraints)
                    

                    #* rewrite the L, K, P, A columns with the initial conditions
                    CSV.write(path*"/$(round(val1; digits = 2))_$(round(val2;digits = 2))_$(round(val3; digits=2)).csv", oscillatory_points_df)
                end
                next!(loopprogress)
                i += 1
            end
        end
    end
    
    results_df = DataFrame(fixed_names[1] => vals1, fixed_names[2] => vals2, fixed_names[3] => vals3, "DF" => fixedDF, "num_oscillatory_points" => num_oscillatory_points_array, 
                        "average_period" => average_periods, "maximum_period"=>maximum_periods, "minimum_period"=>minimum_periods,
                        "average_amplitude" => average_amplitudes, "maximum_amplitude"=>maximum_amplitudes, "minimum_amplitude"=>minimum_amplitudes)
    return results_df
end



#< loop through combinations of parameters and run the function of each triplet

function run_all_triplets(constraints::AllConstraints, prob::ODEProblem; startinput::String="ka1", rangelength = 4, DFrange = [100., 1000., 10000.])
    names = ["ka1", "kb1", "kcat1", "ka2", "kb2", "ka3", "kb3", "ka4", "kb4", "ka7", "kb7", "kcat7", "DF", "L", "K", "P", "A"]
    triplets = collect(combinations(names, 3))
    start_idx = findfirst(x -> x[1] == startinput, triplets)

    mainpath = mkpath("./ROCKFISH_DATA/3Fixed/3FixedParams+ICsFixedDF_CSVs")

    for triplet in triplets[start_idx:end]
        @info triplet
        tripletpath = mkpath(mainpath*"/$(triplet[1])_$(triplet[2])_$(triplet[3])")
        for df in DFrange
            @info df
            results_df = fixed_triplet_csv_maker(constraints, prob, triplet; rangelength=rangelength, fixedDF=df)
            CSV.write(tripletpath*"/DF=$(df).csv", results_df)
        end
    end
end


function run_MAIN()
    tspan = (0., 2000.0)
    fullrn = make_fullrn()
    ogprob = ODEProblem(fullrn, [], tspan, [])

    de = modelingtoolkitize(ogprob)

    ogprobjac = ODEProblem(de, [], tspan, jac=true)

    #* Optimization of parameters to produce data for CSV
    param_constraints = ParameterConstraints(; karange = (1e-3, 1e2), kbrange = (1e-3, 1e3), kcatrange = (1e-3, 1e3), dfrange = (1e2, 2e4))
    ic_constraints = InitialConditionConstraints(; Lrange = (1e-1, 1e2), Krange = (1e-2, 1e2), Prange = (1e-2, 1e2), Arange = (1e-1, 1e2))

    allconstraints = AllConstraints(param_constraints, ic_constraints)

    run_all_triplets(allconstraints, ogprobjac; startinput="ka1", rangelength=4)
end


run_MAIN()







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


