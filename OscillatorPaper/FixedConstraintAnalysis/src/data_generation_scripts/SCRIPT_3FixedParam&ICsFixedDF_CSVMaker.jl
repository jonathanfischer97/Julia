begin 
    # using Plots; #theme(:juno)
    using Catalyst
    using OrdinaryDiffEq, ModelingToolkit
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

    # include("../../UTILITIES/TestBenchPlotUtils.jl")

    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))

    numthreads = Threads.nthreads()
    println("Threads detected: $numthreads")
    numcores = min(1,numthreads÷2)
    BLAS.set_num_threads(numcores)
    FFTW.set_num_threads(numcores)
end





#* Modification to fixed_triplet_csv_maker function
function fixed_triplet_csv_maker(constraints::AllConstraints, ode_prob::ODEProblem; path, rangelength::Int = 4, popsize::Int = 20000) 

    fixed_constraintranges = get_fixed_constraintranges(constraints)

    fixed_names = [conrange.name for conrange in fixed_constraintranges]
    # @info "Fixed inputs: $fixed_names"
    
    fixed_valrange1, fixed_valrange2, fixed_valrange3 = [logrange(constraintrange.min, constraintrange.max, rangelength) for constraintrange in fixed_constraintranges]
    

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

    #* get the fixed DF value
    fixedDF = constraints.DF.fixed_value

    #* make folder to hold all the csv files 
    DFpath = mkpath(path*"/DF=$(round(fixedDF))")
    i = 1

    #* make progress bar 
    # loopprogress = Progress(num_rows, desc ="Looping thru fixed triplets: " , color=:blue)

    initial_population = generate_empty_population(constraints, popsize)

    ga_problem = GAProblem(constraints = constraints, ode_problem = ode_prob)

    for val1 in fixed_valrange1
        for val2 in fixed_valrange2
            for val3 in fixed_valrange3
                #* set the fixed values in the constraints
                set_fixed_values!(fixed_constraintranges, val1, val2, val3)

                # @info get_fixed_indices(fixed_constraints)

                Random.seed!(1234)
            
                generate_population!(initial_population, constraints)
                # @info length(initial_population[1])
            
                oscillatory_points_results = run_GA(ga_problem, initial_population; iterations = 5, show_trace = false)
            
                num_oscillatory_points = length(oscillatory_points_results.fitvals)
                println("Number of oscillatory points: $num_oscillatory_points")


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
                    maximum_periods[i]::Float64 = maximum(oscillatory_points_results.periods)
                    minimum_periods[i]::Float64 = minimum(oscillatory_points_results.periods)

                    average_amplitudes[i]::Float64 = mean(oscillatory_points_results.amplitudes)
                    maximum_amplitudes[i]::Float64 = maximum(oscillatory_points_results.amplitudes)
                    minimum_amplitudes[i]::Float64 = minimum(oscillatory_points_results.amplitudes)
                    
                    #* save the results to the results_df
                    vals1[i] = val1
                    vals2[i] = val2
                    vals3[i] = val3
                    num_oscillatory_points_array[i] = num_oscillatory_points
                    

                    #* rewrite the L, K, P, A columns with the initial conditions
                    save_to_csv(oscillatory_points_results, constraints, DFpath*"/$(fixed_names[1])=$(round(val1; digits = 2))_$(fixed_names[2])=$(round(val2;digits = 2))_$(fixed_names[3])=$(round(val3; digits=2)).csv")
                end
                # next!(loopprogress)
                i += 1
            end
        end
    end
    
    return DataFrame([fixed_names[1] => vals1, fixed_names[2] => vals2, fixed_names[3] => vals3, :DF => fixedDF, :num_oscillatory_points => num_oscillatory_points_array, 
                        :average_period => average_periods, :maximum_period => maximum_periods, :minimum_period => minimum_periods,
                        :average_amplitude => average_amplitudes, :maximum_amplitude => maximum_amplitudes, :minimum_amplitude => minimum_amplitudes])
end


#< loop through combinations of parameters/ics and run the function of each combination
# function run_all_triplets(constraints::AllConstraints, prob::ODEProblem; rootpath::String, startinput::Symbol=:ka1, rangelength = 4, DFrange = [100., 1000., 10000.], popsize=20000)
#     names = [:ka1, :kb1, :kcat1, :ka2, :kb2, :ka3, :kb3, :ka4, :kb4, :ka7, :kb7, :kcat7, :L, :K, :P, :A]
#     triplets = collect(combinations(names, 3))
#     start_idx = findfirst(x -> x[1] == startinput, triplets)

#     #* mark DF as fixed in the constraints. Fixed value won't be assigned until loop
#     df_constraintrange = constraints.DF
#     df_constraintrange.isfixed = true


#     for triplet in triplets[start_idx:end]
#         tripletpath = mkpath(rootpath*"/$(triplet[1])_$(triplet[2])_$(triplet[3])")
#         summarypath = mkpath(tripletpath*"/SummaryResults")
#         rawpath = mkpath(tripletpath*"/RawData")
#         for df in DFrange
#             @info "$triplet with DF = $df"
#             df_constraintrange.fixed_value = df
#             set_fixed_constraints!(constraints, triplet)
#             results_df = fixed_triplet_csv_maker(constraints, prob; rangelength=rangelength, fixedDF=df, popsize=popsize, path=rawpath)
#             CSV.write(summarypath*"/Summary_DF=$(round(df)).csv", results_df)
#         end
#         unset_fixed_constraints!(constraints, triplet)
#     end
# end


# function run_MAIN(rangelength=4, popsize=20000)
#     tspan = (0., 2000.0)
#     fullrn = make_fullrn()
#     ogprob = ODEProblem(fullrn, [], tspan, [])

#     de = modelingtoolkitize(ogprob)

#     ogprobjac = ODEProblem(de, [], tspan, jac=true)

#     #* Optimization of parameters to produce data for CSV
#     param_constraints = ParameterConstraints(; karange = (1e-3, 1e2), kbrange = (1e-3, 1e3), kcatrange = (1e-3, 1e3), dfrange = (1e2, 2e4))
#     ic_constraints = InitialConditionConstraints(; Lrange = (1e-1, 1e2), Krange = (1e-2, 1e2), Prange = (1e-2, 1e2), Arange = (1e-1, 1e2))

#     allconstraints = AllConstraints(param_constraints, ic_constraints)

#     rootpath = mkpath("./ROCKFISH_DATA/3Fixed/$PopSize_$popsize")

#     run_all_triplets(allconstraints, ogprobjac; rootpath=rootpath, startinput=:ka1, rangelength=rangelength, popsize=popsize)
# end


# if isempty(ARGS)
#     run_MAIN()
# else
#     run_MAIN(parse(Int, ARGS[1]), parse(Int, ARGS[2]))
# end


function run_all_triplets(constraints::AllConstraints, prob::ODEProblem; start_idx::Int, end_idx::Int, rootpath::String, rangelength=4, DFrange=[100., 1000., 10000.], popsize=20000)
    #* Assert that no other constraints are fixed or assigned fixed values except for DF
    @assert all([constraintrange.isfixed == false && isnan(constraintrange.fixed_value) for constraintrange in constraints if constraintrange.name != :DF])
   
    names = [:ka1, :kb1, :kcat1, :ka2, :kb2, :ka3, :kb3, :ka4, :kb4, :ka7, :kb7, :kcat7, :L, :K, :P, :A]
    triplets = collect(combinations(names, 3))

    df_constraintrange = constraints.DF
    df_constraintrange.isfixed = true

    proglength = length(triplets[start_idx:min(end_idx, length(triplets))])
    tripletprogress = Progress(proglength, desc ="Looping thru triplets: ")


    for triplet in triplets[start_idx:min(end_idx, length(triplets))]
        tripletpath = mkpath(rootpath*"/$(triplet[1])_$(triplet[2])_$(triplet[3])")
        summarypath = mkpath(tripletpath*"/SummaryResults")
        rawpath = mkpath(tripletpath*"/RawData")
        for df in DFrange
            @info "$triplet with DF = $df"
            df_constraintrange.fixed_value = df
            set_fixed_constraints!(constraints, triplet)
            results_df = fixed_triplet_csv_maker(constraints, prob; rangelength=rangelength, popsize=popsize, path=rawpath)
            CSV.write(summarypath*"/Summary_DF=$(round(df)).csv", results_df)
        end
        unset_fixed_constraints!(constraints, triplet)
        next!(tripletprogress)
    end
end

function run_MAIN(rangelength=4, popsize=20000; start_idx=1, end_idx=560)

    ogprobjac = make_ODE_problem()

    allconstraints = AllConstraints()

    rootpath = mkpath("./ROCKFISH_DATA/3Fixed/PopSize_$popsize")

    run_all_triplets(allconstraints, ogprobjac; start_idx=start_idx, end_idx=end_idx, rootpath=rootpath, rangelength=rangelength, popsize=popsize)
end


start_idx = parse(Int, ENV["START_IDX"])
end_idx = parse(Int, ENV["END_IDX"])

run_MAIN(parse(Int, ARGS[1]), parse(Int, ARGS[2]); start_idx, end_idx)






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


