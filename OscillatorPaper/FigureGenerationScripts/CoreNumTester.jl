begin 
    using Catalyst
    using OrdinaryDiffEq, ModelingToolkit
    using Statistics

    using Evolutionary, FFTW
    using Random
    using Distributions
    using DataFrames#, DataFrameMacros
    using CSV

    using LinearAlgebra
    using ProgressMeter

    # using Combinatorics

    include("../../UTILITIES/EvolutionaryOverloads.jl")

    # import the Catalyst model "fullrn"
    include("../../UTILITIES/ReactionNetwork.jl")

    # import the cost function and other evaluation functions
    include("../../UTILITIES/EvaluationFunctions.jl")
    # using .EvaluationFunctions

    # import the genetic algorithm and associated functions
    include("../../UTILITIES/GA_functions.jl")

    # include("../../UTILITIES/TestBenchPlotUtils.jl")

    # const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))

    numthreads = Threads.nthreads()
    numcores = max(1,numthreads÷2)
    BLAS.set_num_threads(1)
    FFTW.set_num_threads(1)
end

println("Number of threads: ", numthreads)
println("Number of cores: ", numcores)


"""Function loops through 4D grid of different initial conditions, letting all parameters be freely optimized, and saves the results to a csv file"""
function fixed_quadruplet_ic_searcher(paramconstraints::ParameterConstraints, icconstraints::InitialConditionConstraints, prob::ODEProblem; rangelength::Int = 4, fixedDF::Float64=1000., popsize::Int=20000)
    #* get the ranges of the fixed initial conditions that will be looped through
    icranges = [logrange(constraints.min, constraints.max, rangelength) for constraints in icconstraints]

    icnames = [constraints.name for constraints in icconstraints]

    #* construct the AllConstraints object
    allconstraints = AllConstraints(paramconstraints, icconstraints)

    #* mark the fixed constraints, the fixed values will be assigned later in the loop
    set_fixed_constraints!(allconstraints, icnames)

    #* unpack the fixed constraintranges in a vector
    fixed_constraintranges = get_fixed_constraintranges(allconstraints)

    #* set the fixed DF value
    allconstraints.DF.isfixed = true
    allconstraints.DF.fixed_value = fixedDF

    #* number of rows in the dataframe 
    num_rows = rangelength^length(icnames)

    #* initialize vectors to hold the results, will later be made into dataframe and written to CSV
    icvals1 = Vector{Float64}(undef, num_rows)
    icvals2 = Vector{Float64}(undef, num_rows)
    icvals3 = Vector{Float64}(undef, num_rows)
    icvals4 = Vector{Float64}(undef, num_rows)
    num_oscillatory_points_array = Vector{Int}(undef, num_rows)
    average_periods = Vector{Float64}(undef, num_rows)
    maximum_periods = Vector{Float64}(undef, num_rows)
    minimum_periods = Vector{Float64}(undef, num_rows)
    average_amplitudes = Vector{Float64}(undef, num_rows)
    maximum_amplitudes = Vector{Float64}(undef, num_rows)
    minimum_amplitudes = Vector{Float64}(undef, num_rows)


    #* initialize counter
    i = 1

    #* make progress bar 
    # loopprogress = Progress(num_rows, desc ="Looping thru fixed ICs: " , color=:red)

    #* make path for the raw data for this particular DF value
    # DFpath = mkpath(path*"/DF=$(round(fixedDF))")

    #* initialize the population, where the length of each individual is the number of constraints minus the fixed ones
    initial_population = generate_empty_population(allconstraints, popsize)
 
    #* make the GA problem
    ga_problem = GAProblem(constraints = allconstraints, ode_problem = prob)

    #* loop through each ic range and run the GA on each set of initial conditions after remaking the problem with them
    for icval1 in icranges[1]
        for icval2 in icranges[2]
            for icval3 in icranges[3]
                for icval4 in icranges[4]

                    #* set the fixed values in the fixed_constraintranges vector
                    set_fixed_values!(fixed_constraintranges, icval1, icval2, icval3, icval4)

                    #* set seed for reproducibility
                    Random.seed!(1234)

                    #* run the GA on the new problem
                    generate_population!(initial_population, allconstraints)
                    oscillatory_points_results = run_GA(ga_problem, initial_population; iterations = 5, show_trace = false)

                    #* get the number of oscillatory points
                    num_oscillatory_points = length(oscillatory_points_results.fitvals)
                    # println("Number of oscillatory points for L=$icval1, K=$icval2, P=$icval3, A=$icval4: $num_oscillatory_points")

                    #* if there are no oscillatory points, save the results to the results_df and continue
                    if iszero(num_oscillatory_points)
                        icvals1[i] = icval1
                        icvals2[i] = icval2
                        icvals3[i] = icval3
                        icvals4[i] = icval4
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
                        icvals1[i] = icval1
                        icvals2[i] = icval2
                        icvals3[i] = icval3
                        icvals4[i] = icval4
                        num_oscillatory_points_array[i] = num_oscillatory_points

                        #* save the dataframe of this particular fixed value combination to a csv file, identified by the fixed initial conditions
                        # csv_filestring = DFpath*"/L=$(round(icval1; digits = 2))_K=$(round(icval2;digits = 2))_P=$(round(icval3; digits=2))_A=$(round(icval4; digits=2)).csv"
                        # println("Saving results to $csv_filestring")
                        # CSV.write(csv_filestring, oscillatory_points_df)
                        # save_to_csv(oscillatory_points_results, allconstraints, csv_filestring)
                    end

                    # next!(loopprogress)
                    i += 1
                end
            end
        end
    end

    #* make the results dataframe, which holds all the summary statistics for the whole set of fixed initial conditions
    return DataFrame(icnames[1] => icvals1, icnames[2] => icvals2, icnames[3] => icvals3, icnames[4] => icvals4,
                            :num_oscillatory_points => num_oscillatory_points_array, 
                            :average_period => average_periods, :maximum_period => maximum_periods, :minimum_period => minimum_periods,
                            :average_amplitude => average_amplitudes, :maximum_amplitude => maximum_amplitudes, :minimum_amplitude => minimum_amplitudes)          
end



"""Loops through each fixed value of DF and runs the fixed_quadruplet_ic_searcher function"""
# function loop_4fixedICs_thru_DFvals(paramconstraints::ParameterConstraints, icconstraints::InitialConditionConstraints, prob::ODEProblem; rangelength::Int = 4, DFrange = [100.,1000.,10000.], popsize::Int=20000)
#     rootpath = mkpath("./ROCKFISH_DATA/4Fixed/PopSize_$popsize")
#     summarypath = mkpath(rootpath*"/SummaryResults")
#     rawpath = mkpath(rootpath*"/4FixedICRawSets")
#     for DF in DFrange
#         println("Running DF = $DF")
#         results_df = fixed_quadruplet_ic_searcher(paramconstraints, icconstraints, prob; rangelength=rangelength, fixedDF=DF, popsize=popsize, path=rawpath)
#         CSV.write(summarypath*"/Summary_DF=$(round(DF)).csv", results_df)
#     end
# end



function run_4fixedIC(;rangelength=3, popsize=10000, fixedDF=1000.)

    ogprobjac = make_ODE_problem()

    param_constraints = ParameterConstraints(;karange = (1e-3, 1e2), kbrange = (1e-3, 1e3), kcatrange = (1e-3, 1e3), dfrange = (1e2, 2e4))
    ic_constraints = InitialConditionConstraints(; Lrange = (1e-1, 1e2), Krange = (1e-2, 1e2), Prange = (1e-2, 1e2), Arange = (1e-1, 1e2))

    #* make paths
    # rootpath = mkpath("./ROCKFISH_DATA/4Fixed/PopSize_$popsize")
    # summarypath = mkpath(rootpath*"/SummaryResults")
    # rawpath = mkpath(rootpath*"/4FixedICRawSets")

    results_df = fixed_quadruplet_ic_searcher(param_constraints, ic_constraints, ogprobjac; rangelength=rangelength, fixedDF=fixedDF, popsize=popsize)
    # CSV.write(summarypath*"/Summary_DF=$(round(fixedDF)).csv", results_df)
    return nothing
    # loop_4fixedICs_thru_DFvals(param_constraints, ic_constraints, ogprobjac; rangelength=rangelength, DFrange = [100.,1000.,10000.], popsize=popsize)
end

# @time run_4fixedIC(5, 10000)

# println(time())

time_results = @timed run_4fixedIC(;popsize = parse(Int, ARGS[1]))

# time_df = DataFrame(numcores = numcores, time = time_results.t, bytes = time_results.bytes, gctime = time_results.gctime)
# CSV.write("OscillatorPaper/FigureGenerationScripts/TimerResults", time_df; append=true)


println("Population size: ", parse(Int, ARGS[1]))
println("Time to run: ", time_results.time)
println("Allocations: ", time_results.bytes)
println("GC time: ", time_results.gctime, "\n")




# using Plots 

# corelist = [1,5,9,13,17,21,25,29,33,37,41,45]
# timelist = [49998.015669374, 11992.902370869, 7919.690304564, 6804.007481667, 4605.985745292, 4661.490750618, 3923.057727081, 3676.451899478, 3620.190313874, 3443.095502602, 3310.522324744, 3341.668692622]
# gclist = [2056.171880943, 1677.28337679, 1932.153616315, 2223.68263324, 1459.319987119, 1911.251829509, 1583.697396817, 1571.539846504, 1658.581128171, 1583.471377401, 1513.257647106, 1570.993172414]

# plot(corelist, timelist, label="Time", xlabel="Number of cores", ylabel="Time (s)", title="Time to run 4FixedICSearcher.jl", legend=:bottomright)
# scatter!(corelist, timelist, label="")
# plot!(corelist, gclist, label="GC time")

"
Number of cores: 1
Population size: 10000
Time to run: 49998.015669374
Allocations: 15845712497323
GC time: 2056.171880943

Number of cores: 5
Population size: 10000
Time to run: 11992.902370869
Allocations: 15845692774267
GC time: 1677.28337679

Number of cores: 9
Population size: 10000
Time to run: 7919.690304564
Allocations: 15845666705931
GC time: 1932.153616315

Number of cores: 13
Population size: 10000
Time to run: 6804.007481667
Allocations: 15845590577947
GC time: 2223.68263324

Number of cores: 17
Population size: 10000
Time to run: 4605.985745292
Allocations: 15845545320075
GC time: 1459.319987119

Number of cores: 21
Population size: 10000
Time to run: 4661.490750618
Allocations: 15845545762443
GC time: 1911.251829509

Number of cores: 25
Population size: 10000
Time to run: 3923.057727081
Allocations: 15845614097227
GC time: 1583.697396817

Number of cores: 29
Population size: 10000
Time to run: 3676.451899478
Allocations: 15845530792251
GC time: 1571.539846504

Number of cores: 33
Population size: 10000
Time to run: 3620.190313874
Allocations: 15845543419707
GC time: 1658.581128171

Number of cores: 37
Population size: 10000
Time to run: 3443.095502602
Allocations: 15845462953019
GC time: 1583.471377401

Number of cores: 41
Population size: 10000
Time to run: 3310.522324744
Allocations: 15845511154699
GC time: 1513.257647106

Number of cores: 45
Population size: 10000
Time to run: 3341.668692622
Allocations: 15845516164411
GC time: 1570.993172414
"