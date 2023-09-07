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

    # using StaticArrays
    # using Cthulhu
    # using JET
    # using MultivariateStats, UMAP, TSne, StatsPlots
    # using GlobalSensitivity, QuasiMonteCarlo
    using LinearAlgebra
    using ProgressMeter
    # using BifurcationKit, Setfield, ForwardDiff, Parameters; const BK = BifurcationKit
    # using OrderedCollections
    using Combinatorics
    # using LazySets, Polyhedra
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


    # const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))

    numthreads = Threads.nthreads()
    numcores = numthreads÷2
    BLAS.set_num_threads(numcores)
    FFTW.set_num_threads(numcores)
end









#* Function loops through 4D grid of different initial conditions, letting all parameters be freely optimized, and saves the results to a csv file
function fixed_quadruplet_ic_searcher(paramconstraints::ParameterConstraints, icconstraints::InitialConditionConstraints, prob::ODEProblem; rangelength::Int = 4, fixedDF::Float64=1000., popsize::Int=20000)
    #* get the ranges of the initial conditions
    icranges = [logrange(constraints.min, constraints.max, rangelength) for constraints in icconstraints]

    icnames = [constraints.name for constraints in icconstraints]

    allconstraints = AllConstraints(paramconstraints, icconstraints)

    set_fixed_constraints!(allconstraints, icnames)
    allconstraints.DF.isfixed = true
    allconstraints.DF.fixed_value = fixedDF

    num_rows = rangelength^length(icnames)


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


    i = 1

    #* make progress bar 
    loopprogress = Progress(num_rows, desc ="Looping thru fixed ICs: " , color=:red)

    mainrawpath = mkpath("./ROCKFISH_DATA/4Fixed/4FixedICRawSets")

    initial_population = generate_empty_population(allconstraints, popsize)

    ga_problem = GAProblem(allconstraints, prob)

    #* loop through each ic range and run the GA on each set of initial conditions after remaking the problem with them
    for icval1 in icranges[1]
        for icval2 in icranges[2]
            for icval3 in icranges[3]
                for icval4 in icranges[4]
                    icvals = [icval1, icval2, icval3, icval4]
                    # @info icvals

                    set_fixed_values!(icconstraints, icval1, icval2, icval3, icval4)
                    
                    #* make new GA problem with new initial conditions
                    ga_problem.fitness_function = make_fitness_function(allconstraints, prob)

                    #* set seed for reproducibility
                    Random.seed!(1234)

                    #* run the GA on the new problem
                    generate_population!(initial_population, constraints)
                    oscillatory_points_results = run_GA(ga_problem, initial_population; iterations = 5)

                    #* get the number of oscillatory points
                    num_oscillatory_points = length(oscillatory_points_results.fitvals)

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
                        
                        #* make dataframe from oscillatory_points_results
                        oscillatory_points_df = make_ga_dataframe(oscillatory_points_results, allconstraints)


                        innerrawpath = mkpath(mainrawpath*"/$(round(icval1; digits = 2))_$(round(icval2;digits = 2))_$(round(icval3; digits=2))_$(round(icval4; digits=2))")

                        CSV.write(innerrawpath*"/Raw_DF=$(round(fixedDF)).csv", oscillatory_points_df)
                    end
                    next!(loopprogress)
                    i += 1
                end
            end
        end
    end
    results_df = DataFrame(icnames[1] => icvals1, icnames[2] => icvals2, icnames[3] => icvals3, icnames[4] => icvals4,
                            :num_oscillatory_points => num_oscillatory_points_array, 
                            :average_period => average_periods, :maximum_period => maximum_periods, :minimum_period => minimum_periods,
                            :average_amplitude => average_amplitudes, maximum_amplitude => maximum_amplitudes, :minimum_amplitude => minimum_amplitudes)
    return results_df                
end




function loop_4fixedICs_thru_DFvals(paramconstraints::ParameterConstraints, icconstraints::InitialConditionConstraints, prob::ODEProblem; rangelength::Int = 4, DFrange = [100.,1000.,10000.])
    for DF in DFrange
        results_df = fixed_quadruplet_ic_searcher(paramconstraints, icconstraints, prob; rangelength=rangelength, fixedDF=DF)
        CSV.write("./ROCKFISH_DATA/4Fixed/SummaryResults/Summary_DF=$(round(fixedDF)).csv", results_df)
    end
end



function run_4fixedIC()
    tspan = (0., 2000.0)
    fullrn = make_fullrn()
    ogprob = ODEProblem(fullrn, [], tspan, [])

    de = modelingtoolkitize(ogprob)

    ogprobjac = ODEProblem(de, [], tspan, jac=true)


    #* Optimization of parameters to produce data for CSV
    param_constraints = define_parameter_constraints(;karange = (1e-3, 1e2), kbrange = (1e-3, 1e3), kcatrange = (1e-3, 1e3), dfrange = (1e2, 2e4))
    ic_constraints = define_initialcondition_constraints(; Lrange = (1e-1, 1e2), Krange = (1e-2, 1e2), Prange = (1e-2, 1e2), Arange = (1e-1, 1e2))

    # fixed_quadruplet_ic_searcher(param_constraints, ic_constraints, ogprobjac; rangelength=4, fixedDF=fixedDF)
    loop_4fixedICs_thru_DFvals(param_constraints, ic_constraints, ogprobjac; rangelength=10, DFrange = [100.,1000.,10000.])
end

run_4fixedIC()