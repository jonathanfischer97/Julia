begin 
    using Plots; #theme(:juno)
    using Catalyst
    using DifferentialEquations, ModelingToolkit
    using Statistics
    # using Peaks
    # using FindPeaks1D
    using Evolutionary, FFTW
    using Random
    using Distributions
    using DataFrames#, DataFrameMacros
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


# @code_warntype make_fitness_function(eval_param_fitness, ogprobjac; fitidx = 4)

# @btime solve($ogprob, saveat = 0.1, save_idxs = 1)
# @btime solve($ogprobjac, saveat = 0.1, save_idxs = 1)

# newprob = remake(ogprob, p = ogprob.p .* 1.0)
# @btime solve($newprob, Rosenbrock23(), saveat = 0.1, save_idxs = 1)

# newprobjac = remake(ogprobjac, p = ogprobjac.p .* 1.0)
# @btime solve($newprobjac, Rosenbrock23(), saveat = 0.1, save_idxs = 1)



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
param_constraints = define_parameter_constraints(ogprob; karange = (1e-3, 1e2), kbrange = (1e-2, 1e3), kcatrange = (1e-2, 1e3), dfrange = (1e2, 1e5))
ic_constraints = define_initialcondition_constraints(ogprob; Lrange = (1e-2, 1e2), Krange = (1e-2, 1e2), Prange = (1e-2, 1e2), Arange = (1e-2, 1e2))


#* Function loops through 4D grid of different initial conditions, letting all parameters be freely optimized, and saves the results to a csv file
function fixed_quadruplet_ic_searcher(paramconstraints::ParameterConstraints, icconstraints::InitialConditionConstraints, prob::ODEProblem; rangelength = 4)
    #* get the ranges of the initial conditions
    icranges = [logrange(constraints.min, constraints.max, rangelength) for constraints in icconstraints.ranges]

    icnames = [constraints.name for constraints in icconstraints.ranges]

    num_rows = rangelength^length(icnames)

    # results_df = DataFrame(icnames[1] => Vector{Float64}(undef, num_rows), icnames[2] => Vector{Float64}(undef, num_rows), icnames[3] => Vector{Float64}(undef, num_rows),icnames[4] => Vector{Float64}(undef, num_rows),
    #                         "num_oscillatory_points" => Vector{Int}(undef, num_rows), 
    #                         "average_period" => Vector{Float64}(undef, num_rows), "maximum_period"=>Vector{Float64}(undef, num_rows), "minimum_period"=>Vector{Float64}(undef, num_rows),
    #                         "average_amplitude" => Vector{Float64}(undef, num_rows), "maximum_amplitude"=>Vector{Float64}(undef, num_rows), "minimum_amplitude"=>Vector{Float64}(undef, num_rows))

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

    path = mkpath("./OscillatorPaper/FigureGenerationScripts/4FixedICRawSets/")


    #* loop through each ic range and run the GA on each set of initial conditions after remaking the problem with them
    for icval1 in icranges[1]
        for icval2 in icranges[2]
            for icval3 in icranges[3]
                for icval4 in icranges[4]
                    icvals = [icval1, icval2, icval3, icval4]
                    @info icvals

                    #* remake the problem with the new initial conditions
                    newprob = remake(prob, u0 = [icvals; prob.u0[5:end]])
                    
                    #* make new GA problem with new initial conditions
                    ga_problem = GAProblem(paramconstraints, newprob)

                    #* set seed for reproducibility
                    Random.seed!(1234)

                    #* run the GA on the new problem
                    oscillatory_points_results = run_GA(ga_problem; population_size = 5000, iterations = 5)

                    #* get the number of oscillatory points
                    num_oscillatory_points = length(oscillatory_points_results.population)

                    #* if there are no oscillatory points, save the results to the results_df and continue
                    if iszero(num_oscillatory_points)
                        # results_df[i, :] = (icval1, icval2, icval3, icval4, 0, NaN, NaN, NaN, NaN, NaN, NaN)
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
                        # results_df[i, :] = (icval1, icval2, icval3, icval4, num_oscillatory_points, average_period, maximum_period, minimum_period, average_amplitude, maximum_amplitude, minimum_amplitude)
                        icvals1[i] = icval1
                        icvals2[i] = icval2
                        icvals3[i] = icval3
                        icvals4[i] = icval4
                        num_oscillatory_points_array[i] = num_oscillatory_points
                        

                    
                        #* make dataframe from oscillatory_points_results
                        oscillatory_points_df = make_df(oscillatory_points_results)
                        
                        #* split parameter values into separate columns and add initial conditions
                        split_dataframe!(oscillatory_points_df, newprob)

                        #* rewrite the L, K, P, A columns with the initial conditions
                        oscillatory_points_df.L .= icval1
                        oscillatory_points_df.K .= icval2
                        oscillatory_points_df.P .= icval3
                        oscillatory_points_df.A .= icval4

                        CSV.write(path*"$(round(icval1; digits = 2))_$(round(icval2;digits = 2))_$(round(icval3; digits=2))_$(round(icval4; digits=2)).csv", oscillatory_points_df)
                    end
                    next!(loopprogress)
                    i += 1
                end
            end
        end
    end
    results_df = DataFrame(icnames[1] => icvals1, icnames[2] => icvals2, icnames[3] => icvals3, icnames[4] => icvals4,
                            "num_oscillatory_points" => num_oscillatory_points_array, 
                            "average_period" => average_periods, "maximum_period"=>maximum_periods, "minimum_period"=>minimum_periods,
                            "average_amplitude" => average_amplitudes, "maximum_amplitude"=>maximum_amplitudes, "minimum_amplitude"=>minimum_amplitudes)
    CSV.write("./OscillatorPaper/FigureGenerationScripts/4FixedICs.csv", results_df)
    return results_df                
end

df = fixed_quadruplet_ic_searcher(param_constraints, ic_constraints, ogprobjac; rangelength=3)





scatter3d(df.Kinase, df.Phosphatase, df.AP2, color = df.num_oscillatory_points, zlims = (0, 100), colorbar_title = "Number of Oscillatory Points", title = "Number of Oscillatory Points vs. Initial Conditions", xlabel = "L", ylabel = "K", zlabel = "P", markersize = df.maximum_period, markerstrokewidth = 0, colormap = :viridis, camera = (30, 30),
            xscale = :log10, yscale = :log10)  