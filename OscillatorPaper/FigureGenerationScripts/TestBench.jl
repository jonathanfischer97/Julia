begin 
    using Plots; #theme(:juno)
    # using Compose
    using Catalyst
    using DifferentialEquations
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


# debuglogger = ConsoleLogger(stderr, Logging.Debug)
# infologger = ConsoleLogger(stderr, Logging.Info)
# global_logger(infologger)

fullrn = make_fullrn()
ogprob = ODEProblem(fullrn, [], (0.,2000.0), [])
ogsol = solve(ogprob, saveat=0.1)

fftdata = getFrequencies(ogsol[1,:])
fft_peakindexes, fft_peakvals = findmaxima(fftdata,1) #* get the indexes of the peaks in the fft
fft_peakindexes, peakprops = findpeaks1d(fftdata; height = 0.0, distance = 1) #* get the indexes of the peaks in the fft

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



function test_fixedparam(param1::String, param2::String, param3::String, constraints::ConstraintType, prob::ODEProblem; fixed_values = [0.001, 0.01, 0.001])
    variable_constraints = deepcopy(constraints)
    # fixedtrip = [x for x in variable_constraints.ranges if x.name == param1 || x.name == param2 || x.name == param3]
    filter!(x -> x.name != param1 && x.name != param2 && x.name != param3, variable_constraints.ranges)

    fixed_ga_problem = GAProblem(variable_constraints, prob)
    fixedtrip_idxs = find_indices(param1, param2, param3, constraints.ranges) 

    make_fitness_function_closure(evalfunc,prob) = make_fitness_function_with_fixed_inputs(evalfunc, prob, fixed_values, fixedtrip_idxs)


    oscillatory_points_df = run_GA(fixed_ga_problem, make_fitness_function_closure; population_size = 5000, iterations = 5) 
    num_oscillatory_points = length(oscillatory_points_df.ind)
    @info num_oscillatory_points

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
    # split_dataframe!(oscillatory_points_df, prob)
    return oscillatory_points_df
end

param_triplet = ["ka2", "kb2", "ka4"]
testfixed_df = test_fixedparam(param_triplet..., param_constraints, ogprob)
CSV.write("OscillatorPaper/FigureGenerationScripts/testbench.csv", testfixed_df)


function plot_everything(df::DataFrame, prob::ODEProblem)
    for i in 1:5:nrow(df)
        p = plotboth(i, df, prob)
        savefig(p, "OscillatorPaper/FigureGenerationScripts/TestbenchPlots/plot_$(i).png")
    end
end

plot_everything(testfixed_df, ogprob)


plotboth(row) = plotboth(row, testfixed_df, ogprob)
plotboth(1)

for i in 1:50:nrow(testfixed_df)
    plotboth(i)
end

test_fitness(row) = eval_param_fitness(testfixed_df.ind[row], ogprob)
test_fitness(2)
#########################