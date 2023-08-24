begin 
    using Plots; #theme(:juno)
    # using Compose
    using Catalyst
    using DifferentialEquations, ModelingToolkit
    using Statistics
    # using Peaks
    using FindPeaks1D
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

    include("../../UTILITIES/TestBenchPlotUtils.jl")

    # include("../../UTILITIES/UnitTools.jl")


    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))
end


# debuglogger = ConsoleLogger(stderr, Logging.Debug)
# infologger = ConsoleLogger(stderr, Logging.Info)
# global_logger(infologger)

tspan = (0., 2000.0)
fullrn = make_fullrn()
# ogsys = convert(ODESystem, fullrn)
# @unpack L, K, P, A = ogsys
ogprob = ODEProblem(fullrn, [:L => 10.0], tspan, [])
de = modelingtoolkitize(ogprob)

ogprobjac = ODEProblem(de, [], tspan, jac=true)




ogjacsol = solve(ogprobjac, saveat=0.1)
plot(ogjacsol)

maxpeaks = findextrema(ogjacsol[1,:]; height=1e-2, distance=2)
minpeaks = findextrema(ogjacsol[1,:]; height=-1e-2, distance=2, find_maxima=false)
minpeaks = findextrema(flip_about_mean(ogjacsol[1,:]); height=-1e-2, distance=2, find_maxima=false)





@btime CostFunction($ogjacsol)
@code_warntype CostFunction(ogjacsol)


@btime ogjacsol[1:100]
plot(ogjacsol)

flipped = flip_about_mean(ogjacsol[1,:])
plot(flipped)
flip_about_mean!(ogjacsol[1,:])

ogjacsol[1:end] .= 0.1

@btime solve(ogprob, saveat=0.1, save_idxs=4)
@btime solve(ogprobjac, saveat=0.1, save_idxs=4)



# ogsol = solve(ogprob, saveat=0.1, save_idxs=4)

# fftdata = getFrequencies(ogsol[1,:])
# frequencies_per_minute!(ogsol.t, fftdata)
# fftdata

# cld(length(fftdata), 1000)
# fft_peakindexes, fft_peakvals = findmaxima(fftdata,1) #* get the indexes of the peaks in the fft
# fft_peakindexes, peakprops = findpeaks1d(fftdata; height = 0.0, distance = 1) #* get the indexes of the peaks in the fft

param_constraints = define_parameter_constraints(ogprobjac; karange = (1e-3, 1e2), kbrange = (1e-2, 1e3), kcatrange = (1e-2, 1e3), dfrange = (1e2, 1e5))
ic_constraints = define_initialcondition_constraints(ogprobjac; Lrange = (1e-2, 1e2), Krange = (1e-2, 1e2), Prange = (1e-2, 1e2), Arange = (1e-2, 1e2))

allconstraints = AllConstraints(param_constraints, ic_constraints)

testgaproblem = GAProblem(allconstraints, ogprobjac)

generate_population(testgaproblem.constraints, 1000)

# Modification to make_fitness_function_with_fixed_inputs function
function make_fitness_function_with_fixed_inputs(evalfunc::Function, prob::ODEProblem, fixed_input_triplet::Vector{Float64}, triplet_idxs::Tuple{Int, Int, Int}; fitidx = 4)
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

        return evalfunc(new_input, prob; idx = fitidx)
    end
    return fitness_function
end


# Modification to make_fitness_function_with_fixed_inputs function
function make_fitness_function_with_fixed_inputs_bothparamsIC(evalfunc::Function, prob::ODEProblem, fixed_input_triplet::Vector{Float64}, triplet_idxs::Tuple{Int, Int, Int}, fixedDF=1000.; fitidx = 4)
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

        return evalfunc(new_param_input, new_ic_input, prob; idx = fitidx)
    end
    return fitness_function
end



function test_fixedparam(param1::String, param2::String, param3::String, constraints::ConstraintType, prob::ODEProblem; fixed_values = [0.001, 0.01, 0.001], fixedDF=1000.)
    Random.seed!(1234)
    variable_constraints = deepcopy(constraints)
    # fixedtrip = [x for x in variable_constraints.ranges if x.name == param1 || x.name == param2 || x.name == param3]
    filter!(x -> x.name != param1 && x.name != param2 && x.name != param3, variable_constraints.ranges)
    # filter!(x -> x.name != param1 && x.name != param2 && x.name != param3, variable_constraints.icranges)


    #* remove DF from range because it'll be fixed
    filter!(x -> x.name != "DF", variable_constraints.ranges)

    fixed_ga_problem = GAProblem(variable_constraints, prob)
    fixedtrip_idxs = find_indices(param1, param2, param3, constraints.ranges) 

    # make_fitness_function_closure(evalfunc,prob; fitidx) = make_fitness_function_with_fixed_inputs(evalfunc, prob, fixed_values, fixedtrip_idxs; fitidx)
    make_fitness_function_closure(evalfunc,prob; fitidx) = make_fitness_function_with_fixed_inputs_bothparamsIC(evalfunc, prob, fixed_values, fixedtrip_idxs, fixedDF; fitidx)



    oscillatory_points_df = make_df(run_GA(fixed_ga_problem, make_fitness_function_closure; population_size = 10000, iterations = 5)) 
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
        insert!(ind, 13, fixedDF)
    end
    # return oscillatory_points_df
    #* split parameter values into separate columns and add initial conditions
    # split_dataframe!(oscillatory_points_df, prob)
    return oscillatory_points_df
end

param_triplet = ["ka1", "kcat1", "kcat7"]
testfixed_df = test_fixedparam(param_triplet..., allconstraints, ogprob; fixed_values = [0.1,1000.,1000.])
split_dataframe!(testfixed_df, ogprob)
CSV.write("OscillatorPaper/FigureGenerationScripts/test.csv", testfixed_df)


plot_everything(testfixed_df, ogprob; setnum=13, label="Params&Inits", jump=50)

testfixed_df = CSV.read("OscillatorPaper/FigureGenerationScripts/test.csv", DataFrame)
testsol = solve(remake(ogprob, p = testfixed_df.ind[401]), saveat=0.1)

newprob = remake(ogprob, u0 = [[97.23405302393168, 76.78427415397275, 8.928095154587805, 22.611353860124073]; zeros(12)],p = [0.1, 249.16082412239868, 1000.0, 4.969336136315573, 115.12839147825775, 3.399727998677895, 50.40303569894044, 0.001, 0.01, 34.13205544050995, 903.6369655144391, 1000.0, 1000.0])
testsol = solve(newprob, saveat=0.1)
plotsol(testsol)
plotboth(401, testfixed_df, ogprob)


plotboth(row) = plotboth(row, testfixed_df, ogprob)
plotboth(1)

test_fitness(row) = eval_param_fitness(testfixed_df.ind[row], ogprob)
test_fitness(2)
#########################

#* OBSERVATIONS FROM THE TESTBENCH
#* No false positives or negatives currently with no fixed param GA optimization
#* Bias towards high frequency I think
#* Last half STD check seems necessary, not sure about trim
#* Not much difference between Rosenbrock23 and Tsit5, or dt = 0.1 and dt = 0.01
#* Peak height threshold seems to be important


#< Regular GA testbench
# Random.seed!(1234)
test_gaproblem = GAProblem(param_constraints, ogprob)
test_results = run_GA(test_gaproblem; population_size = 5000, iterations = 5, fitidx = 4)
@code_warntype run_GA(test_gaproblem; population_size = 5000, iterations = 5, fitidx = 4)
typeof(test_results)
test_results_df = trace_to_df(test_results)


function foo(gaprob::GAProblem)
    test_results = run_GA(gaprob; population_size = 5000, iterations = 5, fitidx = 4)
    results_df = trace_to_df(test_results)
    mean_fitness = mean(results_df.fit)
    mean_period = mean(results_df.per)
    mean_amplitude = mean(results_df.amp)
    return mean_fitness, mean_period, mean_amplitude
end

@code_warntype foo(test_gaproblem)

# avg_fitness = mean(test_results.fit)
# avg_period = mean(test_results.per)
# avg_amplitude = mean(test_results.amp)

# plotboth(row) = plotboth(row, test_results, ogprob)
# plotboth(1338)

# testprob = remake(ogprob, p = test_results.ind[1338])
# testsol = solve(testprob, saveat=0.1)
# fftdata = getFrequencies(testsol[4,:])
# normalize_time_series!(fftdata)
# plot(fftdata, xlims=(0,100))
# plotfft(testsol; fitidx=4)
# fft_peakindexes, peakprops = findpeaks1d(fftdata; height = 1e-3, distance = 2) #* get the indexes of the peaks in the fft
# fft_peakvals = peakprops["peak_heights"]


function testbench(param_constraints::ConstraintType, prob::ODEProblem)
    test_gaproblem = GAProblem(param_constraints, prob)
    Random.seed!(1234)
    test_results = run_GA(test_gaproblem; population_size = 5000, iterations = 5, fitidx = 4, show_trace=true)

    avg_fitness = mean(test_results.fitvals)
    # @info "Average fitness: $avg_fitness"
    avg_period = mean(test_results.periods)
    # @info "Average period: $avg_period"
    avg_amplitude = mean(test_results.amplitudes)
    # @info "Average amplitude: $avg_amplitude"
    test_results_df = make_df(test_results)
    return test_results_df#, avg_fitness, avg_period, avg_amplitude
end

#* now testing Rosenbrock23 and new peak finder in getPerAmp for the ringing solutions
@code_warntype testbench(param_constraints, ogprobjac)
test_results_df = testbench(param_constraints, ogprobjac)

@btime testbench($param_constraints, $ogprob)
@btime testbench($param_constraints, $ogprobjac)

plot_everything(test_results_df, ogprob; setnum=13, label="Params&Inits", jump = 10)



#* measure A in solution vs A membrane 
#* quadruplet fixed search initial conditions 


# get data 
testdf = CSV.read("OscillatorPaper/FigureGenerationScripts/test.csv", DataFrame)

#combine all parameter columns into one column of vectors






