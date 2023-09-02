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
    # using OrderedCollections
    # using Combinatorics
    # using LazySets, Polyhedra
    
    using ColorSchemes, Plots.PlotMeasures
    default(lw = 2, size = (1000, 600), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)
    # plotlyjs()
    # import CairoMakie as cm 
    # gr()
    # push!(LOAD_PATH, "../../UTILITIES")

    #* import the overloads for Evolutionary.jl
    include("../../UTILITIES/EvolutionaryOverloads.jl")

    #* import the Catalyst model "fullrn"
    include("../../UTILITIES/ReactionNetwork.jl")

    #* import the cost function and other evaluation functions
    include("../../UTILITIES/EvaluationFunctions.jl")
    # using .EvaluationFunctions

    #* import the genetic algorithm and associated functions
    include("../../UTILITIES/GA_functions.jl")

    #* import the plotting functions
    include("../../UTILITIES/TestBenchPlotUtils.jl")

    # include("../../UTILITIES/UnitTools.jl")

    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))
end


# debuglogger = ConsoleLogger(stderr, Logging.Debug)
# infologger = ConsoleLogger(stderr, Logging.Info)
# global_logger(infologger)

function make_ODE_problem()
    tspan = (0., 2000.0)

    fullrn = make_fullrn()

    ogprob = ODEProblem(fullrn, [], tspan, [])

    de = modelingtoolkitize(ogprob)

    ogprobjac::ODEProblem = ODEProblem(de, [], tspan, jac=true)
    return ogprobjac, ogprob
end




# ogjacsol = solve(ogprobjac, saveat=0.1, save_idxs= [6, 9, 10, 11, 12, 15, 16])

# tstart = cld(length(ogjacsol[1,:]), 1000)

# std(ogjacsol[1,end-tstart:end])
# # plot(ogjacsol)

# maxpeaks = findextrema(ogjacsol[1,:]; height=1e-2, distance=2)
# # minpeaks = findextrema(ogjacsol[1,:]; height=-1e-2, distance=2, find_maxima=false)
# # minpeaks = findextrema(flip_about_mean(ogjacsol[1,:]); height=-1e-2, distance=2, find_maxima=false)


# @btime CostFunction($ogjacsol)
# @code_warntype CostFunction(ogjacsol)


ogprobjac, ogprob = make_ODE_problem();



param_constraints = define_parameter_constraints(; karange = (1e-3, 1e2), kbrange = (1e-3, 1e3), kcatrange = (1e-3, 1e3), dfrange = (1e2, 2e4))
ic_constraints = define_initialcondition_constraints(; Lrange = (1e-1, 1e2), Krange = (1e-2, 1e2), Prange = (1e-2, 1e2), Arange = (1e-1, 1e2))

allconstraints = AllConstraints(param_constraints, ic_constraints)

fix_value!(allconstraints.DF, 1000.0)

gaproblem = GAProblem(allconstraints, ogprobjac)




fitfunc = make_fitness_function(gaproblem)

@code_warntype make_fitness_function(gaproblem)


@btime generate_population(allconstraints, 5000)

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


# Modification to make_fitness_function_with_fixed_inputs function
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


        insert!(new_input, 13, fixedDF)

        return evalfunc(new_input, prob)
    end
    return fitness_function
end



function test_fixedparam(param1::String, param2::String, param3::String, constraints::ConstraintType, prob::ODEProblem; fixed_values = [0.001, 0.01, 0.001], fixedDF=1000.)
    variable_constraints = deepcopy(constraints)
    # fixedtrip = [x for x in variable_constraints.ranges if x.name == param1 || x.name == param2 || x.name == param3]
    filter!(x -> x.name != param1 && x.name != param2 && x.name != param3, variable_constraints.ranges)
    # filter!(x -> x.name != param1 && x.name != param2 && x.name != param3, variable_constraints.icranges)


    #* remove DF from range because it'll be fixed
    filter!(x -> x.name != "DF", variable_constraints.ranges)

    fixed_ga_problem = GAProblem(variable_constraints, prob)
    fixedtrip_idxs = find_indices(param1, param2, param3, constraints.ranges) 

    # make_fitness_function_closure(evalfunc,prob; fitidx) = make_fitness_function_with_fixed_inputs(evalfunc, prob, fixed_values, fixedtrip_idxs; fitidx)
    make_fitness_function_closure(evalfunc,prob) = make_fitness_function_with_fixed_inputs_bothparamsIC(evalfunc, prob, fixed_values, fixedtrip_idxs, fixedDF)


    Random.seed!(1234)

    oscillatory_points_df = make_ga_dataframe(run_GA(fixed_ga_problem, make_fitness_function_closure; population_size = 10000, iterations = 5), prob, fixedDF) 
    num_oscillatory_points = length(oscillatory_points_df.fit)
    @info num_oscillatory_points


    return oscillatory_points_df
end

param_triplet = ["L", "K", "P"]
testfixed_df = test_fixedparam(param_triplet..., allconstraints, ogprob; fixed_values = [100.0,1.,1.])

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

#!NOTES 
#* Need to play around with PM mutation scheme. Not working as well as BGA, but has benefits
#* Adapt some sort of PM scheme for the initial population generation. Need to make sure it's diverse enough
#* Need to make this GA more rigorous as far as what it tells about a regime. Having more oscillatory points but where each point isn't that different from the next isn't a good comparitive metric.
#* DE instead of GA? Want global optimization and search, not trapped in local minima
function testbench(constraints::ConstraintType, prob::ODEProblem)
    test_gaproblem = GAProblem(constraints, prob)
    Random.seed!(1234)
    test_results = run_GA(test_gaproblem; population_size = 5000, iterations = 5, show_trace=true)

    # avg_fitness = mean(test_results.fitvals)
    # @info "Average fitness: $avg_fitness"
    # avg_period = mean(test_results.periods)
    # @info "Average period: $avg_period"
    # avg_amplitude = mean(test_results.amplitudes)
    # @info "Average amplitude: $avg_amplitude"
    if length(test_results.fitvals) == 0
        @info "No oscillatory points found."
    else 
        test_results_df = make_ga_dataframe(test_results, prob)
        return test_results_df#, avg_fitness, avg_period, avg_amplitude
    end
end

#* now testing Rosenbrock23 and new peak finder in getPerAmp for the ringing solutions
@code_warntype testbench(param_constraints, ogprobjac)

ogprobjac = remake(ogprobjac, u0 = [[100., 0.2, 0.2, 4.64]; zeros(12)])

test_results_df = testbench(param_constraints, ogprobjac)

test_results_df = testbench(allconstraints, ogprobjac)


function testBGA(valrange::Vector, m::Int = 2)
    prob = 1.0 / m
    function mutation(recombinant::T;
                      rng::AbstractRNG=Random.default_rng()
                     ) where {T <: AbstractVector}
        d = length(recombinant)
        @assert length(valrange) == d "Range matrix must have $(d) columns"
        δ = zeros(m)
        for i in 1:length(recombinant)
            for j in 1:m
                δ[j] = (rand(rng) < prob) ? δ[j] = 2.0^(-j) : 0.0
            end
            if rand(rng, Bool)
                recombinant[i] += sum(δ)*valrange[i]
            else
                recombinant[i] -= sum(δ)*valrange[i]
            end
        end
        return recombinant
    end
    return mutation
end

valrange = fill(2.0, 13)

mutationfunc = testBGA(valrange, 1)

params = copy(ogprob.p)

mutationfunc(params)


plotboth(test_results_df[9,:], ogprob)

@btime testbench($param_constraints, $ogprobjac)

plot_everything(test_results_df, ogprob; setnum=16, label="Testing", jump = 10)


stats_df = describe(test_results_df)
show(stats_df, allrows=true)

using StatsPlots

@df stats_df plot(:min, :max)

test_results_df.amp_percentage = test_results_df.amp./test_results_df.A


# split_dataframe!(test_results_df, ogprobjac)

CSV.write("OscillatorPaper/FigureGenerationScripts/high_amp_Amem.csv", test_results_df)

@btime testbench($param_constraints, $ogprob)
@btime testbench($param_constraints, $ogprobjac)

plot_everything(test_results_df, ogprob; setnum=15, label="TestingSTDWindow", jump = 10)


plotboth(test_results_df[1,:], ogprob)



#* measure A in solution vs A membrane 
#* quadruplet fixed search initial conditions 


# get data 
testdf = CSV.read("OscillatorPaper/FigureGenerationScripts/test.csv", DataFrame)

#combine all parameter columns into one column of vectors



p = [param for param in test_results_df[1763, Between(:ka1, :DF)]]
u0 = [ic for ic in test_results_df[1763, Between(:L, :A)]]

reprob = remake(ogprob, p = p, u0 = [u0; zeros(length(ogprob.u0) - length(u0))])

sol = solve(reprob, Rosenbrock23(), saveat=0.1, save_idxs = [6, 9, 10, 11, 12, 15, 16])

Amem = map(sum, sol.u)

findextrema(Amem; height=1e-2, distance=2)
findextrema(Amem; height=0.0, distance=2, find_maxima=false)

CostFunction(sol)
plot(sol)




testdf = CSV.read("/Users/jonathanfischer/Desktop/PhD_ThesisWork/Julia/OscillatorPaper/FigureGenerationScripts/DF=100.0.csv", DataFrame)

plot_everything(testdf, ogprob; setnum=15, label="TestingSTDWindow", jump = 10)