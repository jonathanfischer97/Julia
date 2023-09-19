begin 
    using Plots; #theme(:juno)
    using Catalyst
    using OrdinaryDiffEq, ModelingToolkit
    using Statistics
    using Evolutionary, FFTW
    using Random
    using Distributions
    using DataFrames#, DataFrameMacros
    using CSV
    using StaticArrays
    using BenchmarkTools, ProgressMeter

    # using JET

    using LinearAlgebra

    # using Setfield
    
    using ColorSchemes, Plots.PlotMeasures
    default(lw = 2, size = (1000, 600), dpi = 200, bottom_margin = 12px, left_margin = 16px, top_margin = 10px, right_margin = 8px)

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


    const SHOW_PROGRESS_BARS = parse(Bool, get(ENV, "PROGRESS_BARS", "true"))

    FFTW.set_num_threads(18)
end

AutoTsit5(Rodas5P())

# debuglogger = ConsoleLogger(stderr, Logging.Debug)
# infologger = ConsoleLogger(stderr, Logging.Info)
# global_logger(infologger)
@unpack ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, DF = osys
@unpack L, K, P, A, Lp, LpA, LK, LpP, LpAK, LpAP, LpAKL, LpAPLp, AK, AP, AKL, APLp = osys

fullrn = make_fullrn()

# oprob = ODEProblem(fullrn, [], (0., 1000.), ())

osys = convert(ODESystem, fullrn)

oprob = ODEProblem(osys, [], (0., 2000.), ())

osol = solve(oprob, AutoTsit5(Rodas5P()), saveat=0.1)

@btime solve_for_fitness_peramp(oprob, [6, 9, 10, 11, 12, 15, 16])
plotboth(osol)


tspans = logrange(10, 100000, 20)
costvals = []

for tend in tspans
    reprob = remake(oprob, tspan = (0., round(tend)))
    sol = solve(reprob, AutoTsit5(Rodas5P()), saveat=0.1, save_idxs=[6, 9, 10, 11, 12, 15, 16])
    cost = CostFunction(sol)
    push!(costvals, cost[1])
end

plot(tspans, costvals, xaxis=:log10, xlabel="tspan", ylabel="Cost", label="Cost vs tspan")





@btime solve($oprob, AutoTsit5(Rodas5P()))
@btime solve($oprob, Rosenbrock23())
@btime solve($oprob, Rodas5P())

osol[L]

ogprobjac = make_ODE_problem();

sol = solve(ogprobjac, Rodas5())
sol[L]


paramconstraints = ParameterConstraints(; karange = (1e-3, 1e2), kbrange = (1e-3, 1e3), kcatrange = (1e-3, 1e3), dfrange = (1e2, 2e4))
ic_constraints = InitialConditionConstraints(; Lrange = (1e-1, 1e2), Krange = (1e-2, 1e2), Prange = (1e-2, 1e2), Arange = (1e-1, 1e2))

#* calculate the minimum tspan in order to capture largest period oscillations
minimum_tspan = 1/(0.1 * 0.1) #multiply slowest forward rate by lowest concentration and take reciprocal to get seconds 

allconstraints = AllConstraints()
gaproblem = GAProblem(allconstraints, ogprobjac)
# @report_opt GAProblem(allconstraints, ogprobjac)

# @code_warntype GAProblem(allconstraints, ogprobjac)


initial_population = generate_population(allconstraints, 50000)
ga_results = run_GA(gaproblem, initial_population; iterations = 5)
ga_df = make_ga_dataframe(ga_results, allconstraints)
desc_df = describe(ga_df)
show(desc_df, allrows=true)
save_to_csv(ga_results, allconstraints, "adaptive_tspan_test.csv")

function test(gaproblem)
    Random.seed!(1234)
    initial_population = generate_empty_population(gaproblem.constraints, 3000)
    population = generate_population_stratified!(initial_population, gaproblem.constraints,10)
    run_GA(gaproblem, population; iterations = 5)
end

ga_results = test(gaproblem)

@btime test($gaproblem)

@bprofile test(gaproblem)

#* testing GA output for when fixed K=0.1, P=0.1, A=0.1, DF=100
ogprobjac = make_ODE_problem(100.);

gaprob = GAProblem(ode_problem = ogprobjac)
gaconstraints = gaprob.constraints

gaconstraints.DF.isfixed = true
gaconstraints.DF.fixed_value = 1000.0

gaconstraints.K.isfixed = true
gaconstraints.K.fixed_value = 0.01

gaconstraints.P.isfixed = true
gaconstraints.P.fixed_value = 0.1

gaconstraints.A.isfixed = true
gaconstraints.A.fixed_value = 0.1

gaprob
Random.seed!(1234)

initial_population = generate_population(gaconstraints, 10000)

ga_results = run_GA(gaprob, initial_population; iterations = 5)

#! population generation testing ##########
Random.seed!(1234)
loguniform_pop = generate_population(gaproblem.constraints, 5000)

initial_population = generate_empty_population(gaproblem.constraints, 5000)
stratified_pop = generate_population_stratified!(initial_population, gaproblem.constraints,3)

loguniform_df = make_pop_dataframe(loguniform_pop, gaproblem.constraints)
stratified_df = make_pop_dataframe(stratified_pop, gaproblem.constraints)


loguniform_stats = describe(loguniform_df)
show(loguniform_stats, allrows=true)
filter!(row -> row[1] != :DF, loguniform_stats)
CSV.write("loguniform_stats.csv", loguniform_df)

stratified_stats = describe(stratified_df)
show(stratified_stats, allrows=true)
filter!(row -> row[1] != :DF, stratified_stats)
CSV.write("stratified_stats.csv", stratified_df)

using StatsPlots

function plot_boxplots(df1::DataFrame, df2::DataFrame, title1::String, title2::String)
    p = plot(layout=(1,2), size=(1000, 500))
    @df df1 bar!(p[1], cols(16), title=title1, label = "P", yscale=:log10)
    @df df2 bar!(p[2], cols(16), title=title2, label = "P", yscale=:log10)
    display(p)
end


function plot_scatter(df1::DataFrame, df2::DataFrame)
    p = plot(size=(1000, 500), legend=:topright)
    for col in names(df1)
        if col != :variable
            scatter!(p, ["Min", "Max", "Median"], [df1[col, :min], df1[col, :max], df1[col, :median]], label="$col - Method 1")
            scatter!(p, ["Min", "Max", "Median"], [df2[col, :min], df2[col, :max], df2[col, :median]], label="$col - Method 2", shape=:diamond)
        end
    end
    xlabel!(p, "Statistic")
    ylabel!(p, "Value")
    title!(p, "Comparison of Min, Max, and Median Values")
    display(p)
end


plot_boxplots(loguniform_df, stratified_df, "Loguniform Population", "Stratified Population")



function shannon_diversity(population::Vector{Vector{Float64}}, nbins::Int)
    n_params = length(population[1])
    shannon_indices = zeros(n_params)
    
    for i in 1:n_params
        param_values = [indiv[i] for indiv in population]
        counts = fit(Histogram, param_values, nbins:nbins).weights
        p = counts ./ sum(counts)
        shannon_indices[i] = -sum(p .* log.(p))
    end
    
    return shannon_indices
end





using Profile
Profile.print(format=:tree, mincount=100)

save_to_csv(ga_results, allconstraints, "gentest.csv")



obj = Evolutionary.EvolutionaryObjective(threaded_fitness_function, test_input, zeros(3,1))

fv = zeros(3,1)
value(obj, fv, test_input)

@btime value($obj, $fv, $test_input)


#< START
fixed_inputs = (L = 100.0, K = 1.0, P = 1.0, A = 10.0)

set_fixed_constraints!(gaproblem.constraints, fixed_inputs...)


#* test objective function 
fitness_function = make_fitness_function(gaproblem.constraints, ogprobjac)

test_input = [ogprobjac.p; ogprobjac.u0[1:4]]

@btime fitness_function($test_input)

@code_warntype fitness_function(test_input)

fitness_function(test_input)


###############
threaded_fitness_function = make_fitness_function_threaded(gaproblem.constraints, ogprobjac)
@btime threaded_fitness_function($test_input)

inplace_fitness_function = make_fitness_function_inplace(gaproblem.constraints, ogprobjac)
@btime inplace_fitness_function($test_input)

testpop = ga_results.population

F = zeros(Float64, (3, length(testpop)))

function testfunc(F,pop)
    Threads.@threads for i in 1:length(pop)
        F[:,i] .= threaded_fitness_function(pop[i])
    end
    F
end

BLAS.set_num_threads(18)
FFTW.set_num_threads(18)

@btime testfunc($F, $testpop)

"""Non-threaded: 36.736 s (161490695 allocations: 16.52 GiB)"""





"""
    Takes a named tuple of fixed inputs and a GAProblem and sets the constraints of the GAProblem to the fixed inputs.
    Returns `DataFrame` of the optimization results.
"""
function test_fixedparam(gaprob::GAProblem; fixed_inputs)

    constraints = gaprob.constraints

    set_fixed_constraints!(constraints; fixed_inputs)


    Random.seed!(1234)

    initial_population = generate_population(constraints, 5000)

    ga_results = run_GA(gaprob, initial_population; iterations = 5)

    oscillatory_points_df = make_ga_dataframe(ga_results, constraints) 
    num_oscillatory_points = nrow(oscillatory_points_df)
    @info num_oscillatory_points

    return oscillatory_points_df
end


fixed_inputs = (L = 100.0, K = 1.0, P = 1.0, A = 10.0)

testfixed_df = test_fixedparam(gaproblem; fixed_inputs)

@btime test_fixedparam($gaproblem; fixed_inputs)
#*  112.566 s (1484360277 allocations: 196.73 GiB)

#* when eval is broadcasted
#* 118.019 s (1484282110 allocations: 196.73 GiB)





set_fixed_constraints!(allconstraints, fixed_inputs)

fitfunc = gaproblem.fitness_function

input = ogprob.p

fitfunc(input)
#########################

#* OBSERVATIONS FROM THE TESTBENCH
#* No false positives or negatives currently with no fixed param GA optimization
#* Bias towards high frequency I think
#* Last half STD check seems necessary, not sure about trim
#* Not much difference between Rosenbrock23 and Tsit5, or dt = 0.1 and dt = 0.01
#* Peak height threshold seems to be important





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

