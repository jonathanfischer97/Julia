begin 
    using Plots 
    using Catalyst
    using DifferentialEquations
    using Statistics
    using Peaks
    using Evolutionary, FFTW
    using Random
    using Distributions
    # using DataFrames
    # using Unitful
    # using Unitful: µM, nm, s
    using StaticArrays
    using BenchmarkTools, Profile, ProgressMeter
    # using MultivariateStats, UMAP, TSne, StatsPlots
    # using GlobalSensitivity, QuasiMonteCarlo
    using LinearAlgebra
    # using BifurcationKit, Setfield, ForwardDiff, Parameters; const BK = BifurcationKit
    using ColorSchemes, Plots.PlotMeasures
    using OrderedCollections
    using Combinatorics
    # using LazySets, Polyhedra
    default(lw = 2, size = (1000, 600))
    # plotlyjs()
    # gr()
    # using Base.Threads
end

# import the Catalyst model "fullrn"
include("../../UTILITIES/ReactionNetwork.jl")

# import the cost function and other evaluation functions
include("../../UTILITIES/EvaluationFunctions.jl")

# import the genetic algorithm and associated functions
include("../../UTILITIES/GA_functions.jl")

#! Solve model for arbitrary oscillatory parameters and initial conditions
begin
    #? Parameter list
    psym = [:ka1 => 0.009433439939827041, :kb1 => 2.3550169939427845, :kcat1 => 832.7213093872278, :ka2 => 12.993995997539924, :kb2 => 6.150972501791291,
            :ka3 => 1.3481451097940793, :kb3 => 0.006201726090609513, :ka4 => 0.006277294665474662, :kb4 => 0.9250191811994848, :ka7 => 57.36471615394549, 
            :kb7 => 0.04411989797898752, :kcat7 => 42.288085868394326, :DF => 3631.050539219606]
    p = [x[2] for x in psym]
        
    #? Initial condition list
    usym = [:L => 3.0,  :K => 0.5, :P => 0.3, :A => 2.0,:Lp => 0.0, :LpA => 0.0, :LK => 0.0, 
            :LpP => 0.0, :LpAK => 0.0, :LpAP => 0.0, :LpAKL => 0.0, :LpAPLp => 0.0, :AK => 0.0, :AP => 0.0, 
            :AKL => 0.0, :APLp => 0.0]
    u0 = [x[2] for x in usym]

    #? Timespan for integration
    tspan = (0., 100.)

    #? Create ODE problem and solve
    fullprob = ODEProblem(fullrn, u0, tspan, p)
    # sol = solve(fullprob, saveat=0.1, save_idxs=1)

    #? Plot the results
    # plot(sol)
end



begin
    ## parameter constraint ranges ##
    ka_min, ka_max = 0.001, 10. #uM^-1s^-1
    kb_min, kb_max = 0.001, 1000.0 #s^-1
    kcat_min, kcat_max = 0.001, 1000.0 #s^-1

    param_values = OrderedDict(
        "ka1" => Dict("min" => ka_min, "max" => ka_max),
        "kb1" => Dict("min" => kb_min, "max" => kb_max),
        "kcat1" => Dict("min" => kcat_min, "max" => kcat_max),
        "ka2" => Dict("min" => ka_min, "max" => ka_max),
        "kb2" => Dict("min" => kb_min, "max" => kb_max),
        "ka3" => Dict("min" => ka_min, "max" => ka_max),
        "kb3" => Dict("min" => kb_min, "max" => kb_max),
        "ka4" => Dict("min" => ka_min, "max" => ka_max),
        "kb4" => Dict("min" => kb_min, "max" => kb_max),
        "ka7" => Dict("min" => ka_min, "max" => ka_max),
        "kb7" => Dict("min" => kb_min, "max" => kb_max),
        "kcat7" => Dict("min" => kcat_min, "max" => kcat_max),
        "y" => Dict("min" => 10., "max" => 20000.)
        # "L" => Dict("min" => 0.1, "max" => 3.0),
        # "Lp" => Dict("min" => 0.1, "max" => 3.0),
        # "K" => Dict("min" => 0.1, "max" => 1.0),
        # "P" => Dict("min" => 0.1, "max" => 1.0),
        # "A" => Dict("min" => 0.1, "max" => 3.0),
    );
end

begin
    ic_values = OrderedDict(
        "L" => Dict("min" => 0.1, "max" => 3.0),
        "Lp" => Dict("min" => 0.1, "max" => 3.0),
        "K" => Dict("min" => 0.1, "max" => 1.0),
        "P" => Dict("min" => 0.1, "max" => 1.0),
        "A" => Dict("min" => 0.1, "max" => 3.0),
    );
end

"""Function that generates population of log-uniform sampled random parameter values"""
function generate_population(param_values::OrderedDict, n::Int)
    params = keys(param_values)
    n_params = length(params)
    population = Matrix{Float64}(undef, n_params, n)
    for i in 1:n
        for (ind, param) in enumerate(params)
            min_val = log(param_values[param]["min"])
            max_val = log(param_values[param]["max"])
            min_val < max_val ? population[ind, i] = exp(rand(Uniform(min_val, max_val))) : population[ind, i] = exp(min_val)
        end
    end
    return [collect(population[:, i]) for i in 1:n]
end

# using Evolutionary: value, value!, minimizer, initial_state, update_state!, trace!, evaluate!

#! OVERRIDES FOR Evolutionary.jl ##
"""Trace override function"""
function Evolutionary.trace!(record::Dict{String,Any}, objfun, state, population, method::GA, options)
    record["staterecord"] = [(ind=population[i], fit=state.fitpop[i]) for i in eachindex(population)]
end

"""Show override function to prevent printing large arrays"""
function Evolutionary.show(io::IO, t::Evolutionary.OptimizationTraceRecord)
    print(io, lpad("$(t.iteration)",6))
    print(io, "   ")
    print(io, lpad("$(t.value)",14))
    for (key, value) in t.metadata
        if !isa(value, AbstractArray)
            print(io, "\n * $key: $value")
        end
    end
    return
end







#! SINGLE RUN GENETIC ALGORITHM ##
"""Wrapper function to create a fitness function that includes your ODE problem as a constant"""
function make_fitness_function(func::Function, prob::ODEProblem)
    function fitness_function(p::Vector{Float64})
        return func(p, prob)
    end
    return fitness_function
end



fitness_function = make_fitness_function(eval_ic_fitness, fullprob) # Create a fitness function that includes your ODE problem as a constant

# using Debugger
#! Optimization block
function testGAfunc(evalfunc, fitness_function_factory, prob)
    population_size = 1000
    pop = generate_population(ic_values, population_size)

    myconstraints = BoxConstraints([ic_values[p]["min"] for p in keys(ic_values)], [ic_values[p]["max"] for p in keys(ic_values)])
    opts = Evolutionary.Options(abstol=1e-2, reltol=1.00, successive_f_tol = 5, iterations=5, store_trace = true, 
            show_trace=true, show_every=1, parallelization=:thread)
    common_range = 0.5; valrange = fill(common_range, length(ic_values))
    mthd = GA(populationSize = population_size, selection = tournament(Int(population_size/10)),
            crossover = TPX, crossoverRate = 0.5,
            mutation  = BGA(valrange, 2), mutationRate = 0.7)

    # Progress Bar
    progress = ProgressThresh(opts.abstol, "Converging: ")

    # Optimization
    result = Evolutionary.optimize(fitness_function, myconstraints, mthd, pop, opts)

    # Get the individual, fitness, and extradata of the population
    record = reduce(vcat,[gen.metadata["staterecord"] for gen in result.trace])

    # Filter out individuals with fitness values less than 0.1
    # fitpops = filter(x -> x.fit < -0.1, record)

    # # Get the best solution
    # newp = result.minimizer
    # newsol = solve(remake(fullprob, u0=[newp,zeros(11)]))


    # #plot the results
    # plot(newsol, xlabel="Time (s)", ylabel="Concentration (mM)", title="Optimized Model")
end

testGAfunc(eval_ic_fitness, fitness_function, fullprob)


ic_constraints = define_initialcondition_constraints()

ga_problem = GAProblem(ic_constraints, fullprob)

run_GA(ga_problem; population_size=10000, iterations = 5)


