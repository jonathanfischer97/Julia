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
# include("../../UTILITIES/EvaluationFunctions.jl")

# import the genetic algorithm and associated functions
# include("../../UTILITIES/GA_functions.jl")

#! Solve model for arbitrary oscillatory parameters and initial conditions
begin
    #? Parameter list
    psym = [:ka1 => 0.009433439939827041, :kb1 => 2.3550169939427845, :kcat1 => 832.7213093872278, :ka2 => 12.993995997539924, :kb2 => 6.150972501791291,
            :ka3 => 1.3481451097940793, :kb3 => 0.006201726090609513, :ka4 => 0.006277294665474662, :kb4 => 0.9250191811994848, :ka7 => 57.36471615394549, 
            :kb7 => 0.04411989797898752, :kcat7 => 42.288085868394326, :DF => 3631.050539219606]
    p = [x[2] for x in psym]
        
    #? Initial condition list
    usym = [:L => 0.0, :Lp => 3.0, :K => 0.5, :P => 0.3, :A => 2.0, :LpA => 0.0, :LK => 0.0, 
            :LpP => 0.0, :LpAK => 0.0, :LpAP => 0.0, :LpAKL => 0.0, :LpAPLp => 0.0, :AK => 0.0, :AP => 0.0, 
            :AKL => 0.0, :APLp => 0.0]
    u0 = [x[2] for x in usym]

    #? Timespan for integration
    tspan = (0., 100.)

    #? Create ODE problem and solve
    fullprob = ODEProblem(fullrn, u0, tspan, p)
    # sol = solve(prob, saveat=0.1, save_idxs=1)

    #? Plot the results
    # plot(sol)
end

#! Helper functions for cost function ## 
begin
    """Get summed difference of peaks in the frequency domain"""
    function getDif(indexes::Vector{Int}, arrayData::Vector{Float64}) 
        idxarrLen = length(indexes)
        if idxarrLen < 2
            return 0.0
        end
        sum_diff = @inbounds sum(arrayData[indexes[i]] - arrayData[indexes[i+1]] for i in 1:(idxarrLen-1))
        sum_diff += arrayData[indexes[end]]
        return sum_diff
    end

    """Get summed average standard deviation of peaks in the frequency domain"""
    function getSTD(peakindxs::Vector{Int}, arrayData::Vector{Float64}, window_ratio::Float64) #get average standard deviation of fft peak indexes
        arrLen = length(arrayData)
        window = max(1, round(Int, window_ratio * arrLen))
        sum_std = @inbounds sum(std(arrayData[max(1, ind - window):min(arrLen, ind + window)]) for ind in peakindxs)
        return sum_std / length(peakindxs)
    end

    """Return normalized FFT of solution vector"""
    function getFrequencies(y::Vector{Float64})
        res = abs.(rfft(y))
        return res ./ cld(length(y), 2) #normalize amplitudes
    end

    """Calculates the period and amplitude of each individual in the population"""
    function getPerAmp(sol::ODESolution)
        # Find peaks and calculate amplitudes and periods
        indx_max, vals_max = findmaxima(sol.u, 1)
        indx_min, vals_min = findminima(sol.u, 1)

        if length(indx_max) < 2 || length(indx_min) < 2
            return 0., 0.
        else
            # Calculate amplitudes and periods
            @inbounds amps = [(vals_max[i] - vals_min[i])/2 for i in 1:min(length(indx_max), length(indx_min))]
            @inbounds pers = [sol.t[indx_max[i+1]] - sol.t[indx_max[i]] for i in 1:(length(indx_max)-1)]

            # Calculate means of amplitudes and periods
            amp = mean(amps)
            per = mean(pers)

            return per, amp
        end
    end

    """Cost function to be plugged into eval_fitness wrapper"""
    function CostFunction(Y::ODESolution)
        #get the fft of the solution
        fftData = getFrequencies(Y.u)
        fftindexes = findmaxima(fftData,1)[1] #get the indexes of the peaks in the fft
        timeindexes = findmaxima(Y.u,10)[1] #get the times of the peaks in the fft
        if isempty(fftindexes) || length(timeindexes) < 2 #if there are no peaks, return 0
            return 0.0
        end
        std = getSTD(fftindexes, fftData, 0.0001) #get the standard deviation of the peaks
        diff = getDif(fftindexes, fftData) #get the difference between the peaks


        # Return cost, period, and amplitude as a tuple
        return -std + diff
    end
end


"""Cost function that catches errors and returns 1.0 if there is an error"""
function eval_fitness_catcherrors(p::Vector{Float64},  prob::ODEProblem)
    Y = nothing
    try 
        Y = solve(remake(prob, view(u0,1:5)=p), saveat=0.1, save_idxs=1, maxiters=100000, verbose=false)
        if Y.retcode in (ReturnCode.Unstable, ReturnCode.MaxIters) #|| any(x==1 for array in isnan.(Y) for x in array) || any(x==1 for array in isless.(Y, 0.0) for x in array)
            return 1.0
        end
    catch e 
        if e isa DomainError #catch domain errors
            return 1.0
        else
            rethrow(e) #rethrow other errors
        end
    end
    fitness = CostFunction(Y)
    return -fitness
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

# Create a PeriodAmplitudes instance
# tracker = PeriodAmplitudes()

fitness_function = make_fitness_function(eval_fitness_catcherrors, fullprob) # Create a fitness function that includes your ODE problem as a constant

# using Debugger
#! Optimization block
begin
    population_size = 1000
    pop = generate_population(ic_values, population_size)

    myconstraints = BoxConstraints([ic_values[p]["min"] for p in keys(ic_values)], [ic_values[p]["max"] for p in keys(ic_values)])
    opts = Evolutionary.Options(abstol=1e-2, reltol=1.00, successive_f_tol = 5, iterations=10, store_trace = true, 
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
    fitpops = filter(x -> x.fit < -0.1, record)

    # Get the best solution
    newp = result.minimizer
    newsol = solve(remake(fullprob, u0=vcat(newp,zeros(11))))


    #plot the results
    plot(newsol, xlabel="Time (s)", ylabel="Concentration (mM)", title="Optimized Model")
end