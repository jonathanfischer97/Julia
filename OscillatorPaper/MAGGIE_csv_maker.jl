using Plots 
using Catalyst
using DifferentialEquations
using Peaks
using Statistics
using BenchmarkTools
using CSV
using DataFrames


"""Full oscillator model for comparison"""
osc_rn = @reaction_network osc_rn begin
    (ka1,kb1), L + K <--> LK
    kcat1, LK --> Lp + K 
    (ka2,kb2), Lp + A <--> LpA 
    (ka3,kb3), LpA + K <--> LpAK  
    (ka1*y,kb1), LpAK + L <--> LpAKL
    kcat1, LpAKL --> Lp + LpAK  
    (ka7,kb7), Lp + P <--> LpP 
    kcat7, LpP --> L + P
    (ka4,kb4), LpA + P <--> LpAP 
    (ka7*y,kb7), Lp + LpAP <--> LpAPLp
    kcat7, LpAPLp --> L + LpAP
end

    
#parameter list
"""ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, y"""
p = [0.05485309578515125, 19.774627209108715, 240.99536193310848, 
    1.0, 0.9504699043910143, 
    41.04322510426121, 192.86642772763489,
    0.19184180144850807, 0.12960624157489123, 
    0.6179131289475834, 3.3890271820244195, 4.622923709012232,750]

#initial condition list
u0 = [0., 0.3, 0.0, 3.0, 0.6, 0.0, 0., 0., 0., 0., 0., 0.] #working 
u0 = [0.01, 3.0, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.6, 0.3, 0.2]
L, Lp, LpA, LpAP, LK, LpAK, LpAKL, LpP, LpAPLp, A, K, P = u0

umap = [:L => L, :Lp => Lp, :LpA => LpA, :LpAP => LpAP, :LK => LK, :LpAK => LpAK, :LpAKL => LpAKL, :LpP => LpP, :LpAPLp => LpAPLp, :A => A, :K => K, :P => P]



#timespan for integration
tspan = (0., 100.);


#solve the full ODEs

prob = ODEProblem(osc_rn, u0, tspan, p)
sol = solve(prob)

#plot the results
plot(sol, lw=2, title="Full Oscillator Model", xlabel="Time", ylabel="Concentration")




## Genetic algorithm to find the best parameters for the reduced model ## 
using Evolutionary
using FFTW

## COST FUNCTIONS and dependencies 
function getDif(indexes::Vector{Int}, arrayData::Vector{Float64}) #get summed difference between fft peak indexes
    arrLen = length(indexes)
    if arrLen < 2
        return 0.0
    end
    sum_diff = @inbounds sum(arrayData[indexes[i]] - arrayData[indexes[i+1]] for i in 1:(arrLen-1))
    sum_diff += arrayData[indexes[end]] #add the last element
    return sum_diff
end

function getSTD(peakindxs::Vector{Int}, arrayData::Vector{Float64}, window_ratio::Float64) #get average standard deviation of fft peak indexes
    window = max(1, round(Int, window_ratio * length(arrayData)))
    sum_std = @inbounds sum(std(arrayData[max(1, ind - window):min(length(arrayData), ind + window)]) for ind in peakindxs)
    return sum_std / length(peakindxs)
end


function getFrequencies(y::Vector{Float64})
    res = abs.(rfft(y))
    return res ./ cld(length(y), 2) #normalize amplitudes
end



function eval_fitness(p::Vector{Float64},  prob::ODEProblem, idxs; dt = 0.1)
    Y = solve(remake(prob, p=p), Tsit5(), save_idxs=idxs, saveat=dt)
    # if any(isnan.(Y.u)) || any(isless.(Y.u, 0.0))
    #     return 0.0
    # end
    fftData = getFrequencies(Y.u)
    indexes = findmaxima(fftData, 1)[1]
    if isempty(indexes)
        return 0.
    end
    std = getSTD(indexes, fftData, 0.1)
    diff = getDif(indexes, fftData)
    return -std - diff
end


#function closer for fitness function
function make_fitness_function(prob::ODEProblem; idxs = 1) # Create a fitness function that includes your ODE problem as a constant
    function fitness_function(p::Vector{Float64})
        return eval_fitness(p, prob, idxs)  
    end
    return fitness_function
end

fitness_function = make_fitness_function(prob) # Create a fitness function that includes your ODE problem as a constant



## parameter constraint ranges ##
ka_min, ka_max = 0.0001, 100.
kb_min, kb_max = 0.0001, 500.0
kcat_min, kcat_max = 0.0001, 500.0


param_values = Dict(
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
    "y" => Dict("min" => 100., "max" => 5000.)
);



#Optimization parameters
constraints = BoxConstraints([param_values[p]["min"] for p in keys(param_values)], [param_values[p]["max"] for p in keys(param_values)])

function loguniform_population(n, param_values)
    lower_bounds = [param_values[p]["min"] for p in keys(param_values)]
    upper_bounds = [param_values[p]["max"] for p in keys(param_values)]
    population = Vector{Vector{Float64}}(undef, n)
    for i in 1:n
        population[i] = [10^(rand() * (log10(upper_bounds[j]) - log10(lower_bounds[j])) + log10(lower_bounds[j])) for j in eachindex(lower_bounds)]
    end
    return population
end


initial_population = loguniform_population(5000, param_values)



opts = Evolutionary.Options(show_trace=true,show_every=1, store_trace=true, iterations=20, parallelization=:thread)

common_range = 0.5
valrange = fill(common_range, 13)
mthd = GA(populationSize = 5000, selection = tournament(500),
          crossover = TPX, crossoverRate = 0.5,
          mutation  = BGA(valrange, 2), mutationRate = 0.9)


# function callback(trace)

#Optimization
result = Evolutionary.optimize(fitness_function, constraints, mthd, opts)

newp = result.minimizer
newsol = solve(remake(prob, p=newp), saveat=0.1)
plot(newsol)
eval_fitness(newp, prob, 1)

trace = result.trace
sorted_trace = sort(trace, by = x -> x.value, rev = true) # Sort the trace by fitness value in descending order
top_10_fittest = sorted_trace[1:10] # Extract the top 10 fittest parameter values

parameter_names = collect(keys(param_values))
top_10_fittest_table = [transpose([p.minimizer for p in top_10_fittest])...] # Create a table with the top 10 fittest parameter values
CSV.write("top_10_fittest_parameters.csv", DataFrame(top_10_fittest_table, Symbol.(parameter_names))) # Write the table to a CSV file





