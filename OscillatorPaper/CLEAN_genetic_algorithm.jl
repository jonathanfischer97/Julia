using Catalyst
using DifferentialEquations
using Plots
using Evolutionary
using StaticArrays
using BenchmarkTools
using FFTW
# using Zygote
using Statistics
using Random
using Distributions
using Peaks

Threads.nthreads()




## COST FUNCTIONS and dependencies 
function getDif(indexes::Vector{Int}, arrayData::Vector{Float64})::Float64 #get difference between peaks 
    arrLen = length(indexes)
    sum = 0
    for (i, ind) in enumerate(indexes)
        if i == arrLen
            break 
        end
        sum += arrayData[ind] - arrayData[indexes[i+1]]
    end
    sum += arrayData[indexes[end]] 
    return sum
end
    
function getSTD(indexes::Vector{Int}, arrayData::Vector{Float64}, window::Int)::Float64 #get standard deviation of fft peak indexes
    numPeaks = length(indexes)
    arrLen = length(arrayData)
    sum = 0
    for ind in indexes 
        minInd = max(1, ind - window)
        maxInd = min(arrLen, ind+window)
        sum += std(arrayData[minInd:maxInd])
    end
    sum = sum/numPeaks 
    return sum
end

 
function getFrequencies(y::Vector{Float64}) #normally would slice y by interval for a sampling rate but not here
    #fft sample rate: 1 sample per 5 minutes
    res = broadcast(abs,rfft(y)) #broadcast abs to all elements of rfft return array. Think .= syntax does the same 
    #normalize the amplitudes
    res/cld(10000, 2)
end




# parse string input into vector of floats
function parse_input(input::AbstractString)::Vector{Float64}
    output = Float64[]
    input = replace(input, r"[\{\}]"=> "") # remove curly braces
    items = split(input, ", ")
    for s in items
        val = split(s, ": ")[2]
        push!(output, parse(Float64, val))
    end
    return output
end

## PROBLEM SETUP
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

    #previously hidden NERDSS reactions
    (ka2,kb2), Lp + AK <--> LpAK
    (ka2*y,kb2), Lp + AKL <--> LpAKL
    (ka2,kb2), Lp + AP <--> LpAP
    (ka2*y,kb2), Lp + APLp <--> LpAPLp
    (ka3,kb3), A + K <--> AK
    (ka4,kb4), A + P <--> AP
    (ka3,kb3), A + LK <--> AKL
    (ka4,kb4), A + LpP <--> APLp
    (ka3*y,kb3), LpA + LK <--> LpAKL
    (ka4*y,kb4), LpA + LpP <--> LpAPLp
    (ka1,kb1), AK + L <--> AKL #binding of kinase to lipid
    kcat1, AKL --> Lp + AK #phosphorylation of lipid
    (ka7,kb7), AP + Lp <--> APLp #binding of phosphatase to lipid
    kcat7, APLp --> L + AP #dephosphorylation of lipid
end ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 ka7 kb7 kcat7 y 

input = "{'ka1': 0.01642097502199415, 'kb1': 0.0643048008980449, 'kcat1': 286.6382253193995, 'ka2': 1.0, 'kb2': 0.39569337786534897, 'ka3': 0.024784025572903687, 'kb3': 0.5393197910059361, 'ka4': 0.03281183067630924, 'kb4': 0.2897657637531564, 'ka7': 0.11450246770405478, 'kb7': 0.0028126177315505618, 'kcat7': 1.2733781341040291, 'VA': 2650.5049775102034, 'L': 2.107360215531841, 'Lp': 2.624942515606095, 'K': 0.15746804082956445, 'P': 0.7958264308210892, 'A': 2.785231513223431}"
parsed = parse_input(input)

p = parsed[1:end-5]
u0 = [:L => parsed[end-4], :Lp => parsed[end-3], :K => parsed[end-2], :P => parsed[end-1], :A => parsed[end], :LK => 0., :LpA => 0., :LpAK => 0., :LpAP => 0., :LpAPLp => 0., :LpAKL => 0., :LpP => 0., :AK => 0., :AP => 0., :AKL => 0., :APLp => 0.]


tspan = (0.,100.)

oprob = ODEProblem(osc_rn, u0, tspan, p)
osol = solve(oprob, saveat = 0.01)
plot(osol)











ka_min, ka_max = 0.0001, 1.0
kb_min, kb_max = 0.0001, 1000.0
kcat_min, kcat_max = 0.0001, 1000.0

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
    "y" => Dict("min" => 100., "max" => 10000.)
)




















### SANITY CHECK with simple cost function### 
function test_cost_function(p::Vector{Float64}, prob::ODEProblem, target_ss::Vector{Float64})
    Y = solve(remake(prob, p=p))
    ss_values = Y[1, end] # Get the steady-state values of the ODEProblem solution
    cost = sum((ss_values .- target_ss).^2) # Calculate the sum of squared differences
    return cost
end

function make_test_fitness_function(prob::ODEProblem, target_ss::Vector{Float64})
    function fitness_function(p::Vector{Float64})
        return test_cost_function(p, prob, target_ss)
    end
    return fitness_function
end

target_ss = [2.107] # Replace with your desired target steady-state values
test_fitness_function = make_test_fitness_function(oprob, target_ss)



constraints = BoxConstraints([param_values[p]["min"] for p in keys(param_values)], [param_values[p]["max"] for p in keys(param_values)])
opts = Evolutionary.Options(show_trace=true,show_every=1, store_trace=true, iterations=20, parallelization=:thread)

common_range = 0.5
valrange = fill(common_range, 13)

mthd = GA(populationSize = 5000, selection = tournament(500;select=argmax),
crossover = TPX, crossoverRate = 0.5,
mutation  = BGA(valrange, 2), mutationRate = 0.75)

result = Evolutionary.optimize(test_fitness_function, constraints, mthd, opts)

newp = result.minimizer
newsol = solve(remake(oprob2, p=newp))
newsol[1, end] # Get the steady-state values of the ODEProblem solution






## REAL FITNESS FUNCTION ##
function eval_fitness(p::Vector{Float64}, prob::ODEProblem)
    Y = solve(remake(prob, p=p), saveat = 0.01)
    p1 = Y[1,:]
    fftData = getFrequencies(p1) #get Fourier transform of p1
    indexes = findmaxima(fftData)[1] #get indexes of local maxima of fft
    if length(indexes) == 0 #if no peaks, return 0
        return 0.
    end
    std = getSTD(indexes, fftData, 1) #get standard deviation of fft peak indexes
    diff = getDif(indexes, fftData) #get difference between peaks
    return std + diff 
end

function make_fitness_function(prob::ODEProblem)
    function fitness_function(p::Vector{Float64})
        return eval_fitness(p, prob)  
    end
    return fitness_function
end

fitness_function = make_fitness_function(oprob) # Create a fitness function that includes your ODE problem as a constant




#Real optimization
constraints = BoxConstraints([param_values[p]["min"] for p in keys(param_values)], [param_values[p]["max"] for p in keys(param_values)])
opts = Evolutionary.Options(show_trace=true,show_every=1, store_trace=true, iterations=10, parallelization=:threads)

common_range = 0.7
valrange = fill(common_range, 13)
mthd = GA(populationSize = 10000, selection = tournament(1000;select=argmax),
          crossover = TPX, crossoverRate = 0.5,
          mutation  = BGA(valrange, 2), mutationRate = 0.75)
result = Evolutionary.optimize(fitness_function, p, mthd, opts)



trace_vals = [x.value for x in result.trace]
trace_time = [x.iteration for x in result.trace]

plot(trace_time,trace_vals, xlabel="Iteration", ylabel="Cost", label="Cost", title="Cost vs Iteration", legend=:topleft)