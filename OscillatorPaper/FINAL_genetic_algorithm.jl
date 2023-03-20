using Catalyst
using DifferentialEquations
using Plots
using Evolutionary
# using StaticArrays
using BenchmarkTools
using FFTW
# using Zygote
using Statistics
using Random
using Distributions
using Peaks
using LinearAlgebra





## COST FUNCTIONS and dependencies 
function getDif(indexes::Vector{Int}, arrayData::Vector{Float64})::Float64 #get difference between fft peak indexes
    arrLen = length(indexes)
    if arrLen < 2
        return 0.0
    end
    sum_diff = @inbounds sum(arrayData[indexes[i]] - arrayData[indexes[i+1]] for i in 1:(arrLen-1))
    sum_diff += arrayData[indexes[end]]
    return sum_diff
end


function getSTD(indexes::Vector{Int}, arrayData::Vector{Float64}, window::Int)::Float64 #get standard deviation of fft peak indexes
    sum_std = @inbounds sum(std(arrayData[max(1, ind - window):min(length(arrayData), ind + window)]) for ind in indexes)
    return sum_std / length(indexes)
end


 
function getFrequencies(y::Vector{Float64}, tlen::Int)::Vector{Float64} #get fft of ODE solution
    res = broadcast(abs,rfft(y)) #broadcast abs to all elements of rfft return array. 
    res./cld(tlen, 2) #normalize the amplitudes
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


input = "{'ka1': 3.6396514622137452, 'kb1': 10.311889380044056, 'kcat1': 361.5709556626289, 'ka2': 89.4797815238423, 'kb2': 500.0, 'ka3': 98.3482118725635, 'kb3': 500.0, 'ka4': 9.746019998479238, 'kb4': 121.60536126369898, 'ka7': 72.19871144097954, 'kb7': 124.66417463084926, 'kcat7': 80.47863850796283, 'VA': 1500}"

input = "{'ka1': 0.01642097502199415, 'kb1': 0.0643048008980449, 'kcat1': 286.6382253193995, 'ka2': 1.0, 'kb2': 0.39569337786534897, 'ka3': 0.024784025572903687, 'kb3': 0.5393197910059361, 'ka4': 0.03281183067630924, 'kb4': 0.2897657637531564, 'ka7': 0.11450246770405478, 'kb7': 0.0028126177315505618, 'kcat7': 1.2733781341040291, 'VA': 2650.5049775102034, 'L': 2.107360215531841, 'Lp': 2.624942515606095, 'K': 0.15746804082956445, 'P': 0.7958264308210892, 'A': 2.785231513223431}"
parsed = parse_input(input)

p = parsed
p = parsed[1:end-5]
u0 = [:L => parsed[end-4], :Lp => parsed[end-3], :K => parsed[end-2], :P => parsed[end-1], :A => parsed[end], :LK => 0., :LpA => 0., :LpAK => 0., :LpAP => 0., :LpAPLp => 0., :LpAKL => 0., :LpP => 0., :AK => 0., :AP => 0., :AKL => 0., :APLp => 0.]
# u0 = [:L => 0., :Lp => 3.0, :K => 0.2, :P => 0.3, :LK => 0., :A => 0.5, :LpA => 0., :LpAK => 0., :LpAP => 0., :LpAPLp => 0., :LpAKL => 0., :LpP => 0., :AK => 0., :AP => 0.]


tspan = (0.,100.)


u0 = [0.0, 0.3, 0.0, 3.0, 0.9, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0]
oprob = ODEProblem(osc_rn, u0, tspan, p)
osol = solve(oprob)
plot(osol)



p2 = p./2
osol2 = solve(remake(oprob, p=p2))
plot(osol2)


osolshort = solve(remake(oprob, tspan=(0.,10.)))
osollong = solve(remake(oprob, tspan=(0.,100.)))





ka_min, ka_max = 0.0, 1.
kb_min, kb_max = 0.00, 100.0
kcat_min, kcat_max = 0.00, 500.0


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
)







## REAL FITNESS FUNCTION ##
function eval_fitness(p::Vector{Float64}, prob::ODEProblem)::Float64
    Y = solve(remake(prob, p=p), saveat=0.01, save_idxs = 5)
    # p1 = Y[5,:] #get first species
    fftData = getFrequencies(Y.u, length(Y.t)) #get Fourier transform of p1
    indexes = findmaxima(fftData)[1] #get indexes of local maxima of fft
    realpeaks = length(findmaxima(Y.u, 10)[1]) #get indexes of local maxima of p1
    if length(indexes) == 0 || realpeaks < 4  #if no peaks, return 0
        return 0.
    end
    std = getSTD(indexes, fftData, 1) #get standard deviation of fft peak indexes
    println(std)
    diff = getDif(indexes, fftData) #get difference between peaks
    println(diff)
    return - (std*10) - diff 
end

# NEW EVALUATION FUNCTION
function normalize_time_series(y::Vector{Float64})::Vector{Float64}
    return y / maximum(abs.(y))
end

function autocorrelation(y::Vector{Float64})::Vector{Float64}
    n = length(y)
    y_mean = mean(y)
    y_centered = y .- y_mean
    autocorr = [dot(y_centered[1:(n - k)], y_centered[(k + 1):n]) for k in 0:(n - 1)]
    return autocorr ./ autocorr[1]
end

function peak_count(y::Vector{Float64}, threshold::Float64)::Int
    abs_threshold = threshold * maximum(y)
    peaks = findmaxima(y)
    return length(filter(x -> x > abs_threshold, peaks[2]))
end

function eval_fitness2(p::Vector{Float64}, prob::ODEProblem, threshold::Float64 = 0.1)::Float64
    Y = solve(remake(prob, p=p), saveat=0.01, save_idxs = 5)

    Y_normalized = normalize_time_series(Y.u)
    autocorr = autocorrelation(Y_normalized)
    println(autocorr)
    fitness = peak_count(autocorr, threshold)
    return -fitness
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
opts = Evolutionary.Options(show_trace=true,show_every=1, store_trace=true, iterations=20, parallelization=:thread)

common_range = 0.5
valrange = fill(common_range, 13)
mthd = GA(populationSize = 5000, selection = tournament(500;select=argmin),
          crossover = TPX, crossoverRate = 0.5,
          mutation  = BGA(valrange, 2), mutationRate = 0.9, ɛ = 0.1)
result = Evolutionary.optimize(fitness_function, constraints, mthd, opts)

newp = result.minimizer
newsol = solve(remake(oprob, p=newp))
# plot(newsol, vars=(1,2), title="ODE Solution", xlabel="Time", ylabel="Concentration", label="Concentration")
plot(newsol)

newsol1 = newsol[1,:]
newfft = abs.(rfft(newsol1))
optfft = abs.(rfft(osol[1,:]))

#subplots with time plot and fft plot
default(size=(800,600), dpi=300)
p1 = plot(newsol.t, newsol1, title="ODE Solution", xlabel="Time", ylabel="Concentration", label = "New")
plot!(osol.t, osol[1,:], label="Optimal")
p2 = plot(newfft, title="Fourier Transform", xlabel="Frequency", ylabel="Amplitude", xlim = (0.,20.), label="New")
plot!(p2, optfft, label="Optimal")
plot(p1, p2, layout=(2,1))


plot(newfft, title="Fourier Transform", xlabel="Frequency", ylabel="Amplitude", label="New")
plot!(optfft, label="Optimal")




results.trace.entries[end].population


eval_fitness(result.minimizer,oprob)

trace_vals = [x.value for x in result.trace]
trace_time = [x.iteration for x in result.trace]

plot(trace_time,trace_vals, xlabel="Iteration", ylabel="Cost", label="Cost", title="Cost vs Iteration", legend=:topleft)