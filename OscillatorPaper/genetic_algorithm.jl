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

## made by ChatGPT
function getSTD(indexes::Vector{Int}, arrayData::Vector{Float64})
    numPeaks = length(indexes)
    arrLen = length(arrayData)
    sum = 0
    for ind in indexes
        # Determine the frequency corresponding to the current peak
        freq = (ind - 1) / arrLen
        
        # Compute the window size as a function of the frequency
        windowSize = 0.1 / freq  # Scale the window size by the reciprocal of the frequency
        
        # Define the bounds of the window
        minInd = max(1, Int(ind - windowSize))
        maxInd = min(arrLen, Int(ind + windowSize))
        
        # Compute the standard deviation over the window
        sum += std(arrayData[minInd:maxInd])
    end
    return sum / numPeaks
end

function getSamplingRate(max_peaks::Int, time_interval::Float64)
    max_frequency = max_peaks / time_interval
    sampling_rate = 2 * max_frequency
    return sampling_rate
end

 
function getFrequencies(y::Vector{Float64}) #normally would slice y by interval for a sampling rate but not here
    #fft sample rate: 1 sample per 5 minutes
    res = broadcast(abs,rfft(y)) #broadcast abs to all elements of rfft return array. Think .= syntax does the same 
    #normalize the amplitudes
    res/cld(10000, 2)
end


#homemade peakfinder
function findlocalmaxima(signal::Vector) 
    inds = Int[]
    # buff = Zygote.Buffer(inds) #Zygote.Buffer creates mutable array that is invisible to autodifferention, aka no matrix operations allowed
    if length(signal)>1
        if signal[1]>signal[2]
            push!(inds,1)
        end
        for i=2:length(signal)-1
            if signal[i-1]<signal[i]>signal[i+1]
                push!(inds,i)
            end
        end
        if signal[end]>signal[end-1]
            push!(inds,length(signal))
        end
    end
    return inds #return copy of mutable buffer array to make immutable (essential for Zygote pullback)
  end

using StatsBase

function find_peak_indexes(data::AbstractArray{T}, thres::Union{Real, Nothing}=0.5, min_dist::Int=1) where T<:Real
    length(data) < 3 && return Int[] # Return an empty array for data with length < 3
    thres = isnothing(thres) ? 0.0 : thres # If threshold is not provided, set it to 0.0
    dx = 1
    y = data .- thres
    dy = diff(y) |> Vector{T} # Ensure that dy is always a 1-dimensional vector
    min_peak_distance = max(1, min_dist)
    zero_crossings = findall(dy .> 0) .+ 1
    zero_crossings_diff = diff(zero_crossings)
    candidate_peaks = zero_crossings[zero_crossings_diff >= min_peak_distance]
    if isempty(candidate_peaks)
        return Int[]
    end
    peak_properties = [
        (p - 1 + argmax(y[p:min(p+dx, end)]), y[p - 1 + argmax(y[p:min(p+dx, end)])])
        for p in candidate_peaks
    ]
    peak_properties = filter(p -> p[1] > 0, peak_properties)
    return sort([p[1] for p in peak_properties])
end







function getCorr(fftData::Vector{Float64})
    n = length(fftData)
    x = 1:n
    y = fftData
    cosine = cos.(2 * π * x / n)
    mean_y = mean(y)
    mean_cos = mean(cosine)
    std_y = std(y)
    std_cos = std(cosine)
    corr = sum((y .- mean_y) .* (cosine .- mean_cos)) / ((n-1) * std_y * std_cos)
    return corr
end

using Statistics

function findpeaks(x::Vector{T}, prominence::Union{Nothing, T}=nothing, width::Union{Nothing, T}=nothing) where T<:Real
    # Smooth the data using a Gaussian filter
    window_size = Int(4 * ceil(prominence / 4) + 1)
    window = gausswin(window_size, prominence / 4)
    y = conv(x, window, mode="same")

    # Find local maxima in the smoothed data
    maxima = findall(2: length(y)-1) do i
        y[i-1] < y[i] >= y[i+1]
    end

    # Filter the maxima by prominence and width
    if prominence != nothing
        peak_values = x[maxima]
        peak_prominence = peak_values .- maximum([x[max(i-1,1):min(i+1,length(x))] for i in maxima])
        maxima = maxima[findall(p -> p >= prominence, peak_prominence)]
    end

    if width != nothing
        halfwidth = div(width, 2)
        maxima = filter(i -> i - halfwidth > 0 && i + halfwidth <= length(x), maxima)
        maxima = filter(i -> i - halfwidth < 2 || i + halfwidth > length(x)-1 || maximum(x[i-halfwidth:i+halfwidth]) >= x[i], maxima)
    end

    return maxima
end

function gausswin(n::Int, alpha::Float64=2.5)
    if n < 1
        error("n must be greater than 0")
    end
    if alpha < 0
        error("alpha must be non-negative")
    end
    if n == 1
        return [1.0]
    end
    half = div(n, 2)
    x = collect(-half:half)
    w = exp.(-0.5 * (x / (alpha * half)).^2)
    return w ./ sum(w)
end


#current testing cost function, USE 
function get_fitness(p::Vector{Float64})::Float64
    Y = solve(remake(oprob2, p=p),saveat=0.001)
    p1 = Y[1,:]
    # println(length(p1))
    fftData = getFrequencies(p1) #get Fourier transform of p1
    indexes = findmaxima(fftData)[1] #get indexes of local maxima of fft
    # println(indexes)
    # println("Peak Number: ", findmaxima(p1, 400)[1])
    if length(indexes) == 0 #if no peaks, return 0
        return 0.
    end
    std = getSTD(indexes, fftData, 1) #get standard deviation of fft peak indexes
    diff = getDif(indexes, fftData) #get difference between peaks
    return std + diff
end

function testcost2(p::Vector{Float64})
    Y = Array(solve(oprob, tspan = (0., 20.), p = p, saveat = 0.001))
    p1 = Y[1,:]
    fftData = getFrequencies1(p1) #get Fourier transform of p1
    indexes = findPeaks(fftData, 0.05) #get indexes of local maxima of fft
    if length(indexes) == 0 #if no peaks, return 0
        return 0
    end
    std = corrPeakStd(indexes, fftData, 1) #get standard deviation of fft peak indexes
    diff = corrPeakDiff(indexes, fftData) #get difference between peaks
    return std + diff
end

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

oprob2 = ODEProblem(osc_rn, u0, tspan, p)
osol2 = solve(oprob2, saveat = 0.01)
plot(osol2)
get_fitness(p)



signal = osol[1,:]

pks,vals = findmaxima(signal,3)
for idx in pks
    println("Peak at time $(osol.t[idx])")
end

getFrequencies1(signal)















uSA = SA[0., 0.2, 0., 3.0, 3.0, 0., 0., 0., 0.3, 0., 0., 0., 0., 0., 0., 0.]
osys = convert(ODESystem, osc_rn)
pmap = symmap_to_varmap(osc_rn, p)
oprob = ODEProblem{false}(osys, uSA, tspan, pmap)

@benchmark osol = solve(oprob, saveat = 0.01)


plot(osol)


# generate candidate parameter vector
function generateCandidate(param_values::Dict{String, Dict{String}})
    candidate = Float64[]
    for p in param_values
        param_range = p[2]
        push!(candidate, rand(Uniform(param_range["min"], param_range["max"])))
    end
    return candidate
end




ka_min, ka_max = 0.0, 100.0
kb_min, kb_max = 0.0, 500.0
kcat_min, kcat_max = 0.0, 500.0

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


constraints = BoxConstraints([param_values[p]["min"] for p in keys(param_values)], [param_values[p]["max"] for p in keys(param_values)])
opts = Evolutionary.Options(show_trace=true,show_every=1, store_trace=true, iterations=20, parallelization=:thread)
mthd = GA(populationSize = 5000, selection = tournament(2;select=argmax),
          crossover = TPX, crossoverRate = 0.5,
          mutation  = uniform(), mutationRate = 0.7)

result = Evolutionary.optimize(get_fitness, constraints, mthd, opts)



opts = Evolutionary.Options(show_trace=true,show_every=10, store_trace=true, iterations=10, parallelization=:thread)
mthd = GA(populationSize = 5000, selection = roulette,
crossover = TPX, crossoverRate = 0.5,
mutation  = PLM(), mutationRate = 0.9)

result = Evolutionary.optimize(fitness_function, constraints, mthd, opts)

newp = result.minimizer
newsol = solve(remake(oprob2, p=newp))
plot(newsol)













result = Evolutionary.optimize(testcost, constraints, GA(populationSize = Int(5000), crossoverRate = 0.5, crossover = TPX, mutationRate = 0.90, mutation = inversion,  metrics = [Evolutionary.AbsDiff(1e-2)]))
newsol = solve(remake(oprob, p = result.minimizer), Tsit5(), saveat = 0.001)
plot(newsol)


newp= [1.5088004628136753, 474.74790524323765, 413.77725974011156, 53.161947464834014, 450.48565960154116, 97.63166395821786, 343.31846212387194, 4.858456659722766, 129.45559308306565, 96.4450517730593, 67.91459290725088, 125.82747629586815, 1.0446379594697963/0.001]
newp = [:ka1 => newp[1], :kb1 => newp[2], :kcat1 => newp[3], 
:ka2 => newp[4], :kb2 => newp[5], 
:ka3 => newp[6], :kb3 => newp[7],
:ka4 => newp[8], :kb4 => newp[9], 
:ka7 => newp[10], :kb7 => newp[11], :kcat7 => newp[12], :y => newp[13]]
newp = symmap_to_varmap(osc_rn, newp)
oprob = ODEProblem{false}(osys, u0, tspan, newp)
osol = solve(oprob, saveat = 0.001)
plot(osol)

newsol = solve(remake(oprob, p = newp), Tsit5(), saveat = 0.001)
plot(newsol)







function eval_fitness(p::Vector{Float64}, oprob::ODEProblem)
    Y = solve(remake(oprob, p=p), saveat = 0.01)
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

function make_fitness_function(oprob::ODEProblem)
    function fitness_function(p::Vector{Float64})
        return eval_fitness(p, oprob)  
    end
    return fitness_function
end

fitness_function = make_fitness_function(oprob2) # Create a fitness function that includes your ODE problem as a constant


### SANITY CHECK ### 
function test_cost_function(p::Vector{Float64}, oprob::ODEProblem, target_ss::Vector{Float64})
    Y = solve(remake(oprob, p=p))
    ss_values = Y[1, end] # Get the steady-state values of the ODEProblem solution
    cost = sum((ss_values .- target_ss).^2) # Calculate the sum of squared differences
    return cost
end

function make_test_fitness_function(oprob::ODEProblem, target_ss::Vector{Float64})
    function fitness_function(p::Vector{Float64})
        return test_cost_function(p, oprob, target_ss)
    end
    return fitness_function
end







constraints = BoxConstraints([param_values[p]["min"] for p in keys(param_values)], [param_values[p]["max"] for p in keys(param_values)])
opts = Evolutionary.Options(show_trace=true,show_every=1, store_trace=true, iterations=20, parallelization=:thread)

common_range = 0.5
valrange = fill(common_range, 13)

mthd = GA(populationSize = 5000, selection = tournament(500;select=argmax),
          crossover = TPX, crossoverRate = 0.5,
          mutation  = BGA(valrange, 2), mutationRate = 0.75)

target_ss = [2.107] # Replace with your desired target steady-state values
test_fitness_function = make_test_fitness_function(oprob2, target_ss)
result = Evolutionary.optimize(test_fitness_function, constraints, mthd, opts)

#Real fitness function
result = Evolutionary.optimize(fitness_function, constraints, mthd, opts)


newp = result.minimizer
newsol = solve(remake(oprob2, p=newp))
newsol[1, end] # Get the steady-state values of the ODEProblem solution

trace_vals = [x.value for x in result.trace]
trace_time = [x.iteration for x in result.trace]

plot(trace_time,trace_vals, xlabel="Iteration", ylabel="Cost", label="Cost", title="Cost vs Iteration", legend=:topleft)