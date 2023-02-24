using Catalyst
using DifferentialEquations
using Plots
using Evolutionary
using StaticArrays
using BenchmarkTools
using FFTW
using Zygote
using Statistics
using Random
using Distributions




## COST FUNCTIONS and dependencies 
function getDif(indexes::Vector{Int}, arrayData::Vector{Float64}) 
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
    
function getSTD(indexes::Vector{Int}, arrayData::Vector{Float64}, window::Int) #get standard deviation of fft peak indexes
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

 
function getFrequencies1(y) #normally would slice y by interval for a sampling rate but not here
    y1 = y
    #fft sample rate: 1 sample per 5 minutes
    res = broadcast(abs,rfft(y1)) #broadcast abs to all elements of rfft return array. Think .= syntax does the same 
    #normalize the amplitudes
    norm_res = res/cld(250, 2)
    return norm_res #smallest integer larger than or equal to. Rounding up
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
function testcost(p::Vector{Float64})
    Y = solve(remake(oprob2, p=p), tspan = (0., 20.), saveat = 0.001)
    p1 = Y[1,:]
    fftData = getFrequencies1(p1) #get Fourier transform of p1
    indexes = findlocalmaxima(fftData) #get indexes of local maxima of fft
    if length(indexes) == 0 #if no peaks, return 0
        return 0
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
    (ka2,kb2), AP + Lp <--> APLp #binding of phosphatase to lipid
end ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 ka7 kb7 kcat7 y 

p = [:ka1 => 0.05485309578515125, :kb1 => 19.774627209108715, :kcat1 => 240.99536193310848, 
    :ka2 => 1.0, :kb2 => 0.9504699043910143, 
    :ka3 => 41.04322510426121, :kb3 => 192.86642772763489,
    :ka4 => 0.19184180144850807, :kb4 => 0.12960624157489123, 
    :ka7 => 0.6179131289475834, :kb7 => 3.3890271820244195, :kcat7 => 4.622923709012232, :y => 750.]

input = "{'ka1': 71.1660475072189, 'kb1': 1.5159502159283402, 'kcat1': 13.424687589058786, 'ka2': 0.01, 'kb2': 0.5917927793632454, 'ka3': 2.8460552266571644, 'kb3': 0.01, 'ka4': 59.51538294987245, 'kb4': 0.13776483888525468, 'ka7': 23.590647636435467, 'kb7': 1.3220997829155607, 'kcat7': 36.90013394947226, 'VA': 525.170659247689}"
p = parse_input(input)

u0 = [:L => 0., :Lp => 3.0, :K => 0.2, :P => 0.3, :LK => 0., :A => 2.0, :LpA => 0., :LpAK => 0., :LpAP => 0., :LpAPLp => 0., :LpAKL => 0., :LpP => 0., :AK => 0., :AP => 0., :AKL => 0., :APLp => 0.]


tspan = (0.,20.)

oprob2 = ODEProblem(osc_rn, u0, tspan, p)
osol = solve(oprob2, Tsit5(), saveat = 0.001)
plot(osol)






u0 = SA[0., 0.2, 0., 3.0, 3.0, 0., 0., 0., 0.3, 0., 0., 0.]
osys = convert(ODESystem, osc_rn)
p = symmap_to_varmap(osc_rn, p)
oprob = ODEProblem{false}(osys, u0, tspan, p)

osol = solve(oprob, Tsit5(), saveat = 0.001)


plot(osol)



## GENETIC ALGORITHM
# define Candidate and Model structs
struct Candidate
    params::Vector{Float64}
    fitness::Float64
end

struct Model
    rn::ReactionSystem  # a Catalyst.jl reaction network
    param_values::Dict{String, Dict{String, Float64}}
    cost_function::Function
    population::Vector{Candidate}
    best_candidate::Candidate
end

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
    "y" => Dict("min" => 100, "max" => 1500)
)

generateCandidate(param_values)


result = Evolutionary.optimize(testcost, generateCandidate(param_values), GA(populationSize = 5000, crossoverRate = 0.5, crossover = TPX, mutationRate = 0.9, mutation = PLM,  metrics = [Evolutionary.AbsDiff(1e-4)]))
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
osol = solve(oprob, Tsit5(), saveat = 0.001)
plot(osol)

newsol = solve(remake(oprob, p = newp), Tsit5(), saveat = 0.001)
plot(newsol)