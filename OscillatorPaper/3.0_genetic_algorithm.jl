using Catalyst
using DifferentialEquations
using Plots
using Evolutionary
using FFTW
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
    # (ka2,kb2), Lp + AK <--> LpAK
    # (ka2*y,kb2), Lp + AKL <--> LpAKL
    # (ka2,kb2), Lp + AP <--> LpAP
    # (ka2*y,kb2), Lp + APLp <--> LpAPLp
    # (ka3,kb3), A + K <--> AK
    # (ka4,kb4), A + P <--> AP
    # (ka3,kb3), A + LK <--> AKL
    # (ka4,kb4), A + LpP <--> APLp
    # (ka3*y,kb3), LpA + LK <--> LpAKL
    # (ka4*y,kb4), LpA + LpP <--> LpAPLp
    # (ka1,kb1), AK + L <--> AKL #binding of kinase to lipid
    # kcat1, AKL --> Lp + AK #phosphorylation of lipid
    # (ka7,kb7), AP + Lp <--> APLp #binding of phosphatase to lipid
    # kcat7, APLp --> L + AP #dephosphorylation of lipid
end 


p = [0.05485309578515125, 19.774627209108715, 240.99536193310848, 
    1.0, 0.9504699043910143, 
    41.04322510426121, 192.86642772763489,
    0.19184180144850807, 0.12960624157489123, 
    0.6179131289475834, 3.3890271820244195, 4.622923709012232,750]
    
    
u0 = [0.0, 0.3, 0.0, 3.0, 0.6, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0]

tspan = (0.,100.)
oprob = ODEProblem(osc_rn, u0, tspan, p)
osol = solve(oprob, save_idxs= 12)
plot(osol)


## REAL FITNESS FUNCTION ##
function eval_fitness(p::Vector{Float64}, prob::ODEProblem, idxs)::Float64
    Y = solve(remake(prob, p=p), saveat=0.1, save_idxs = idxs)
    # p1 = Y[5,:] #get first species
    fftData = getFrequencies(Y.u, length(Y.t)) #get Fourier transform of p1
    indexes = findmaxima(fftData)[1] #get indexes of local maxima of fft
    # realpeaks = length(findmaxima(Y.u, 10)[1]) #get indexes of local maxima of p1
    if length(indexes) == 0 #|| realpeaks < 4  #if no peaks, return 0
        return 0.
    end
    std = getSTD(indexes, fftData, 1) #get standard deviation of fft peak indexes
    # println(std)
    diff = getDif(indexes, fftData) #get difference between peaks
    # println(diff)
    return - std - diff 
end


function make_fitness_function(prob::ODEProblem; idxs = [5])::Function # Create a fitness function that includes your ODE problem as a constant
    function fitness_function(p::Vector{Float64})
        return eval_fitness(p, prob, idxs)  
    end
    return fitness_function
end


fitness_function = make_fitness_function(oprob) # Create a fitness function that includes your ODE problem as a constant



ka_min, ka_max = 0.0, 100.
kb_min, kb_max = 0.00, 500.0
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



#Real optimization
constraints = BoxConstraints([param_values[p]["min"] for p in keys(param_values)], [param_values[p]["max"] for p in keys(param_values)])
opts = Evolutionary.Options(show_trace=true,show_every=1, store_trace=true, iterations=30, parallelization=:thread)

common_range = 0.5
valrange = fill(common_range, 13)
mthd = GA(populationSize = 10000, selection = tournament(1000;select=argmin),
          crossover = TPX, crossoverRate = 0.5,
          mutation  = BGA(valrange, 2), mutationRate = 0.9)
result = Evolutionary.optimize(fitness_function, constraints, mthd, opts)

newp = result.minimizer
newsol = solve(remake(oprob, p=newp))
newsol1 = newsol[1,:]
newsol1dense = solve(remake(oprob, p=newp), saveat=0.1, save_idxs = 1)
getSTD(newsol1, 1)
getSTD(newsol1dense.u, 100)
# plot(newsol, vars=(1,2), title="ODE Solution", xlabel="Time", ylabel="Concentration", label="Concentration")
plot(newsol)
#save newp vector to text file, in text not binary
open("newp.txt", "w") do io
    writedlm(io, newp)
end





results.trace.entries[end].population


eval_fitness(result.minimizer,oprob)

trace_vals = [x.value for x in result.trace]
trace_time = [x.iteration for x in result.trace]

plot(trace_time,trace_vals, xlabel="Iteration", ylabel="Cost", label="Cost", title="Cost vs Iteration", legend=:topleft)