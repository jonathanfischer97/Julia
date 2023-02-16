using Catalyst
using DifferentialEquations
using Plots
using Evolutionary
using StaticArrays
using BenchmarkTools
using FFTW
using Zygote
using Statistics




## COST FUNCTIONS and dependencies 
function getDif(indexes, arrayData)
    arrLen = length(indexes)
    sum = 0
    for (i, ind) in enumerate(indexes)
        if i == arrLen
            break 
        end
        sum += arrayData[ind] - arrayData[indexes[i+1]]
    end
    sum += arrayData[indexes[end]] 
    return sum #return statement is unneeded but just covering bases 
end
    
function getSTD(indexes, arrayData, window) #get standard deviation of fft peak indexes
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
 
function getFrequencies1(y) #normally would slice y by interval for a sampling rate but not here
    y1 = y
    #fft sample rate: 1 sample per 5 minutes
    res = broadcast(abs,rfft(y1)) #broadcast abs to all elements of rfft return array. Think .= syntax does the same 
    #normalize the amplitudes
    norm_res = res/cld(1000, 2)
    return norm_res #smallest integer larger than or equal to. Rounding up
end


#homemade peakfinder
function findlocalmaxima(signal::Vector)
    inds = Int[]
    buff = Zygote.Buffer(inds) #Zygote.Buffer creates mutable array that is invisible to autodifferention, aka no matrix operations allowed
    if length(signal)>1
        if signal[1]>signal[2]
            push!(buff,1)
        end
        for i=2:length(signal)-1
            if signal[i-1]<signal[i]>signal[i+1]
                push!(buff,i)
            end
        end
        if signal[end]>signal[end-1]
            push!(buff,length(signal))
        end
    end
    return copy(buff) #return copy of mutable buffer array to make immutable (essential for Zygote pullback)
  end

#current testing cost function, USE 
function testcost(p)
    Y = Array(solve(oprob, tspan = (0., 20.), p = p, Tsit5(), saveat = 0.001))
    p1 = Y[1,:]
    fftData = getFrequencies1(p1)
    indexes = findlocalmaxima(fftData)[1]
    if length(indexes) == 0
        return 0
    end
    std = getSTD(indexes, fftData, 1)
    diff = getDif(indexes, fftData)
    return std + diff
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
end ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 ka7 kb7 kcat7 y 

p = [:ka1 => 0.05485309578515125, :kb1 => 19.774627209108715, :kcat1 => 240.99536193310848, 
    :ka2 => 1.0, :kb2 => 0.9504699043910143, 
    :ka3 => 41.04322510426121, :kb3 => 192.86642772763489,
    :ka4 => 0.19184180144850807, :kb4 => 0.12960624157489123, 
    :ka7 => 0.6179131289475834, :kb7 => 3.3890271820244195, :kcat7 => 4.622923709012232, :y => 750.]

u0 = [:L => 0., :Lp => 3.0, :K => 0.2, :P => 0.3, :LK => 0., :A => 0.6, :LpA => 0., :LpAK => 0., :LpAP => 0., :LpAPLp => 0., :LpAKL => 0., :LpP => 0.]

u0 = SA[0., 3., 0.2, 0.3, 0., 0.6, 0., 0., 0., 0., 0., 0.]

tspan = (0.,20.)

oprob = ODEProblem(osc_rn, u0, tspan, p)
@btime osol = solve(oprob, Tsit5(), saveat = 0.001)






osys = convert(ODESystem, osc_rn)
p = symmap_to_varmap(osc_rn, p)
oprob = ODEProblem{false}(osys, u0, tspan, p)

@btime osol = solve(oprob, Tsit5(), saveat = 0.001)


plot(osol)



## GENETIC ALGORITHM
# make initial population of 5000 random parameter sets
pop = [rand(1:500, 14) for i in 1:5000]
p_init = [x[2] for x in p]
result = Evolutionary.optimize(testcost, p_init, GA(populationSize = 5000, crossoverRate = 0.5, crossover = TPX, mutationRate = 0.9, mutation = inversion,  metrics = [Evolutionary.AbsDiff(1e-4)]))


newsol = solve(remake(oprob, p = result.minimizer), Tsit5(), saveat = 0.001)
plot(newsol)