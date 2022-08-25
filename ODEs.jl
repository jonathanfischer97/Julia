using DifferentialEquations
using Plots
#using Distributed
using Catalyst
using DiffEqFlux, Flux
using ForwardDiff
using ReverseDiff
using Zygote
using Statistics 
using FFTW
using Peaks
#using Findpeaks
using PyCall
#np = pyimport("numpy")
#peak = pyimport("peakutils")

#define system of ODEs, out of place (returns new array allocation)
#would be my_ode_solver!(du, u, p, t) for in place, where du is interpreted as a preallocated array that is mutated
function my_ode_solver(u::Array, p::Array, t)
    #parameter assignment, unpacked from p
    ka1,kb1,kcat1,ka2,kb2,ka3,kb3,ka4,kb4,ka5,kb5,kcat5,y = p
    ka6,kb6,kcat6 = ka1,kb1,kcat1
    ka7,kb7,kcat7 = ka5,kb5,kcat5

    #initial condition assignment, unpacked from u
    L,Lp,K,P,LK,A,LpA,LpAK,LpAP,LpAPLp,LpAKL,LpP = u

    dL = (kb1*LK) - (ka1*L*K) + (kcat5*LpAPLp) + (kb6*LpAKL) - (y*ka6*LpAK*L) + (kcat7*LpP)
    dLp = (kcat1*LK) + (kb2*LpA) - (ka2*Lp*A) + (kb5*LpAPLp) - (y*ka5*Lp*LpAP) + (kcat6*LpAKL) - (ka7*Lp*P) + (kb7*LpP)
    dK = (kb1*LK) - (ka1*L*K) + (kcat1*LK) + (kb3*LpAK) - (ka3*LpA*K)
    dP = (kb4*LpAP) - (ka4*LpA*P) - (ka7*Lp*P) + (kb7*LpP) + (kcat7*LpP)
    dLK = (ka1*L*K) - (kb1*LK) - (kcat1*LK)
    dA = (kb2*LpA) - (ka2*Lp*A)
    dLpA = (ka2*Lp*A) - (kb2*LpA) + (kb3*LpAK) - (ka3*LpA*K) + (kb4*LpAP) - (ka4*LpA*P)
    dLpAK = (ka3*LpA*K) - (kb3*LpAK) + (kb6*LpAKL) - (y*ka6*LpAK*L) + (kcat6*LpAKL)
    dLpAP = (ka4*LpA*P) - (kb4*LpAP) + (kb5*LpAPLp) - (y*ka5*LpAP*Lp) + (kcat5*LpAPLp)
    dLpAPLp = (y*ka5*LpAP*Lp) - (kb5*LpAPLp) - (kcat5*LpAPLp)
    dLpAKL = (y*ka6*LpAK*L) - (kb6*LpAKL) - (kcat6*LpAKL)
    dLpP = (ka7*Lp*P) - (kb7*LpP) - (kcat7*LpP)
    [dL,dLp,dK,dP,dLK,dA,dLpA,dLpAK,dLpAP,dLpAPLp,dLpAKL,dLpP]
end

#parameter list
p = [0.05485309578515125, 19.774627209108715, 240.99536193310848, 
    1.0, 0.9504699043910143, 
    41.04322510426121, 192.86642772763489,
    0.19184180144850807, 0.12960624157489123, 
    0.6179131289475834, 3.3890271820244195, 4.622923709012232,750]

#initial condition list
u0 = [0., 3.0, 0.2, 0.3, 0., 0.6, 0., 0., 0., 0., 0., 0.]

#timespan for integration
tspan = (0., 100.)

#construct ODE problem from constructor, {false} means out of place 
prob = ODEProblem{false}(my_ode_solver,u0,tspan,p)

sol = solve(prob, Tsit5()) #solve with non-stiff Tsit5 alg
plot(sol) #time series plot example
plot(sol, vars=(3,4), leg=false) #phase plot example

#calculate gradient with respect to initial conditions, parameters through autodiff of the integrator
du01,dp1 = Zygote.gradient((u0,p)->sum(solve(prob,Tsit5(),u0=u0,p=p,saveat=0.1,sensealg=ForwardDiffSensitivity())),u0,p)





## COST FUNCTIONS and dependencies 
function getDif(indexes, arrayData)
    arrLen = length(indexes)
    sum = 0
    for (i, ind) in enumerate(indexes)
        if i == arrLen - 1
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

##This cost function is the goal but doesn't yet work, IGNORE
function costTwo(y)
    Y = Array(solve(prob, tspan = (0., 100.), p = y, Tsit5()))
    p1 = Y[1,:]
    fftData = getFrequencies1(p1)

    indexes, _ = findmaxima(fftData) 
    if length(indexes) == 0
        return 0
    end
    std = getSTD(indexes, fftData, 1)
    diff = getDif(indexes, fftData)
    return std + diff
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
    Y = Array(solve(prob, tspan = (0., 100.), p = p, Tsit5()))
    p1 = Y[1,:]
    #fftData = getFrequencies1(p1)
    indexes = findlocalmaxima(p1)[1]
    if length(indexes) == 0
        return 0
    end
    std = getSTD(indexes, p1, 1)
    #diff = getDif(indexes, p1) #PROBLEM
    return std #+ diff
    #indexes[1]
    #Tuple(x)[1]
    
end

dcostTwo = Zygote.gradient(testcost, p)






#IGNORE 
function optimize(init,tend)
    function costTwo(y)
        Y = solve(remake(oprob, tspan = (0.,tend), p = y), Tsit5())
        p1 = Y[1,:]
        fftData = getFrequencies(p1)

        indexes, _ = findmaxima(fftData) 
        if length(indexes) == 0
            return 0
        end
        std = getSTD(indexes, fftData, 1)
        diff = getDif(indexes, fftData)
        std + diff
    end
    return DiffEqFlux.sciml_train(costTwo,init,ADAM(0.1),maxiters = 100)
end

#didn't work
estimate = optimize([0.,3.,0.2,0.3,0.,0.6,0.,0.,0.,0.,0.,0.],10.).minimizer

#worked but returned all NaN
du01,dp1 = gradient((u0,p)->sum(solve(oprob,Tsit5(),u0=u0,p=p,saveat=0.1,sensealg=ReverseDiffAdjoint())),u0,p)

#didn't work
ForwardSensitivity(oprob, u0, tspan, p)

#didn't work
forwardprob = ODEForwardSensitivityProblem(oprob, u0, tspan, p)










