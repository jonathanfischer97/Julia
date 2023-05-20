using DifferentialEquations
using Plots
using Catalyst

theme(:solarized)

#changed to in-place definition 
function my_ode_solver(du, u, p, t)
    #parameter assignment, unpacked from p
    ka1,kb1,kcat1,ka2,kb2,ka3,kb3,ka4,kb4,kcat4,ka5,kb5,kcat5,ka6,kb6,ka7,kb7,kcat7,y = p

    #initial condition assignment, unpacked from u
    L,Lp,K,LK,A,LpA,LpAK,LpAKL,P,LpP,LpAP,LpAPLp = u

    du[1] = (kb1*LK) - (ka1*L*K) + (kb4*LpAKL) - (y*ka4*LpAK*L) + (kcat5*LpP) + (kcat7*LpAPLp)
    du[2] = (kcat1*LK) - (ka2*Lp*A) + (kb2*LpA) + (kcat4*LpAKL) + (kb5*LpP) - (ka5*Lp*P) - (y*ka7*Lp*LpAP) + (kb7*LpAPLp)
    du[3] = (kb1*LK) - (ka1*L*K) + (kcat1*LK) - (ka3*LpA*K) + (kb3*LpAK)
    du[4] = (ka1*L*K) - (kb1*LK) - (kcat1*LK)
    du[5] = (kb2*LpA) - (ka2*Lp*A)
    du[6] = (ka2*Lp*A) - (kb2*LpA) + (kb3*LpAK) - (ka3*LpA*K) - (ka6*LpA*P) + (kb6*LpAP)
    du[7] = (ka3*LpA*K) - (kb3*LpAK) + (kb4*LpAKL) - (y*ka4*LpAK*L) + (kcat4*LpAKL)
    du[8] = (y*ka4*LpAK*L) - (kb4*LpAKL) - (kcat4*LpAKL)
    du[9] = (kb5*LpP) + (kcat5*LpP) - (Lp*P*ka5) -(ka6*LpA*P) + (kb6*LpAP)
    du[10] = (ka5*Lp*P) - (kb5*LpP) - (kcat5*LpP)
    du[11] = (ka6*LpA*P) - (kb6*LpAP) - (y*ka7*Lp*LpAP) + (kb7*LpAPLp) + (kcat7*LpAPLp)
    du[12] = (y*ka7*Lp*LpAP) - (kb7*LpAPLp) - (kcat7*LpAPLp)
    #[dL,dLp,dK,dP,dLK,dA,dLpA,dLpAK,dLpAP,dLpAPLp,dLpAKL,dLpP]
end

#parameter list
p = [0.47375415262252124, 70.82403936369272, 300.346311110198, 0.624675949351274, 10, 0.001249351898702548, 
    70, 0.47375415262252124, 70.82403936369272, 300.346311110198, 0.03748055696107644, 500, 10, 0.0624675949351274, 
    50, 0.03748055696107644, 500, 10, 1500.0]

#initial condition list
u0 = [200., 50., 100., 0., 150., 0., 0., 0., 100., 0., 0., 0.]

#timespan for integration
tspan = (0., 10.)

#construct ODE problem from constructor, {false} means out of place 

prob = ODEProblem(my_ode_solver,u0,tspan,p)

sol = solve(remake(prob, tspan = (0.,10.)), Tsit5()) #solve with non-stiff Tsit5 alg
plot(sol,size = (1000,700)) #time series plot example


##REACTION NETWORK
param_symbols = [:ka1,:kb1,:kcat1,:ka2,:kb2,:ka3,:kb3,:ka4,:kb4,:kcat4,:ka5,:kb5,:kcat5,:ka6,:kb6,:ka7,:kb7,:kcat7,:y]
pmap = [x[1] => x[2] for x in zip(param_symbols,p)]

u_symbols = [:L,:Lp,:K,:LK,:A,:LpA,:LpAK,:LpAKL,:P,:LpP,:LpAP,:LpAPLp]
umap = [y[1] => y[2] for y in zip(u_symbols,u0)]


rn = @reaction_network begin
    (ka1,kb1), L + K <--> LK
    kcat1, LK --> Lp + K 
    (ka2,kb2), Lp + A <--> LpA 
    (ka3,kb3), LpA + K <--> LpAK  
    (ka4*y,kb4), LpAK + L <--> LpAKL
    kcat4, LpAKL --> Lp + LpAK  
    (ka5,kb5), Lp + P <--> LpP 
    kcat5, LpP --> L + P
    (ka6,kb6), LpA + P <--> LpAP 
    (ka7*y,kb7), Lp + LpAP <--> LpAPLp
    kcat7, LpAPLp --> L + LpAP 
end ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 kcat4 ka5 kb5 kcat5 ka6 kb6 ka7 kb7 kcat7 y 

oprob = ODEProblem(rn, umap, tspan, pmap)
osol = solve(oprob, Tsit5())
plot(osol,linewidth = 1.5, size = (1100,700))


# COST FUNCTIONS and dependencies 
function getDif(indexes, arrayData)
    arrLen = length(indexes)
    println(1)
    sum = 0
    for (i, ind) in enumerate(indexes)
        println(i)
        if i == arrLen
            println(3)
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
    fftData = getFrequencies1(p1)
    indexes = findlocalmaxima(fftData)[1]
    if length(indexes) == 0
        return 0
    end
    std = getSTD(indexes, fftData, 1)
    diff = getDif(indexes, fftData)
    return std + diff
end

dcostTwo = Zygote.gradient(testcost, p)


## for bifurcation analysis
bif_par = :ka2 
p_span = (0.1, 1.)
plot_var = :L

p_bstart = Dict(copy(pmap))
p_bstart[bif_par] = p_span[1]

bifprob = ODEProblem(rn, umap, (0.0,0.0), p_bstart; jac = true)
F = (umap,pmap) -> bifprob.f(umap, pmap, 0)
J = (umap,pmap) -> bifprob.f.jac(umap, pmap, 0)

# get S and X as a symbolic variables
@unpack ka1, L = rn

# find their indices in oprob.p and oprob.u0 respectively
bif_idx  = findfirst(isequal(ka1), parameters(rn))
plot_idx = findfirst(isequal(L), species(rn))

using BifurcationKit, Plots, LinearAlgebra, Setfield

bprob = BifurcationProblem(F, bifprob.u0, bifprob.p, (@lens _[bif_idx]);
                           recordFromSolution = (x, p) -> x[plot_idx], J = J)

bopts = ContinuationPar(dsmax = 0.05,          # Max arclength in PACM.
                        dsmin = 1e-4,          # Min arclength in PACM.
                        ds=0.001,              # Initial (positive) arclength in PACM.
                        maxSteps = 100000,     # Max number of steps.
                        pMin = p_span[1],      # Min p-val (if hit, the method stops).
                        pMax = p_span[2],      # Max p-val (if hit, the method stops).
                        detectBifurcation = 3) # Value in {0,1,2,3}

bf = bifurcationdiagram(bprob, PALC(), 2, (args...) -> bopts)