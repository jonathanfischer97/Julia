using DifferentialEquations
using Plots
using Distributed
using Catalyst
using Statistics 
using FFTW
using Peaks
#using ModelingToolkit
#using Graphviz_jll

rn = @reaction_network begin
    ka1, L+K --> LK
    kb1, LK --> L+K
    kcat1, LK --> Lp + K 
    ka2, Lp+A --> LpA
    kb2, LpA --> Lp+A 
    ka3, LpA+K --> LpAK 
    kb3, LpAK --> LpA+K 
    ka4, LpA+P --> LpAP 
    kb4, LpAP --> LpA+P 
    ka5, Lp+P --> LpP 
    kb5, LpP --> Lp+P
    kcat5, LpP --> L+P
    ka6, LpAK+L --> LpAKL
    kb6, LpAKL --> LpAK+L 
    kcat6, LpAKL --> Lp+LpAK 
    ka7, Lp+LpAP --> LpAPLp
    kb7, LpAPLp --> Lp+LpAP 
    kcat7, LpAPLp --> L+LpAP 
end ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 ka5 kb5 kcat5 ka6 kb6 kcat6 ka7 kb7 kcat7 

y = 608.5222549026549


p = [:ka1 => 0.5399439956871153, :kb1 => 32.49731125365913, :kcat1 => 459.0713154241908, :ka2 => 0.6703697205723796, :kb2 => 98.55386776137301, :ka3 => 66.94761309866233, :kb3 => 187.88476784879947,
    :ka4 => 0.732829082837601, :kb4 => 0.5932803922080772, :ka5 => 0.8833353263653272, :kb5 => 79.8349667706061, :kcat5 => 98.51067194455857, :ka6 => 0.5399439956871153*y, :kb6 => 32.49731125365913,
    :kcat6 => 459.0713154241908, :ka7 => 0.8833353263653272*y, :kb7 => 79.8349667706061, :kcat7 => 98.51067194455857]

tspan = (0., 100.)

u0 = [:L => 8.973816043747753, :Lp => 0., :K => 0.6170991018130821, :P => 0.7934153177539267, :LK => 0., :A => 2.6856696379362606, :LpA => 0., :LpAK => 0., :LpAP => 0., :LpAPLp => 0., :LpAKL => 0., :LpP => 0.]

#solve ODEs 
oprob = ODEProblem(rn, u0, tspan, p)
osol = solve(oprob, Tsit5())

species(rn)

## for bifurcation analysis
bif_par = :ka1 
p_span = (0.1, 1.)
plot_var = :L

p_bstart = copy(p)
p_bstart[bif_par] = p_span[1]

bifprob = ODEProblem(rn, u0, (0.0,0.0), p_bstart; jac = true)
F = (u,p) -> bifprob.f(u, p, 0)
J = (u,p) -> bifprob.f.jac(u, p, 0)

## conserved equations 
conservationlaws(rn)
conservedequations(rn)

Z = complexstoichmat(rn)
B = incidencemat(rn)
incidencematgraph(rn)

#check if weakly reversible
subnets = subnetworks(rn)
isweaklyreversible(rn, subnets)

#complex graph 
complexgraph(rn)

#homotopy continuation
ns = convert(NonlinearSystem,rn)
subs = Dict(Pair.(ModelingToolkit.parameters(ns),last.(p)))
new_eqs = map(eq -> substitute(eq.rhs,subs), equations(ns))

using HomotopyContinuation
sols = real_solutions(as_polynomial((f, x...) -> HomotopyContinuation.solve(collect(x)), new_eqs...))

#deficiency
lcs = linkageclasses(rn)
deficiency(rn)


## COST FUNCTION 
function getDif(indexes, arrayData)
    arrLen = length(indexes)
    sum = 0
    for (i, ind) in enumerate(indexes)
        if i == arrLen - 1
            break 
        end
        sum += arrayData[ind] - arrayData[indexes[i+1]]
    end
    sum += arrayData[indexes[end]] #not sure if i need "indexes[end]" to be a list
end
    
function getSTD(indexes, arrayData, window)
    numPeaks = length(indexes)
    arrLen = length(arrayData)
    sum = 0
    for ind in indexes 
        minInd = max(0, ind - window)
        maxInd = min(arrLen, ind+window)
        sum += std(arrayData[minInd:maxInd])
    end
    sum = sum/numPeaks 
end
 
function getFrequencies(y, jump=10, nS=10000)
    #fft sample rate: 1 sample per 5 minutes
    y = y[1:jump:end]
    res = broadcast(abs,rfft(y))
    #normalize the amplitudes
    res = res/cld(nS, 2) #smallest integer larger than or equal to. Rounding up
end

function costTwo(Y)
    p1 = Y[:,5]
    fftData = getFrequencies(p1)

    indexes = maxima(fftData) 
    if length(indexes) == 0
        return 0
    end
    std = getSTD(indexes, fftData, 1)
    diff = getDif(indexes, fftData)
    cost = std + diff
end

## try the cost function 
costTwo(Array(osol))

getFrequencies(Array(osol))

Array(osol)