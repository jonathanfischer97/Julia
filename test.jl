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

y = 750 #608.5222549026549


p = [:ka1 => 0.05485309578515125, :kb1 => 19.774627209108715, :kcat1 => 240.99536193310848, 
    :ka2 => 1.0, :kb2 => 0.9504699043910143, 
    :ka3 => 41.04322510426121, :kb3 => 192.86642772763489,
    :ka4 => 0.19184180144850807, :kb4 => 0.12960624157489123, 
    :ka5 => 0.6179131289475834, :kb5 => 3.3890271820244195, :kcat5 => 4.622923709012232, 
    :ka6 => 0.05485309578515125*y, :kb6 => 19.774627209108715, :kcat6 => 240.99536193310848, 
    :ka7 => 0.6179131289475834*y, :kb7 => 3.3890271820244195, :kcat7 => 4.622923709012232]

tspan = (0., 100.)

u0 = [:L => 0., :Lp => 3.0, :K => 0.2, :P => 0.3, :LK => 0., :A => 0.6, :LpA => 0., :LpAK => 0., :LpAP => 0., :LpAPLp => 0., :LpAKL => 0., :LpP => 0.]

#solve ODEs 
oprob = ODEProblem(rn, u0, tspan, p)
osol = solve(oprob, Tsit5())
#plot(osol)




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
        minInd = max(1, ind - window)
        maxInd = min(arrLen, ind+window)
        sum += std(arrayData[minInd:maxInd])
    end
    sum = sum/numPeaks 
end
 
function getFrequencies(y, jump=1, nS=1000)
    #fft sample rate: 1 sample per 5 minutes
    y = y[1:jump:end]
    res = broadcast(abs,rfft(y))
    #normalize the amplitudes
    res = res/cld(nS, 2) #smallest integer larger than or equal to. Rounding up
end

function costTwo(Y)
    p1 = Y[1,:]
    fftData = getFrequencies(p1)

    indexes = maxima(fftData) 
    if length(indexes) == 0
        return 0
    end
    std = getSTD(indexes, fftData, 1)
    diff = getDif(indexes, fftData)
    std + diff
end

## try the cost function 
costTwo(Array(osol))

getFrequencies(osol[1,:])




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