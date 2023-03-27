using Plots 
using Catalyst
using DifferentialEquations
using Peaks
using Statistics
# using MKL 
# using Symbolics
# using Latexify

#The oscillator model has 12 variables, or 12 unknowns, and 12 differential
#equations plus 4 mass conservation equations. 

#THIS ONE WORKS WITH 4 VARIABLE UNKNOWNS:
#L, Lp, LpA, and LpAP.

#So we eliminated A, LpAKL, LpAPLp and K using mass conservation--this is
#exact.
#We eliminated LK and LpAK and P and LpP by 4 steady-state approxes. 
#set dLK=0 and dLpAKL=0, solve for LK and LpAK. 
#set dLpP=0 and dK=0, solve for P and LpP

#L+K is reaction 1: ka1, kb1, kcat1
#Lp+A is reaction 2: ka2, kb2
#A+K is reaction 3: ka3, kb3
#A+P is reaction 4: ka4, kb4
#Lp+P is reaction 7: ka7, kb7, kcat7


# """Given the unknowns, calculate the dependent variables"""
# function calc_other_vars(y, p, tots)
#     L, Lp, LpA, LpAP = y
#     ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, DF = p
#     ktot, ptot, atot, liptot = tots

#     P = ((kb7 + kcat7) * (kb1 * kb3 * (L - liptot + Lp + LpA - LpAP + 2 * ptot) +
#         kb3 * kcat1 * (L - liptot + Lp + LpA - LpAP + 2 * ptot) +
#         ka1 * kb3 * L * (ktot + L - liptot + Lp + LpA - LpAP + 2 * ptot) +
#         ka3 * kb1 * LpA * (ktot + L - liptot + Lp + LpA - LpAP + 2 * ptot) +
#         ka3 * kcat1 * LpA * (ktot + L - liptot + Lp + LpA - LpAP + 2 * ptot) +
#         DF * ka1 * ka3 * L * LpA * (2 * ktot + L - liptot + Lp + LpA - LpAP + 2 * ptot))) /
#         ((2 * kb7 + 2 * kcat7 + ka7 * Lp) * (kb3 * (kcat1 + ka1 * L) +
#         ka3 * (kcat1 + DF * ka1 * L) * LpA + kb1 * (kb3 + ka3 * LpA)))

#     LpP = (ka7 *Lp* (kb1* kb3 *(L - liptot + Lp + LpA - LpAP + 2 *ptot) +
#         kb3 *kcat1 *(L - liptot + Lp + LpA - LpAP + 2* ptot) + 
#         ka1 *kb3* L *(ktot + L - liptot + Lp + LpA - LpAP + 2* ptot) +
#         ka3 *kb1 *LpA *(ktot + L - liptot + Lp + LpA - LpAP + 2* ptot) + 
#         ka3* kcat1 *LpA* (ktot + L - liptot + Lp + LpA - LpAP + 2 *ptot) + 
#         DF *ka1* ka3* L *LpA* (2* ktot + L - liptot + Lp + LpA - LpAP + 
#         2 *ptot)))/((2 *kb7 + 2 *kcat7 + ka7* Lp)* (kb3* (kcat1 + ka1* L) + 
#         ka3 *(kcat1 + DF* ka1* L)* LpA + kb1* (kb3 + ka3* LpA)))

#     LK = (ka1 *L* (kb1 *(ktot + L - liptot + Lp + LpA - LpAP - LpP - 2 *P + 2 *ptot) + 
#         kcat1 *(ktot + L - liptot + Lp + LpA - LpAP - LpP - 2* P + 2 *ptot) + 
#         DF *ka1* L *(2 *ktot + L - liptot + Lp + LpA - LpAP - LpP - 2 *P + 
#         2* ptot)))/((kb1 + kcat1)^2 + 2 *DF *ka1 *(kb1 + kcat1)* L + 
#         DF *ka1^2* L^2)

#     LpAK = -(((kb1 + kcat1)* (kb1 *(L - liptot + Lp + LpA - LpAP - LpP - 2 *P + 2* ptot) + 
#         kcat1* (L - liptot + Lp + LpA - LpAP - LpP - 2 *P + 2 *ptot) + 
#         ka1* L* (ktot + L - liptot + Lp + LpA - LpAP - LpP - 2* P + 
#         2 *ptot)))/((kb1 + kcat1)^2 + 2* DF *ka1* (kb1 + kcat1) *L + 
#         DF *ka1^2* L^2))


#     LpAPLp= -LpAP - LpP - P + ptot
#     LpAKL = 0.5*(-L + liptot - LK - Lp - LpA - LpAK + LpAP + LpP + 2*P - 2*ptot)
#     A= 0.5*(2*atot + L - liptot + LK + Lp - LpA - LpAK - LpAP + LpP)
#     K = 0.5* (2*ktot + L - liptot - LK + Lp + LpA - LpAK - LpAP - LpP - 2*P + 2*ptot)

#     return P, LpP, LK, LpAK, LpAPLp, LpAKL, A, K
# end

"""Given the unknowns, calculate the dependent variables"""
function calc_other_vars(y, p, tots)
    L, Lp, LpA, LpAP = y
    ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, DF = p
    ktot, ptot, atot, liptot = tots

    common_sum = kb1 * kb3 * (L - liptot + Lp + LpA - LpAP + 2 * ptot) +
                 kb3 * kcat1 * (L - liptot + Lp + LpA - LpAP + 2 * ptot) +
                 ka1 * kb3 * L * (ktot + L - liptot + Lp + LpA - LpAP + 2 * ptot) +
                 ka3 * kb1 * LpA * (ktot + L - liptot + Lp + LpA - LpAP + 2 * ptot) +
                 ka3 * kcat1 * LpA * (ktot + L - liptot + Lp + LpA - LpAP + 2 * ptot) +
                 DF * ka1 * ka3 * L * LpA * (2 * ktot + L - liptot + Lp + LpA - LpAP + 2 * ptot)

    denominator = (2 * kb7 + 2 * kcat7 + ka7 * Lp) * (kb3 * (kcat1 + ka1 * L) + ka3 * (kcat1 + DF * ka1 * L) * LpA + kb1 * (kb3 + ka3 * LpA))

    P = common_sum / denominator

    LpP = ka7 * Lp * common_sum / denominator

    LK_denom = (kb1 + kcat1)^2 + 2 * DF * ka1 * (kb1 + kcat1) * L + DF * ka1^2 * L^2
    LK = ka1 * L * (kb1 * (ktot + L - liptot + Lp + LpA - LpAP - LpP - 2 * P + 2 * ptot) + kcat1 * (ktot + L - liptot + Lp + LpA - LpAP - LpP - 2 * P + 2 * ptot) + DF * ka1 * L * 
            (2 * ktot + L - liptot + Lp + LpA - LpAP - LpP - 2 * P + 2 * ptot)) / LK_denom

    LpAK = -((kb1 + kcat1) * (kb1 * (L - liptot + Lp + LpA - LpAP - LpP - 2 * P + 2 * ptot) + kcat1 * (L - liptot + Lp + LpA - LpAP - LpP - 2 * P + 2 * ptot) + ka1 * L * 
            (ktot + L - liptot + Lp + LpA - LpAP - LpP - 2 * P + 2 * ptot))) / LK_denom

    LpAPLp= -LpAP - LpP - P + ptot
    LpAKL = 0.5*(-L + liptot - LK - Lp - LpA - LpAK + LpAP + LpP + 2*P - 2*ptot)
    A= 0.5*(2*atot + L - liptot + LK + Lp - LpA - LpAK - LpAP + LpP)
    K = 0.5* (2*ktot + L - liptot - LK + Lp + LpA - LpAK - LpAP - LpP - 2*P + 2*ptot)

    return P, LpP, LK, LpAK, LpAPLp, LpAKL, A, K
end

"""ODE function for the reduced 4th order oscillator model"""
function reduced_oscillator_odes!(dy, y, p, t)
    L, Lp, LpA, LpAP = y #variables in the model
    ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, DF = p #parameters in the model
    tots = p[end-3:end] #total concentrations of the species as constraints

    P, LpP, LK, LpAK, LpAPLp, LpAKL, A, K = calc_other_vars(y, p, tots) #calculate other variables as "constants"

    dy[1] = kb1*LK + kb1*LpAKL + kcat7*LpAPLp + kcat7*LpP - ka1*K*L - ka1*DF*L*LpAK # L
    dy[2] = kcat1*LK + kb2*LpA + kcat1*LpAKL + kb7*LpAPLp + kb7*LpP - ka2*A*Lp - ka7*Lp*P - ka7*DF*Lp*LpAP # Lp
    dy[3] = kb3*LpAK + kb4*LpAP + ka2*A*Lp - kb2*LpA - ka3*K*LpA - ka4*LpA*P # LpA
    dy[4] = kb7*LpAPLp + kcat7*LpAPLp + ka4*LpA*P - kb4*LpAP - ka7*DF*Lp*LpAP # LpAP
end

"""Full oscillator model for comparison"""
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
end

    
#parameter list
"""ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, DF"""
p = [0.05485309578515125, 19.774627209108715, 240.99536193310848, 
    1.0, 0.9504699043910143, 
    41.04322510426121, 192.86642772763489,
    0.19184180144850807, 0.12960624157489123, 
    0.6179131289475834, 3.3890271820244195, 4.622923709012232,750]

#initial condition list
# L, Lp, K, P, A, LpA, LK, LpAK, LpAKL, LpP, LpAP, LpAPLp = u0

# u0 = [0., 3.0, 0.2, 0.3, 0.2, 0.0, 0., 0., 0., 0., 0., 0.] #working 
u0 = [0., 3.0, 0.0, 0.0, 0.0, 0.0, 0., 0., 0., 0.6, 0.3, 0.2]
L, Lp, LpA, LpAP, LK, LpAK, LpAKL, LpP, LpAPLp, A, K, P = u0

umap = [:L => L, :Lp => Lp, :LpA => LpA, :LpAP => LpAP, :LK => LK, :LpAK => LpAK, :LpAKL => LpAKL, :LpP => LpP, :LpAPLp => LpAPLp, :A => A, :K => K, :P => P]

#conserved quantities
ktot=K+LK+LpAK+LpAKL;
ptot=P+LpP+LpAP+LpAPLp;
atot=A+LpA+LpAK+LpAKL+LpAP+LpAPLp;
ltot=L+LK+LpAKL;
lptot=Lp+LpA+LpAK+LpAKL+LpP+LpAP+2*LpAPLp;
liptot=ltot+lptot;
tots = [ktot, ptot, atot, liptot]

#timespan for integration
tspan = (0., 100.);

#solve the reduced ODEs
prob = ODEProblem(reduced_oscillator_odes!, u0[1:4], tspan, vcat(p, tots))
sol = solve(prob)
@time sol = solve(remake(prob, p=vcat(p, tots)))


#solve the full ODEs

prob2 = ODEProblem(osc_rn, umap, tspan, p)
sol2 = solve(prob2, save_idxs=[1,4,6,11])

#plot the results
p1 = plot(sol, label = ["L" "Lp" "LpA" "LpAP"] ,  lw=2, title="Reduced Oscillator Model", xlabel="Time", ylabel="Concentration");
p2 = plot(sol2, label = ["L" "Lp" "LpA" "LpAP"], lw=2, title="Full Oscillator Model", xlabel="Time", ylabel="Concentration");
plot(p1, p2, layout=(2,1), size=(800,800))




## Genetic algorithm to find the best parameters for the reduced model ## 
using Evolutionary
using FFTW

## COST FUNCTIONS and dependencies 
function getDif(indexes::Vector{Int}, arrayData::Vector{Float64})::Float64 #get summed difference between fft peak indexes
    arrLen = length(indexes)
    if arrLen < 2
        return 0.0
    end
    sum_diff = @inbounds sum(arrayData[indexes[i]] - arrayData[indexes[i+1]] for i in 1:(arrLen-1))
    sum_diff += arrayData[indexes[end]]
    return sum_diff
end

function getSTD(indexes::Vector{Int}, arrayData::Vector{Float64}, window::Int)::Float64
    normalized_frequencies = [(ind - 1) / (2 * length(arrayData)) for ind in indexes]
    sum_std = sum(std(normalized_frequencies[max(1, ind - window):min(end, ind + window)]) for ind in 1:length(normalized_frequencies))
    return sum_std / length(normalized_frequencies)
end

function getFrequencies(y::Vector{Float64}, tlen::Int)::Vector{Float64}
    res = abs.(rfft(y))
    return res ./ cld(tlen, 2)
end

function eval_fitness(p::Vector{Float64}, tots::Vector{Float64},  prob::ODEProblem, idxs)::Float64
    Y = solve(remake(prob, p=vcat(p,tots)), save_idxs=idxs)
    println(length(Y.t))
    fftData = getFrequencies(Y.u, length(Y.t))
    indexes = findmaxima(fftData)[1]
    if isempty(indexes)
        return 0.
    end
    std = getSTD(indexes, fftData, 1)
    diff = getDif(indexes, fftData)
    return -std - diff
end






function getSTD(indexes::Vector{Int}, arrayData::Vector{Float64}, window::Int)::Float64 #get average standard deviation of fft peak indexes
    sum_std = @inbounds sum(std(arrayData[max(1, ind - window):min(length(arrayData), ind + window)]) for ind in indexes)
    return sum_std / length(indexes)
end

function getFrequencies(y::Vector{Float64}, tlen::Int)::Vector{Float64} #get fft of ODE solution
    res = broadcast(abs,rfft(y)) #broadcast abs to all elements of rfft return array. 
    res./cld(tlen, 2) #normalize the amplitudes
end

## REAL FITNESS FUNCTION ##
function eval_fitness(p::Vector{Float64}, prob::ODEProblem, idxs)::Float64
    Y = solve(remake(prob, p=p), saveat=0.1, save_idxs = idxs)
    fftData = getFrequencies(Y.u, length(Y.t)) #get Fourier transform of p1
    indexes = findmaxima(fftData)[1] #get indexes of local maxima of fft
    # realpeaks = length(findmaxima(Y.u, 10)[1]) #get indexes of local maxima of p1
    if length(indexes) == 0 #|| realpeaks < 4  #if no peaks, return 0
        return 0.
    end
    std = getSTD(indexes, fftData, 1) #get standard deviation of fft peak indexes
    diff = getDif(indexes, fftData) #get difference between peaks
    return - std - diff 
end







function make_fitness_function(prob::ODEProblem, tots; idxs = 1)::Function # Create a fitness function that includes your ODE problem as a constant
    function fitness_function(p::Vector{Float64})
        return eval_fitness(p, tots, prob, idxs)  
    end
    return fitness_function
end

fitness_function = make_fitness_function(prob, tots) # Create a fitness function that includes your ODE problem as a constant



## parameter constraint ranges ##
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



#Optimization parameters
constraints = BoxConstraints([param_values[p]["min"] for p in keys(param_values)], [param_values[p]["max"] for p in keys(param_values)])
opts = Evolutionary.Options(show_trace=true,show_every=1, store_trace=true, iterations=30, parallelization=:thread, abstol=1e-6, reltol=1e-6)

common_range = 0.5
valrange = fill(common_range, 13)
mthd = GA(populationSize = 5000, selection = tournament(500),
          crossover = TPX, crossoverRate = 0.5,
          mutation  = BGA(valrange, 2), mutationRate = 0.9)


# function callback(trace)

#Optimization
result = Evolutionary.optimize(fitness_function, constraints, mthd, opts)

newp = result.minimizer
newsol = solve(remake(oprob, p=newp))
plot(newsol)