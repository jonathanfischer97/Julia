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
using PyCall

function rate_converter(rate,V)
    rate = rate * 1e6
    nanoV = V*1e9
    rate/nanoV
end

V = 0.0800415

rn = @reaction_network begin
    (ka1,kb1), Lp + A <--> LpA
    (ka2,kb2), LpA + K <--> LpAK 
    (ka3,kb3), LpA + P <--> LpAP  
    (ka4*y,kb4), LpAK + L <--> LpAKL    
    (ka5*y,kb5), Lp + LpAP <--> LpAPLp
end ka1 kb1 ka2 kb2 ka3 kb3 ka4 kb4 ka5 kb5 y 

p = [0.00556479742949202, 8.90829466905371, 0.03554699824428245, 1.9916645419788135, 0.04602260964544391, 3.683718709935799, 0.023863972780133767, 
0.0636702725760359, 0.023863972780133767, 0.31835136288017957, 1500.0]

param_symbols = [:ka1,:kb1,:ka2,:kb2,:ka3,:kb3,:ka4,:kb4,:ka5,:kb5,:y]
pmap = [x[1] => x[2] for x in zip(param_symbols,p)]
# pmap2 = map(x -> rate_converter(x,0.0800415),pmap)

u0 = [200.,300.,150.,100.,100.,0.,0.,0.,0.,0.]

u_symbols = [:L,:Lp,:A,:K,:P,:LpA,:LpAK,:LpAP,:LpAKL,:LpAPLp]
umap = [y[1] => y[2] for y in zip(u_symbols,u0)]


tspan = (0., 10.)
oprob = ODEProblem(rn, umap, tspan, pmap)
osol = solve(oprob, Tsit5())
plot(osol,linewidth = 1.5, size = (1100,700))