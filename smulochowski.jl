using Random 
using DataFrames
using Plots 
using Catalyst
using DifferentialEquations

## Parameter
N = 10                       # maximum cluster size
Vₒ = (4π/3)*(10e-06*100)^3   # volume of a monomers in cm³
Nₒ = 1e-06/Vₒ                # initial conc. = (No. of init. monomers) / bulk volume
uₒ = 10000                   # No. of monomers initially
V = uₒ/Nₒ                    # Bulk volume of system in cm³

integ(x) = Int(floor(x))
n        = integ(N/2)
nr       = N%2 == 0 ? (n*(n + 1) - n) : (n*(n + 1)) # No. of forward reactions

# possible pairs of reactant multimers
pair = []
for i = 2:N
    push!(pair,[1:integ(i/2)  i .- (1:integ(i/2))])
end
pair = vcat(pair...)
vᵢ = @view pair[:,1]   # Reactant 1 indices
vⱼ = @view pair[:,2]   # Reactant 2 indices
volᵢ = Vₒ*vᵢ           # cm⁻³
volⱼ = Vₒ*vⱼ           # cm⁻³
sum_vᵢvⱼ = @. vᵢ + vⱼ  # Product index

# set i to  1 for additive kernel, 2  for constant
i = 1
if i==1
    B = 1.53e03                # s⁻¹
    kv = @. B*(volᵢ + volⱼ)/V  # dividing by volume as its a bi-molecular reaction chain
elseif i==2
    C = 1.84e-04               # cm³ s⁻¹
    kv = fill(C/V, nr)
end

# state variables are X, pars stores rate parameters for each rx
@parameters t
@variables k[1:nr]  (X(t))[1:N]
pars = Pair.(collect(k), kv)

# time-span
if i == 1
    tspan = (0. ,2000.)
elseif i == 2
    tspan = (0. ,350.)
end

 # initial condition of monomers
u₀    = zeros(Int64, N)
u₀[1] = uₒ
u₀map = Pair.(collect(X), u₀)   # map variable to its initial value

# vector to store the Reactions in
# rx = []
basern = @reaction_network rn1 begin
    (ka1,kb1), L + K <--> LK
    kcat1, LK --> Lp + K 
    (ka2,kb2), Lp + A <--> LpA 
    (ka3,kb3), LpA + K <--> LpAK  
    (ka1*y,kb1), LpAK + L <--> LpAKL
    kcat1, LpAKL --> Lp + LpAK  
    (ka5,kb5), Lp + P <--> LpP 
    kcat5, LpP --> L + P
    (ka6,kb6), LpA + P <--> LpAP 
    (ka5*y,kb5), Lp + LpAP <--> LpAPLp
    kcat5, LpAPLp --> L + LpAP 
end ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka5 kb5 kcat5 ka6 kb6 y 

# function make_reaction(; R, P, k, name)
#     @reaction_network $name begin
#         $k, $R --> $P
#     end $k
# end

for n = 1:nr
    kon = k[n]
    # for clusters of the same size, double the rate
    if (vᵢ[n] == vⱼ[n])
        # println("hi")x
        newrn = @reaction_network :$(n) begin
            kon, 2*$(X[vᵢ[n]])  --> $(X[sum_vᵢvⱼ[n]])
        end

        basern = compose(basern, newrn)
        # push!(rx, Reaction(k[n], [X[vᵢ[n]]], [X[sum_vᵢvⱼ[n]]], [2], [1]))
    else
        newrn = @reaction_network :$(n) begin
            kon, $(X[vᵢ[n]]) + $(X[vⱼ[n]]) --> $(X[sum_vᵢvⱼ[n]])
        end

        basern = compose(basern, newrn)

        # push!(rx, Reaction(k[n], [X[vᵢ[n]], X[vⱼ[n]]], [X[sum_vᵢvⱼ[n]]],
        #                    [1, 1], [1]))
    end
end

basern
# @named rs = ReactionSystem(rx, t, collect(X), collect(k))

#parameter list
p = [0.47375415262252124, 70.82403936369272, 300.346311110198, 0.624675949351274, 10, 0.001249351898702548, 
    70, 0.03748055696107644, 500, 10, 0.0624675949351274, 50, 1500.0]

#initial condition list
u0 = [200., 50., 100., 0., 150., 0., 0., 0., 100., 0., 0., 0.]

##REACTION NETWORK
@variables ka1 kb2 kcat1 ka2 kb2 ka3 kb3 ka5 kb5 kcat5 ka6 kb6 ka7 kb7 y
param_symbols = [ka1, kb2, kcat1, ka2, kb2, ka3, kb3, ka5, kb5, kcat5, ka6, kb6, ka7, kb7, y]
# param_symbols = [:ka1,:kb1,:kcat1,:ka2,:kb2,:ka3,:kb3,:ka5,:kb5,:kcat5,:ka6,:kb6,:ka7,:kb7,:y]
pmap = [x[1] => x[2] for x in zip(param_symbols,p)]
push!(pmap,pars...)

@variables L Lp K LK A LpA LpAK LpAKL P LpP LpAP LpAPLp T LpAT
u_symbols = [L, Lp, K, LK, A, LpA, LpAK, LpAKL, P, LpP, LpAP, LpAPLp, T, LpAT]
umap = [y[1] => y[2] for y in zip(u_symbols,u0)]
push!(umap,u₀map...)

#timespan for integration
tspan = (0., 10.)

oprob = ODEProblem(basern, umap, tspan, pmap)
osol = solve(oprob, Tsit5())
plot(osol,linewidth = 1.5, size = (1100,700))
