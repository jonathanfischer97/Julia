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
@parameters k[1:nr]
@variables t (X(t))[1:N]
pars = Pair.(collect(k), kv)

# time-span
if i == 1
    tspan = (0. ,2000.)
elseif i == 2
    tspan = (0. ,350.)
end

 # initial condition of monomers
u₀    = zeros(Int64, N)
u₀[1] = 0 #uₒ
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
    (ka7,kb7), LpA + T <--> LpAT 
end ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka5 kb5 kcat5 ka6 kb6 ka7 kb7 y 

# function make_reaction(; R, P, k, name)
#     @reaction_network $name begin
#         $k, $R --> $P
#     end $k
# end

rx = []
for n = 1:nr
    kon = k[n]
    println(kon)

    # for clusters of the same size, double the rate
    if (vᵢ[n] == vⱼ[n])
        # println("hi")x
        # newrn = @reaction_network :$(n) begin
        #     k[n], 2*$(X[vᵢ[n]])  --> $(X[sum_vᵢvⱼ[n]])
        # end k[n]

        # basern = compose(basern, newrn)
        push!(rx, Reaction(k[n], [X[vᵢ[n]]], [X[sum_vᵢvⱼ[n]]], [2], [1]))
    else
        # newrn = @reaction_network :$(n) begin
        #     k[n], $(X[vᵢ[n]]) + $(X[vⱼ[n]]) --> $(X[sum_vᵢvⱼ[n]])
        # end k[n]

        # basern = compose(basern, newrn)

        push!(rx, Reaction(k[n], [X[vᵢ[n]], X[vⱼ[n]]], [X[sum_vᵢvⱼ[n]]],
                           [1, 1], [1]))
    end
end


basern
@named rs = ReactionSystem(rx, t, collect(X), collect(k))

basern = compose(basern,rs)

#parameter list
p = [0.47375415262252124, 70.82403936369272, 300.346311110198, 0.624675949351274, 10, 0.001249351898702548, 
    70, 0.03748055696107644, 500, 10, 0.0624675949351274, 50, 0.1, 70, 0.4, 100, 1500.0]

#initial condition list
u0 = [200, 50, 100, 0, 150, 0, 0, 0, 100, 0, 0, 0, 500, 500, 2000]

##REACTION NETWORK
oldparams = @parameters ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka5 kb5 kcat5 ka6 kb6 ka7 kb7 ka8 kb8 y
# param_symbols = [ka1, kb2, kcat1, ka2, kb2, ka3, kb3, ka5, kb5, kcat5, ka6, kb6, ka7, kb7, y]
# param_symbols = [:ka1,:kb1,:kcat1,:ka2,:kb2,:ka3,:kb3,:ka5,:kb5,:kcat5,:ka6,:kb6,:ka7,:kb7,:y]
# pmap = [x[1] => x[2] for x in zip(param_symbols,p)]
pmap = [x => y for (x,y) in zip(oldparams,p)]
pars = [rs.k[n] => kvval for (n,kvval) in zip(1:25,kv)]
pars = [x[1] => 1. for x in pars]
push!(pmap,pars...)

oldvars = @variables L(t) Lp(t) K(t) LK(t) A(t) LpA(t) LpAK(t) LpAKL(t) P(t) LpP(t) LpAP(t) LpAPLp(t) T(t) LpAT(t) Y(t)
# u_symbols = [L(t), Lp(t), K(t), LK(t), A(t), LpA(t), LpAK(t), LpAKL(t), P(t), LpP(t), LpAP(t), LpAPLp(t), T(t), LpAT(t)]
umap = [x => y for (x,y) in zip(oldvars,u0)]
u₀map = [rs.X[n] => u for (n,u) in zip(1:10,u₀)]
push!(umap,u₀map...)

#add coupling reaction between LpAT and X[1]
coupling_rn = @reaction_network crn begin
    (ka8,kb8), LpAT + Y <--> $(rs.X[1])
end ka8 kb8

fullrn = extend(coupling_rn, basern)
fullrn = Catalyst.flatten(fullrn)

#timespan for integration
tspan = (0., 1.)

jumpsys = convert(JumpSystem,fullrn)
dprob = DiscreteProblem(jumpsys, umap, tspan, pmap)
jprob = JumpProblem(jumpsys, dprob, Direct(), save_positions=(false,false))
jsol = solve(jprob, SSAStepper(), saveat = tspan[2]/500)

t = jsol.t 

default(size = (1100,700), lw=2, markersize=2, xlabel="Time (sec)")
plot(collect(0:0.002:1.), jsol(t)[16,:], label="X1 (monomers)", markercolor=:blue)
plot!(collect(0:0.002:1.), jsol(t)[17,:], label="X2 (dimers)", markercolor=:orange)
plot!(collect(0:0.002:1.), jsol(t)[18,:], label="X3 (trimers)", markercolor=:purple)
plot!(collect(0:0.002:1.), jsol(t)[25,:], label="X10 (10mers)", markercolor=:red)


oprob = ODEProblem(fullrn, umap, tspan, pmap)
osol = solve(oprob, Tsit5())
plot(osol,linewidth = 1.5, size = (1100,700))
