using Random 
using DataFrames
using Plots 
using Catalyst
using DifferentialEquations

## Parameter
N = 10                       # maximum cluster size
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          

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
sum_vᵢvⱼ = @. vᵢ + vⱼ  # Product index


# state variables are X, pars stores rate parameters for each rx
@parameters k[1:nr]
@variables t (X(t))[1:N]
pars = Pair.(collect(k), 1.)


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
pars = [rs.k[n] => val[2] for (n,val) in zip(1:length(pars),pars)]
push!(pmap,pars...)

oldvars = @variables L(t) Lp(t) K(t) LK(t) A(t) LpA(t) LpAK(t) LpAKL(t) P(t) LpP(t) LpAP(t) LpAPLp(t) T(t) LpAT(t) Y(t)
# u_symbols = [L(t), Lp(t), K(t), LK(t), A(t), LpA(t), LpAK(t), LpAKL(t), P(t), LpP(t), LpAP(t), LpAPLp(t), T(t), LpAT(t)]
umap = [x => y for (x,y) in zip(oldvars,u0)]
u₀map = [rs.X[n] => u for (n,u) in zip(1:N,u₀)]
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
saveinterval = tspan[2]/300
jsol = solve(jprob, SSAStepper(), saveat = saveinterval)

t = jsol.t 

oprob = ODEProblem(fullrn, umap, tspan, pmap)
osol = solve(oprob, Tsit5())

#PLOTTING 
default(size = (1000,700), lw=2, markersize=1, xlabel="Time (sec)", ylabel="Copies")
plot(jsol(t)[16,:], label="Monomers", color = :blue)
plot!(jsol(t)[17,:], label="Dimers", color=:orange)
plot!(jsol(t)[18,:], label="Trimers", color=:purple)
plot!(jsol(t)[25,:], label="10mers", color=:red)

# plot(osol,linewidth = 1.5, size = (1100,700))
plot!(osol(t)[16,:], color = :blue, line = :dot, label = "")
plot!(osol(t)[17,:], color=:orange, line = :dot, label = "")
plot!(osol(t)[18,:], color=:purple, line = :dot, label = "")
plot!(osol(t)[25,:], color=:red, line = :dot, label = "")

title!("10mer Formation on Dynamic Lipid Membrane")

savefig("10mer.png")



