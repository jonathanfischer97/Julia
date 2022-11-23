using Random 
# using DataFrames
using Plots 
using Catalyst
using DifferentialEquations
using Combinatorics

timevar = @variables t

macro string_as_symbolicvar(string)
    # string = s*"(t)"
    string=Symbol(string)
    :(@variables ($string)(t))
    # if setvar
    #     return @eval (@variables ($string)(t))
    # else
    #     return @eval ($string)
    # end
end

function string_as_symbolicvar(string)
    string = Symbol(string)
    @variables ($string)(t)
end

macro varname_as_string(arg)
    string(arg)
end

function varname_as_string(arg)
    string(Symbol(arg))
end



monomers = ["A","B","C","D"]

monovars = [string_as_symbolicvar(str)[1] for str in monomers] 

num_mono_rxns = length(collect(combinations(monomers,2)))

function combine_reactants(dimer)::AbstractString
    dimer[1] * dimer[2]
end

function get_dimers(monomers)
    dimers = []
    for di in collect(permutations(monomers,2))
        push!(dimers, combine_reactants(di))
    end
    return dimers
end

dimers = get_dimers(monomers)

function get_reactions(monomers = monomers)
    dimers = get_dimers(monomers)
    oneplusone = []
    twoplusone = []
    twoplustwo = []
    threeplusone = []
    for complex in collect(permutations(vcat(dimers,monomers),2))
        if length(complex[1])+length(complex[2]) == 2
            push!(oneplusone,combine_reactants(complex))
        elseif length(complex[1])+length(complex[2]) == 3
            push!(twoplusone, combine_reactants(complex))
        elseif length(complex[1])+length(complex[2]) == 4
            push!(twoplustwo, combine_reactants(complex))
        end
    end
    for list in [twoplusone,twoplustwo]
        filter!(x -> allunique(x), list)
    end
    for complex in collect(permutations(vcat(twoplusone,monomers),2))
        push!(threeplusone, combine_reactants(complex))
        filter!(x -> allunique(x), threeplusone)
    end
    return [oneplusone, twoplusone, twoplustwo, threeplusone]
end

lists = get_reactions()

# pars = @parameters kma[1:length(collect(combinations(monomers,2)))] kmb[1:length(collect(combinations(monomers,2)))]

macro string_as_parm(str)
    sym = Symbol(str)
    :(@parameters ($sym))
end

function string_as_parm(str)
    sym = Symbol(str)
    @parameters ($sym)
end

function set_vars(lists)
    varlist1 = []
    varlist2 = []
    varlist3 = []
    varlist4 = []
    metalist = [varlist1,varlist2,varlist3,varlist4]
    for (i,list) in enumerate(lists)
        for mol in list 
            push!(metalist[i],string_as_symbolicvar(mol)[1])
        end
    end
    return metalist
end

var_lists = set_vars(lists)

mono_rx = []
parms = []
for di in var_lists[1]
    # println(di)
    # println(typeof(di))
    reactants = varname_as_string(di)
    # println(typeof(reactants))
    r1 = string_as_symbolicvar(reactants[1])
    # println(r1[1]) 
    r2 = string_as_symbolicvar(reactants[2])
    # println(r2[1])
    ka = string_as_parm("ka"*reactants[1:end-3])
    # println(ka[1])
    kb = string_as_parm("kb"*reactants[1:end-3])
    # println(kb[1])
    # println(typeof(kb[1]))

    push!(mono_rx, Reaction((ka[1],kb[1]),[r1[1],r2[1]],[di]))
    # println("hi")
    push!(parms,ka)
    push!(parms,kb)
end

@named rs = ReactionSystem(mono_rx, t)




symlist = (:a, :b, :c)

@variables $symlist

x = "newvar"

@variables :x

#     if di == AB || di == BA  
#         push!(mono_rx, Reaction((kma[1],kmb[1]), [A,B], di))
        
#     # for clusters of the same size, double the rate
#     if (vᵢ[n] == vⱼ[n])
#         # println("hi")x
#         # newrn = @reaction_network :$(n) begin
#         #     k[n], 2*$(X[vᵢ[n]])  --> $(X[sum_vᵢvⱼ[n]])
#         # end k[n]

#         # basern = compose(basern, newrn)
#         push!(rx, Reaction(k[n], [X[vᵢ[n]]], [X[sum_vᵢvⱼ[n]]], [2], [1]))
#     else
#         # newrn = @reaction_network :$(n) begin
#         #     k[n], $(X[vᵢ[n]]) + $(X[vⱼ[n]]) --> $(X[sum_vᵢvⱼ[n]])
#         # end k[n]

#         # basern = compose(basern, newrn)

#         push!(rx, Reaction(k[n], [X[vᵢ[n]], X[vⱼ[n]]], [X[sum_vᵢvⱼ[n]]],
#                            [1, 1], [1]))
#     end
# end


















trimers = get_trimers(dimers)


twomers = collect(combinations(monomers,2))
# twomers = [(combo[1],combo[2]) for combo in twomers]
trimers = collect(combinations(twomers,2))

function checkcombo(x)
    length(union(x[1],x[2]))==length(vcat(x[1],x[2]))
end

filter!(x -> checkcombo(x), trimers)

num_rxns = (length(twomers) * factorial(4-2)) + length(trimers)

# function string_as_varname(s::AbstractString,v::Any)
#     s=Symbol(s)
#     @eval (($s) = ($v))
# end


# string_as_symbolicvar("AB")


### GET UNIQUE COMBOS ###
function get_combos(monomers)
    combos = Set()
    for mo in monomers
        for mo2 in monomers[2:end]
            push!(combos,Set([mo,mo2]))
        end
    end
    return combos
end

get_combos(monomers)

combos = Set()
i = 0
while i < 2
    for mo in monomers
        i += 1
        println(mo)
        for mo2 in monomers[2:end]
            println([mo,mo2])
            push!(combos,Set([mo,mo2]))
        end
        
    end
end






function generate_reactants(monomers)
    for i in 1:length(monomers)
    twomers = collect(combinations(monomers,2))
    for di in twomers
        combine_reactants(di) |> string_as_symbolicvar
    end
end



function make_rxns(vars, parms)
    rxns = []
    for (var,par) in zip(vars,parms)
        push!(rxns, Reaction(par,var))
    push!(rxns, Reaction(k[n], [X[vᵢ[n]]], [X[sum_vᵢvⱼ[n]]], [2], [1]))

function make_rn()
    

trimer_rn = @reaction_network trirn begin
    (ka_m1, kb_m1), MA + MB <--> MAB 
    (ka_m2, kb_m2), MC + MA <--> MAC 
    (ka_m3, kb_m3), MC + MB <--> MBC 
    (ka_m4, kb_m4), MAB + MC <--> MABC
    (ka_m5, kb_m5), MAC + MB <--> MABC
    (ka_m6, kb_m6), MBC + MA <--> MABC 
end ka_m1 kb_m1 ka_m2 kb_m2 ka_m3 kb_m3 ka_m4 kb_m4 ka_m5 kb_m5 ka_m6 kb_m6 

lipid_rn = @reaction_network liprn begin
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




## Parameter
N = 3                      # maximum cluster size
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          

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
default(size = (1100,700), lw=2, markersize=1, xlabel="Time (sec)")
plot(jsol(t)[16,:], label="X1 (monomers)", color = :blue)
plot!(jsol(t)[17,:], label="X2 (dimers)", color=:orange)
plot!(jsol(t)[18,:], label="X3 (trimers)", color=:purple)
plot!(jsol(t)[25,:], label="X10 (10mers)", color=:red)

# plot(osol,linewidth = 1.5, size = (1100,700))
plot!(osol(t)[16,:], color = :blue, line = :dot, label = "")
plot!(osol(t)[17,:], color=:orange, line = :dot, label = "")
plot!(osol(t)[18,:], color=:purple, line = :dot, label = "")
plot!(osol(t)[25,:], color=:red, line = :dot, label = "")



