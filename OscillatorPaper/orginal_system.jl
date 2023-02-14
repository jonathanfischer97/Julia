using Plots 
using Catalyst
using DifferentialEquations

trimer_rn = @reaction_network trirn begin
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
end ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 ka7 kb7 kcat7 y 

osys = convert(ODESystem, trimer_rn)

for eq in osys.eqs
    println(eq)
end


p = [:ka1 => 0.05485309578515125, :kb1 => 19.774627209108715, :kcat1 => 240.99536193310848, 
    :ka2 => 1.0, :kb2 => 0.9504699043910143, 
    :ka3 => 41.04322510426121, :kb3 => 192.86642772763489,
    :ka4 => 0.19184180144850807, :kb4 => 0.12960624157489123, 
    :ka7 => 0.6179131289475834, :kb7 => 3.3890271820244195, :kcat7 => 4.622923709012232, :y => 750.]

u0 = [:L => 0., :Lp => 3.0, :K => 0.2, :P => 0.3, :LK => 0., :A => 0.6, :LpA => 0., :LpAK => 0., :LpAP => 0., :LpAPLp => 0., :LpAKL => 0., :LpP => 0.]

setdefaults!(trimer_rn, vcat(p,u0))

ns = convert(NonlinearSystem, trimer_rn; remove_conserved=true)

const MT = ModelingToolkit
subs = Dict(MT.parameters(ns) .=> MT.varmap_to_vars([], MT.parameters(ns); defaults=MT.defaults(ns)))

cons_eq = conservedequations(trimer_rn)
new_eqs = map(eq -> substitute(eq.rhs-eq.lhs,subs), [equations(ns)...;cons_eq...])

using HomotopyContinuation
sols = real_solutions(as_polynomial((f, x...) -> HomotopyContinuation.solve(collect(x)), new_eqs...))


## BIFURCATION DIAGRAM ##
p = Dict(:ka1 => 0.05485309578515125, :kb1 => 19.774627209108715, :kcat1 => 240.99536193310848, 
    :ka2 => 1.0, :kb2 => 0.9504699043910143, 
    :ka3 => 41.04322510426121, :kb3 => 192.86642772763489,
    :ka4 => 0.19184180144850807, :kb4 => 0.12960624157489123, 
    :ka7 => 0.6179131289475834, :kb7 => 3.3890271820244195, :kcat7 => 4.622923709012232, :y => 750.)

bif_par = :y 
p_span = (100., 1000.)
plot_var = :Lp

p_bstart = copy(p)
p_bstart[bif_par] = p_span[1]

oprob = ODEProblem(trimer_rn, u0, (0.0,0.0), p_bstart; jac = true)
F = (u,p) -> oprob.f(u, p, 0)
J = (u,p) -> oprob.f.jac(u, p, 0)

# get S and X as a symbolic variables
@unpack y, Lp = trimer_rn

# find their indices in oprob.p and oprob.u0 respectively
bif_idx  = findfirst(isequal(y), Catalyst.parameters(trimer_rn))
plot_idx = findfirst(isequal(Lp), Catalyst.species(trimer_rn))

using BifurcationKit, Setfield