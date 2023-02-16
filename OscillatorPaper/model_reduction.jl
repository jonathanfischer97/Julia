using Plots 
using Catalyst
using DifferentialEquations
using Symbolics
using Latexify

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

for eq in equations(ns)
    println(eq)
    println("\n")
end

const MT = ModelingToolkit
subs = Dict(MT.parameters(ns) .=> MT.varmap_to_vars([], MT.parameters(ns); defaults=MT.defaults(ns)))

cons_eq = conservedequations(trimer_rn)
new_eqs = map(eq -> substitute(eq.rhs-eq.lhs,subs), [equations(ns)...;cons_eq...])

Symbolics.solve_for(ns.eqs[1],:LK)

@variables t L(t) A(t) K(t) P(t) LK(t) Lp(t) LpA(t) LpAK(t) LpAKL(t) LpP(t) LpAP(t) LpAPLp(t)
# Define conservation laws 
Ltot = L + Lp + LpA + LpAK + LpAKL + LpP + LpAP + LpAPLp
Atot = A + LpA + LpAK + LpAP + LpAKL + LpAPLp
Ktot = K + LK + LpAK + LpAKL 
Ptot = P + LpP + LpAP + LpAPLp

# Solve the conservation laws for LpAKL, LpP, LpAP, LpAPLp
sol = Symbolics.solve_for([Ltot,Atot,Ktot,Ptot], [LpAKL,LpP,LpAP,LpAPLp])

# Substitute the conservation laws into the equations
subs = Dict(MT.parameters(ns) .=> MT.varmap_to_vars([], MT.parameters(ns); defaults=MT.defaults(ns)))
reduced_eqs = map(eq -> substitute(eq.rhs-eq.lhs,Dict(LpAKL,LpP,LpAP,LpAPLp .=> cons_eq), equations(osys)))

subs = Dict(LpAKL => cons_eq[1].rhs, LpP => cons_eq[2].rhs, LpAP => cons_eq[3].rhs, LpAPLp => cons_eq[4].rhs)
substitute(osys.eqs[1],Dict(LpAKL => cons_eq[1].rhs))
reduced_eqs = map(eq -> substitute(eq,subs), [equations(osys)...]); deleteat!(reduced_eqs,[8,10,11,12])

for eq in reduced_eqs
    println(eq)
    println("\n")
end

latexify(reduced_eqs) |> render 






@named daesys = ODESystem(vcat(osys.eqs,ODESystem, cons_eq))

simplified_eqs = structural_simplify(daesys)

for eq in daesys.eqs
    println(eq)
    println("\n")
end





@variables t x(t)  # independent and dependent variables
@parameters τ       # parameters
@constants h = 1    # constants
D = Differential(t) # define an operator for the differentiation w.r.t. time

# your first ODE, consisting of a single equation, indicated by ~
@named fol_model = ODESystem(D(x) ~ (h - x) / τ)
@variables t RHS(t)
@named fol_separate = ODESystem([RHS ~ (h - x) / τ,D(x) ~ RHS])

fol_simplified = structural_simplify(fol_separate)
equations(fol_simplified)