using ModelingToolkit, Symbolics, Catalyst

# Define the variables
# @variables L Lp K P LpA LK LpP LpAK LpAP LpAKL LpAPLp
diffvars = @variables t L(t) LpA(t) 
algvars = @variables Lp(t) K(t) P(t) A(t) LK(t) LpP(t) LpAK(t) LpAP(t) LpAKL(t) LpAPLp(t)
parms = @parameters ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 ka7 kb7 kcat7 DF 
masstots = @parameters ktot ptot atot liptot 

"""Full oscillator model"""
osc_rn = @reaction_network osc_rn begin
    (ka1,kb1), L + K <--> LK
    kcat1, LK --> Lp + K 
    (ka2,kb2), Lp + A <--> LpA 
    (ka3,kb3), LpA + K <--> LpAK  
    (ka1*DF,kb1), LpAK + L <--> LpAKL
    kcat1, LpAKL --> Lp + LpAK  
    (ka7,kb7), Lp + P <--> LpP 
    kcat7, LpP --> L + P
    (ka4,kb4), LpA + P <--> LpAP 
    (ka7*DF,kb7), Lp + LpAP <--> LpAPLp
    kcat7, LpAPLp --> L + LpAP
end

full_oscillator_odes = convert(ODESystem,osc_rn)
simplified_oscillator_odes = structural_simplify(full_oscillator_odes)





D = Differential(t)
# Define the differential equations
diffeqs = [D(L) ~ kb1*LK + kb1*LpAKL + kcat7*LpAPLp + kcat7*LpP - ka1*K*L - ka1*DF*L*LpAK,
            D(LpA) ~ kb3*LpAK + kb4*LpAP + ka2*A*Lp - kb2*LpA - ka3*K*LpA - ka4*LpA*P]

# Define the ODE system
@named reduced_oscillator_odes = ODESystem(vcat(algeqs,diffeqs), t)
