using Plots 
using Catalyst
using DifferentialEquations
using Symbolics

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

@variables L A K P LK Lp LpA LpAK LpAKL LpP LpAP LpAPLp
# Define conservation laws 
Ltot = L + Lp + LpA + LpAK + LpAKL + LpP + LpAP + LpAPLp
Atot = A + LpA + LpAK + LpAP + LpAKL + LpAPLp
Ktot = K + LK + LpAK + LpAKL 
Ptot = P + LpP + LpAP + LpAPLp

# Solve the conservation laws for LpAKL, LpP, LpAP, LpAPLp
sol = Symbolics.solve_for([Ltot,Atot,Ktot,Ptot], [LpAKL,LpP,LpAP,LpAPLp])
