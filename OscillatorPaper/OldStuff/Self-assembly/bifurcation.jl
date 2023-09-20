using Plots 
using Catalyst
using DifferentialEquations
# using Combinatorics
using Colors
using ColorSchemes
# using ColorTypes
using LinearAlgebra
using SpecialFunctions
using ProgressMeter
using Evolutionary
# using Latexify

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
    (ka8,kb8), LpA + T <--> X

    (ka1m, kb1m), X + B <--> XB 
    (ka1m, kb1m), X + C <--> XC 
    (ka1m, kb1m), X + D <--> XD 
    (ka4m, kb4m), XB + C <--> XBC
    (ka4m*y, kb4m), XB + XC <--> XXBC
    (ka4m, kb4m), XC + B <--> XBC
    (ka5m, kb5m), XB + D <--> XBD
    (ka5m*y, kb5m), XB + XD <--> XXBD
    (ka5m, kb5m), XD + B <--> XBD
    (ka6m, kb6m), XC + D <--> XCD
    (ka6m*y, kb6m), XC + XD <--> XXCD
    (ka6m, kb6m), XD + C <--> XCD 
    (ka7m, kb7m), XBC + D <--> XBCD 
    (ka7m, kb7m), XXBC + D <--> XXBCD
    (ka7m*y, kb7m), XXBC + XD <--> XXXBCD
    (ka8m, kb8m), XBD + C <--> XBCD 
    (ka8m, kb8m), XXBD + C <--> XXBCD 
    (ka8m*y, kb8m), XXBD + XC <--> XXXBCD
    (ka9m, kb9m), XCD + B <--> XBCD 
    (ka9m, kb9m), XXCD + B <--> XXBCD 
    (ka9m*y, kb9m), XXCD + XB <--> XXXBCD 
    (ka1m*y, kb1m), XBC + X <--> XXBC 
    (ka1m*y, kb1m), XBD + X <--> XXBD 
    (ka1m*y, kb1m), XCD + X <--> XXCD 
    (ka1m*y, kb1m), XBCD + X <--> XXBCD 
    (ka1m*y, kb1m), XXBCD + X <--> XXXBCD 
end ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 ka7 kb7 kcat7 ka8 kb8 ka1m kb1m ka4m kb4m ka5m kb5m ka6m kb6m ka7m kb7m ka8m kb8m ka9m kb9m y 

#parameter list
p = [14.766654049574964, 247.74197239200362, 500.0, 31.051821195245743, 158.0304365645891, 52.425931974326005, 500.0, 29.68965224106596, 
3.1517004441718446, 90.20783657199675, 364.05856621148223, 96.68089672060627, 100.0, 178.4206713875396, 63.913390412700274, 16.670218374002296, 
31.00096833128655, 187.53746678023433, 24.208313355774997, 299.70860654065797, 62.268232479997536, 500.0, 13.217117644835021, 69.69105864551203, 
79.85635210382746, 247.8070093520743, 6.034069612733998, 28.445138804106804, 0.7817627051993913/0.001]

pmap = Dict(:ka1 => 14.766654049574964, :kb1 => 247.74197239200362, :kcat1 => 500.0, :ka2 => 31.051821195245743, :kb2 => 158.0304365645891, :ka3 => 52.425931974326005, :kb3 => 500.0, :ka4 => 29.68965224106596, :kb4 => 3.1517004441718446, :ka7 => 90.20783657199675, :kb7 => 364.05856621148223, :kcat7 => 96.68089672060627, :ka8 => 100.0, :kb8 => 178.4206713875396, :ka1m => 63.913390412700274, :kb1m => 16.670218374002296, :ka4m => 31.00096833128655, :kb4m => 187.53746678023433, :ka5m => 24.208313355774997, :kb5m => 299.70860654065797, :ka6m => 62.268232479997536, :kb6m => 500.0, :ka7m => 13.217117644835021, :kb7m => 69.69105864551203, :ka8m => 79.85635210382746, :kb8m => 247.8070093520743, :ka9m => 6.034069612733998, :kb9m => 28.445138804106804, :y => 0.7817627051993913/0.001)
#initial condition list
u0 = [0. , 0.2, 0. , 3. , 0.9, 0. , 0. , 0. , 0.3, 0. , 0. , 0. , 0.3,
1. , 0.5, 0., 0.5, 0.0 , 0.5 , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
0. , 0. , 0. ]

# pmap = [x => y for (x,y) in zip(parameters(trimer_rn),p)]
# pmap = symmap_to_varmap(trimer_rn, pmap)

umap = [x => y for (x,y) in zip(states(trimer_rn),u0)]
umap = symmap_to_varmap(trimer_rn, umap)

tspan = (0.,20.)

odesys = convert(ODESystem, trimer_rn; combinatoric_ratelaws=true)
oprob = ODEProblem(trimer_rn, umap, tspan, pmap; jac = true, sparse = true)
osol = solve(oprob, Tsit5(), saveat = 0.001)

### Bifurcation analysis ###
bif_par = :ka1           # bifurcation parameter
p_span = (0.1, 20.)    # interval to vary S over
plot_var = :X          # we will plot X vs S


