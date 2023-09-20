using Random 
using Plots 
using Catalyst
using DifferentialEquations
using Combinatorics
using Latexify

trimer_rn = @reaction_network trirn begin
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
    (ka7,kb7), LpA + T <--> X

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
end ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka5 kb5 kcat5 ka6 kb6 ka7 kb7 y ka1m kb1m ka4m kb4m ka5m kb5m ka6m kb6m ka7m kb7m ka8m kb8m ka9m kb9m  

# lipid_rn = @reaction_network liprn begin
#     (ka1,kb1), L + K <--> LK
#     kcat1, LK --> Lp + K 
#     (ka2,kb2), Lp + A <--> LpA 
#     (ka3,kb3), LpA + K <--> LpAK  
#     (ka1*y,kb1), LpAK + L <--> LpAKL
#     kcat1, LpAKL --> Lp + LpAK  
#     (ka5,kb5), Lp + P <--> LpP 
#     kcat5, LpP --> L + P
#     (ka6,kb6), LpA + P <--> LpAP 
#     (ka5*y,kb5), Lp + LpAP <--> LpAPLp
#     kcat5, LpAPLp --> L + LpAP 
#     (ka7,kb7), LpA + T <--> X 
# end ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka5 kb5 kcat5 ka6 kb6 ka7 kb7 y

# fullrn = extend(lipid_rn, trimer_rn)

#parameter list
p = [0.47375415262252124, 70.82403936369272, 300.346311110198, 0.624675949351274, 10, 0.001249351898702548, 
    70, 0.03748055696107644, 500, 10, 0.0624675949351274, 50, 0.1, 70, 1500.0]

#initial condition list
u0 = [200, 50, 100, 0, 150, 0, 0, 0, 100, 0, 0, 0, 500, 500, 2000]

pmap = (:ka1 => 0.47375415262252124, :kb1 => 70.82403936369272, :kcat1 => 300.346311110198, :ka2 => 0.624675949351274, :kb2 => 10, :ka3 => 0.001249351898702548,
        :kb3 => 70, :ka5 => 0.03748055696107644, :kb5 => 500, :kcat5 => 10, :ka6 => 0.0624675949351274, :kb6 => 50, :ka7 => 0.1, :kb7 => 70, :y => 1500.0,
        :ka1m => abs(randn()), :kb1m => abs(randn()*100), :ka4m => abs(randn()), :kb4m => abs(randn()*100), :ka5m => abs(randn()), :kb5m => abs(randn()*100), 
        :ka6m => abs(randn()), :kb6m => abs(randn()*100),
        :ka7m => abs(randn()), :kb7m => abs(randn()*100), :ka8m => abs(randn()), :kb8m => abs(randn()*100), :ka9m => abs(randn()), :kb9m => abs(randn()*100))

pmap = symmap_to_varmap(trimer_rn, pmap)

umap = [:L => 200, :Lp => 50, :K => 100, :LK => 0, :A => 100, :LpA => 0, :LpAK => 0, :LpAKL => 0, :P => 100, :LpP => 0, :LpAP => 0, :LpAPLp => 0,
        :T => 500, :X => 0, :B => 100, :XB => 0, :C => 100, :XC => 0, :D => 100, :XD => 0, :XBC => 0, :XXBC => 0, :XBD => 0, :XXBD => 0, :XCD => 0, :XXCD => 0,
        :XBCD => 0, :XXBCD => 0, :XXXBCD => 0]

umap = symmap_to_varmap(trimer_rn, umap)

tspan = (0.,10.)

oprob = ODEProblem(trimer_rn, umap, tspan, pmap)
osol = solve(oprob, Tsit5())
plot(osol.t, osol[4,:], color = :green, label = "PIP2")
plot!(osol.t, osol[16,:], color = :blue, line = :dot, label = "XB")
plot!(osol.t, osol[18,:], color=:orange, line = :dot, label = "XC")
plot!(osol.t, osol[20,:], color=:purple, line = :dot, label = "XD")
plot!(osol.t, osol[end,:], color=:red, line = :dot, label = "XXXBCD")

title!("10mer Formation on Dynamic Lipid Membrane")

osys = convert(ODESystem, trimer_rn)

latexify(osys.eqs)

osys.eqs[29]