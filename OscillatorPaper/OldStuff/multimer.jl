using Random 
using DataFrames
using Plots 
using Catalyst
using DifferentialEquations

rn1 = @reaction_network begin
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





rn2 = @reaction_network begin
    (ka7,kb7), LpA + T <--> LpAT

end


function make_multimer()
    @reaction_network $name begin
        (ka,kb), A + B <--> C 
    end
end

#parameter list
p = [0.47375415262252124, 70.82403936369272, 300.346311110198, 0.624675949351274, 10, 0.001249351898702548, 
    70, 0.03748055696107644, 500, 10, 0.0624675949351274, 50, 1500.0]

#initial condition list
u0 = [200., 50., 100., 0., 150., 0., 0., 0., 100., 0., 0., 0.]

##REACTION NETWORK
param_symbols = [:ka1,:kb1,:kcat1,:ka2,:kb2,:ka3,:kb3,:ka5,:kb5,:kcat5,:ka6,:kb6,:ka7,:kb7,:y]
pmap = (x[1] => x[2] for x in zip(param_symbols,p))

u_symbols = [:L,:Lp,:K,:LK,:A,:LpA,:LpAK,:LpAKL,:P,:LpP,:LpAP,:LpAPLp]
umap = (y[1] => y[2] for y in zip(u_symbols,u0))

#timespan for integration
tspan = (0., 10.)

oprob = ODEProblem(rn, umap, tspan, pmap)
osol = solve(oprob, Tsit5())
plot(osol,linewidth = 1.5, size = (1100,700))