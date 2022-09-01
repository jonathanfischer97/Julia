using DifferentialEquations
using Plots
using Catalyst

theme(:solarized)

#changed to in-place definition 
function my_ode_solver(du, u, p, t)
    #parameter assignment, unpacked from p
    ka1,kb1,kcat1,ka2,kb2,ka3,kb3,ka4,kb4,kcat4,ka5,kb5,kcat5,ka6,kb6,ka7,kb7,kcat7,y = p

    #initial condition assignment, unpacked from u
    L,Lp,K,LK,A,LpA,LpAK,LpAKL,P,LpP,LpAP,LpAPLp = u

    du[1] = (kb1*LK) - (ka1*L*K) + (kb4*LpAKL) - (y*ka4*LpAK*L) + (kcat5*LpP) + (kcat7*LpAPLp)
    du[2] = (kcat1*LK) - (ka2*Lp*A) + (kb2*LpA) + (kcat4*LpAKL) + (kb5*LpP) - (ka5*Lp*P) - (y*ka7*Lp*LpAP) + (kb7*LpAPLp)
    du[3] = (kb1*LK) - (ka1*L*K) + (kcat1*LK) - (ka3*LpA*K) + (kb3*LpAK)
    du[4] = (ka1*L*K) - (kb1*LK) - (kcat1*LK)
    du[5] = (kb2*LpA) - (ka2*Lp*A)
    du[6] = (ka2*Lp*A) - (kb2*LpA) + (kb3*LpAK) - (ka3*LpA*K) - (ka6*LpA*P) + (kb6*LpAP)
    du[7] = (ka3*LpA*K) - (kb3*LpAK) + (kb4*LpAKL) - (y*ka4*LpAK*L) + (kcat4*LpAKL)
    du[8] = (y*ka4*LpAK*L) - (kb4*LpAKL) - (kcat4*LpAKL)
    du[9] = (kb5*LpP) + (kcat5*LpP) - (Lp*P*ka5) -(ka6*LpA*P) + (kb6*LpAP)
    du[10] = (ka5*Lp*P) - (kb5*LpP) - (kcat5*LpP)
    du[11] = (ka6*LpA*P) - (kb6*LpAP) - (y*ka7*Lp*LpAP) + (kb7*LpAPLp) + (kcat7*LpAPLp)
    du[12] = (y*ka7*Lp*LpAP) - (kb7*LpAPLp) - (kcat7*LpAPLp)
    #[dL,dLp,dK,dP,dLK,dA,dLpA,dLpAK,dLpAP,dLpAPLp,dLpAKL,dLpP]
end

#parameter list
p = [0.47375415262252124, 70.82403936369272, 300.346311110198, 0.624675949351274, 10, 0.001249351898702548, 
    70, 0.47375415262252124, 70.82403936369272, 300.346311110198, 0.03748055696107644, 500, 10, 0.0624675949351274, 
    50, 0.03748055696107644, 500, 10, 1500.0]

#initial condition list
u0 = [200., 50., 100., 0., 150., 0., 0., 0., 100., 0., 0., 0.]

#timespan for integration
tspan = (0., 10.)

#construct ODE problem from constructor, {false} means out of place 

prob = ODEProblem(my_ode_solver,u0,tspan,p)

sol = solve(remake(prob, tspan = (0.,10.)), Tsit5()) #solve with non-stiff Tsit5 alg
plot(sol,size = (1000,700)) #time series plot example


##REACTION NETWORK
param_symbols = [:ka1,:kb1,:kcat1,:ka2,:kb2,:ka3,:kb3,:ka4,:kb4,:kcat4,:ka5,:kb5,:kcat5,:ka6,:kb6,:ka7,:kb7,:kcat7,:y]
pmap = [x[1] => x[2] for x in zip(param_symbols,p)]

u_symbols = [:L,:Lp,:K,:LK,:A,:LpA,:LpAK,:LpAKL,:P,:LpP,:LpAP,:LpAPLp]
umap = [y[1] => y[2] for y in zip(u_symbols,u0)]


rn = @reaction_network begin
    (ka1,kb1), L + K <--> LK
    kcat1, LK --> Lp + K 
    (ka2,kb2), Lp + A <--> LpA 
    (ka3,kb3), LpA + K <--> LpAK  
    (ka4*y,kb4), LpAK + L <--> LpAKL
    kcat4, LpAKL --> Lp + LpAK  
    (ka5,kb5), Lp + P <--> LpP 
    kcat5, LpP --> L + P
    (ka6,kb6), LpA + P <--> LpAP 
    (ka7*y,kb7), Lp + LpAP <--> LpAPLp
    kcat7, LpAPLp --> L + LpAP 
end ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 kcat4 ka5 kb5 kcat5 ka6 kb6 ka7 kb7 kcat7 y 

oprob = ODEProblem(rn, umap, tspan, pmap)
osol = solve(oprob, Tsit5())
plot(osol,linewidth = 1.5, size = (1100,700))