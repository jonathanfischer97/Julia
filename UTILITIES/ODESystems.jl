using DifferentialEquations
function fullmodel_ode!(du, u, p, t)
    L, K, P, A, Lp, LpA, LK, LpP, LpAK, LpAP, LpAKL, LpAPLp, AK, AP, AKL, APLp = u #initial conditions
    ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, y = p #parameters
    du[1] = kb1*AKL + kcat7*APLp + kb1*LK + kb1*LpAKL + kcat7*LpAPLp + kcat7*LpP - ka1*AK*L - ka1*K*L - ka1*y*L*LpAK #* L
    du[2] = kb1*LK + kb3*AK + kcat1*LK + kb3*LpAK - ka3*A*K - ka1*K*L - ka3*K*LpA #* K
    du[3] = kb4*AP + kb4*LpAP + kb7*LpP + kcat7*LpP - ka4*A*P - ka7*Lp*P - ka4*LpA*P #* P
    du[4] = kb2*LpA + kb3*AK + kb3*AKL + kb4*AP + kb4*APLp - ka2*A*Lp - ka3*A*K - ka3*A*LK - ka4*A*LpP - ka4*A*P #* A
    du[5] = kb7*APLp + kcat1*AKL + kcat1*LK + kb2*LpA + kb2*LpAK + kb2*LpAKL + kb2*LpAP + kcat1*LpAKL + kb2*LpAPLp + kb7*LpAPLp + kb7*LpP - ka2*A*Lp - ka2*AK*Lp - ka2*AP*Lp - ka7*AP*Lp - ka7*Lp*P - ka2*y*AKL*Lp - ka2*y*APLp*Lp - ka7*y*Lp*LpAP #* Lp
    du[6] = kb3*LpAK + kb3*LpAKL + kb4*LpAP + kb4*LpAPLp + ka2*A*Lp - kb2*LpA - ka3*K*LpA - ka4*LpA*P - ka3*y*LK*LpA - ka4*y*LpA*LpP #* LpA
    du[7] = kb3*AKL + kb3*LpAKL + ka1*K*L - kb1*LK - kcat1*LK - ka3*A*LK - ka3*y*LK*LpA #* LK
    du[8] = kb4*APLp + kb4*LpAPLp + ka7*Lp*P - kb7*LpP - kcat7*LpP - ka4*A*LpP - ka4*y*LpA*LpP #* LpP
    du[9] = kb1*LpAKL + kcat1*LpAKL + ka2*AK*Lp + ka3*K*LpA - kb2*LpAK - kb3*LpAK - ka1*y*L*LpAK #* LpAK
    du[10] = kb7*LpAPLp + kcat7*LpAPLp + ka2*AP*Lp + ka4*LpA*P - kb2*LpAP - kb4*LpAP - ka7*y*Lp*LpAP #* LpAP
    du[11] = ka1*y*L*LpAK + ka2*y*AKL*Lp + ka3*y*LK*LpA - kb1*LpAKL - kb2*LpAKL - kb3*LpAKL - kcat1*LpAKL #* LpAKL
    du[12] = ka2*y*APLp*Lp + ka7*y*Lp*LpAP + ka4*y*LpA*LpP - kb2*LpAPLp - kb4*LpAPLp - kb7*LpAPLp - kcat7*LpAPLp #* LpAPLp
    du[13] = kb1*AKL + kb2*LpAK + kcat1*AKL + ka3*A*K - kb3*AK - ka1*AK*L - ka2*AK*Lp #* AK
    du[14] = kb2*LpAP + kb7*APLp + kcat7*APLp + ka4*A*P - kb4*AP - ka2*AP*Lp - ka7*AP*Lp #* AP
    du[15] = kb2*LpAKL + ka3*A*LK + ka1*AK*L - kb1*AKL - kb3*AKL - kcat1*AKL - ka2*y*AKL*Lp #* AKL
    du[16] = kb2*LpAPLp + ka7*AP*Lp + ka4*A*LpP - kb4*APLp - kb7*APLp - kcat7*APLp - ka2*y*APLp*Lp #* APLp
    du
end

# function fullmodel(u, p, t)
#         L, K, P, A, Lp, LpA, LK, LpP, LpAK, LpAP, LpAKL, LpAPLp, AK, AP, AKL, APLp = u #initial conditions
#         ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, y = p #parameters
#         @SVector [kb1*AKL + kcat7*APLp + kb1*LK + kb1*LpAKL + kcat7*LpAPLp + kcat7*LpP - ka1*AK*L - ka1*K*L - ka1*y*L*LpAK, #* L
#                  kb1*LK + kb3*AK + kcat1*LK + kb3*LpAK - ka3*A*K - ka1*K*L - ka3*K*LpA, #* K
#                 kb4*AP + kb4*LpAP + kb7*LpP + kcat7*LpP - ka4*A*P - ka7*Lp*P - ka4*LpA*P, #* P
#                 kb2*LpA + kb3*AK + kb3*AKL + kb4*AP + kb4*APLp - ka2*A*Lp - ka3*A*K - ka3*A*LK - ka4*A*LpP - ka4*A*P, #* A
#                 kb7*APLp + kcat1*AKL + kcat1*LK + kb2*LpA + kb2*LpAK + kb2*LpAKL + kb2*LpAP + kcat1*LpAKL + kb2*LpAPLp + kb7*LpAPLp + kb7*LpP - ka2*A*Lp - ka2*AK*Lp - ka2*AP*Lp - ka7*AP*Lp - ka7*Lp*P - ka2*y*AKL*Lp - ka2*y*APLp*Lp - ka7*y*Lp*LpAP, #* Lp
#                 kb3*LpAK + kb3*LpAKL + kb4*LpAP + kb4*LpAPLp + ka2*A*Lp - kb2*LpA - ka3*K*LpA - ka4*LpA*P - ka3*y*LK*LpA - ka4*y*LpA*LpP, #* LpA
#                 kb3*AKL + kb3*LpAKL + ka1*K*L - kb1*LK - kcat1*LK - ka3*A*LK - ka3*y*LK*LpA, #* LK
#                 kb4*APLp + kb4*LpAPLp + ka7*Lp*P - kb7*LpP - kcat7*LpP - ka4*A*LpP - ka4*y*LpA*LpP, #* LpP
#                 kb1*LpAKL + kcat1*LpAKL + ka2*AK*Lp + ka3*K*LpA - kb2*LpAK - kb3*LpAK - ka1*y*L*LpAK, #* LpAK
#                 kb7*LpAPLp + kcat7*LpAPLp + ka2*AP*Lp + ka4*LpA*P - kb2*LpAP - kb4*LpAP - ka7*y*Lp*LpAP, #* LpAP
#                 ka1*y*L*LpAK + ka2*y*AKL*Lp + ka3*y*LK*LpA - kb1*LpAKL - kb2*LpAKL - kb3*LpAKL - kcat1*LpAKL, #* LpAKL
#                 ka2*y*APLp*Lp + ka7*y*Lp*LpAP + ka4*y*LpA*LpP - kb2*LpAPLp - kb4*LpAPLp - kb7*LpAPLp - kcat7*LpAPLp, #* LpAPLp
#                 kb1*AKL + kb2*LpAK + kcat1*AKL + ka3*A*K - kb3*AK - ka1*AK*L - ka2*AK*Lp, #* AK
#                 kb2*LpAP + kb7*APLp + kcat7*APLp + ka4*A*P - kb4*AP - ka2*AP*Lp - ka7*AP*Lp, #* AP
#                 kb2*LpAKL + ka3*A*LK + ka1*AK*L - kb1*AKL - kb3*AKL - kcat1*AKL - ka2*y*AKL*Lp, #* AKL
#                 kb2*LpAPLp + ka7*AP*Lp + ka4*A*LpP - kb4*APLp - kb7*APLp - kcat7*APLp - ka2*y*APLp*Lp] #* APLp
# end

#? Parameter list
psym = [:ka1 => 0.009433439939827041, :kb1 => 2.3550169939427845, :kcat1 => 832.7213093872278, :ka2 => 12.993995997539924, :kb2 => 6.150972501791291,
        :ka3 => 1.3481451097940793, :kb3 => 0.006201726090609513, :ka4 => 0.006277294665474662, :kb4 => 0.9250191811994848, :ka7 => 57.36471615394549, 
        :kb7 => 0.04411989797898752, :kcat7 => 42.288085868394326, :DF => 3631.050539219606]
p = [x[2] for x in psym]
    
#? Initial condition list
usym = [:L => 0.0, :K => 0.5, :P => 0.3, :A => 2.0, :Lp => 3.0, :LpA => 0.0, :LK => 0.0, 
        :LpP => 0.0, :LpAK => 0.0, :LpAP => 0.0, :LpAKL => 0.0, :LpAPLp => 0.0, :AK => 0.0, :AP => 0.0, 
        :AKL => 0.0, :APLp => 0.0]
u0 = [x[2] for x in usym]
# u0static = @SVector [x[2] for x in usym]

#? Timespan for integration
tspan = (0., 100.)

# fullmodel(u,p,t) = fullmodel!(similar(u),u,p,t)

odeprob = ODEProblem(fullmodel_ode!, u0, tspan, p)
# odeprobstatic = ODEProblem(fullmodel, u0static, tspan, p)
# fullprob = ODEProblem(fullrn, u0, tspan, p)

# @benchmark solve($odeprob, saveat=0.1, save_idxs=1)
# @benchmark solve($odeprobstatic, saveat=0.1, save_idxs=1)
# @benchmark solve($fullprob, saveat=0.1, save_idxs=1)