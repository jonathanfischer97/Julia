using DifferentialEquations
using StaticArrays
using BenchmarkTools, Profile

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
    nothing
end

# function fullmodel_ode!(du, u, p, t)
#         du[1] = p[2]*u[15] + p[12]*u[16] + p[2]*u[7] + p[2]*u[11] + p[12]*u[12] + p[12]*u[8] - p[1]*u[13]*u[1] - p[1]*u[2]*u[1] - p[1]*p[13]*u[1]*u[9] #* L
#         du[2] = p[2]*u[7] + p[7]*u[13] + p[3]*u[7] + p[7]*u[9] - p[6]*u[4]*u[2] - p[1]*u[2]*u[1] - p[6]*u[2]*u[6] #* u[2]
#         du[3] = p[9]*u[14] + p[9]*u[10] + p[11]*u[8] + p[12]*u[8] - p[8]*u[4]*u[3] - p[10]*u[5]*u[3] - p[8]*u[6]*u[3] #* u[3]
#         du[4] = p[5]*u[6] + p[7]*u[13] + p[7]*u[15] + p[9]*u[14] + p[9]*u[16] - p[4]*u[4]*u[5] - p[6]*u[4]*u[2] - p[6]*u[4]*u[7] - p[8]*u[4]*u[8] - p[8]*u[4]*u[3] #* u[4]
#         du[5] = p[11]*u[16] + p[3]*u[15] + p[3]*u[7] + p[5]*u[6] + p[5]*u[9] + p[5]*u[11] + p[5]*u[10] + p[3]*u[11] + p[5]*u[12] + p[11]*u[12] + p[11]*u[8] - p[4]*u[4]*u[5] - p[4]*u[13]*u[5] - p[4]*u[14]*u[5] - p[10]*u[14]*u[5] - p[10]*u[5]*u[3] - p[4]*p[13]*u[15]*u[5] - p[4]*p[13]*u[16]*u[5] - p[10]*p[13]*u[5]*u[10] #* u[5]
#         du[6] = p[7]*u[9] + p[7]*u[11] + p[9]*u[10] + p[9]*u[12] + p[4]*u[4]*u[5] - p[5]*u[6] - p[6]*u[2]*u[6] - p[8]*u[6]*u[3] - p[6]*p[13]*u[7]*u[6] - p[8]*p[13]*u[6]*u[8] #* u[6]
#         du[7] = p[7]*u[15] + p[7]*u[11] + p[1]*u[2]*u[1] - p[2]*u[7] - p[3]*u[7] - p[6]*u[4]*u[7] - p[6]*p[13]*u[7]*u[6] #* u[7]
#         du[8] = p[9]*u[16] + p[9]*u[12] + p[10]*u[5]*u[3] - p[11]*u[8] - p[12]*u[8] - p[8]*u[4]*u[8] - p[8]*p[13]*u[6]*u[8] #* u[8]
#         du[9] = p[2]*u[11] + p[3]*u[11] + p[4]*u[13]*u[5] + p[6]*u[2]*u[6] - p[5]*u[9] - p[7]*u[9] - p[1]*p[13]*u[1]*u[9] #* u[9]
#         du[10] = p[11]*u[12] + p[12]*u[12] + p[4]*u[14]*u[5] + p[8]*u[6]*u[3] - p[5]*u[10] - p[9]*u[10] - p[10]*p[13]*u[5]*u[10] #* u[10]
#         du[11] = p[1]*p[13]*u[1]*u[9] + p[4]*p[13]*u[15]*u[5] + p[6]*p[13]*u[7]*u[6] - p[2]*u[11] - p[5]*u[11] - p[7]*u[11] - p[3]*u[11] #* u[11]
#         du[12] = p[4]*p[13]*u[16]*u[5] + p[10]*p[13]*u[5]*u[10] + p[8]*p[13]*u[6]*u[8] - p[5]*u[12] - p[9]*u[12] - p[11]*u[12] - p[12]*u[12] #* u[12]
#         du[13] = p[2]*u[15] + p[5]*u[9] + p[3]*u[15] + p[6]*u[4]*u[2] - p[7]*u[13] - p[1]*u[13]*u[1] - p[4]*u[13]*u[5] #* u[13]
#         du[14] = p[5]*u[10] + p[11]*u[16] + p[12]*u[16] + p[8]*u[4]*u[3] - p[9]*u[14] - p[4]*u[14]*u[5] - p[10]*u[14]*u[5] #* u[14]
#         du[15] = p[5]*u[11] + p[6]*u[4]*u[7] + p[1]*u[13]*u[1] - p[2]*u[15] - p[7]*u[15] - p[3]*u[15] - p[4]*p[13]*u[15]*u[5] #* u[15]
#         du[16] = p[5]*u[12] + p[10]*u[14]*u[5] + p[8]*u[4]*u[8] - p[9]*u[16] - p[11]*u[16] - p[12]*u[16] - p[4]*p[13]*u[16]*u[5] #* APLp
#         nothing
# end

#* Out of place, with StaticArrays
function fullmodel(u, p, t)
        L, K, P, A, Lp, LpA, LK, LpP, LpAK, LpAP, LpAKL, LpAPLp, AK, AP, AKL, APLp = u #initial conditions
        ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, y = p #parameters
        dL = kb1*AKL + kcat7*APLp + kb1*LK + kb1*LpAKL + kcat7*LpAPLp + kcat7*LpP - ka1*AK*L - ka1*K*L - ka1*y*L*LpAK #* L
        dK = kb1*LK + kb3*AK + kcat1*LK + kb3*LpAK - ka3*A*K - ka1*K*L - ka3*K*LpA #* K
        dP = kb4*AP + kb4*LpAP + kb7*LpP + kcat7*LpP - ka4*A*P - ka7*Lp*P - ka4*LpA*P #* P
        dA = kb2*LpA + kb3*AK + kb3*AKL + kb4*AP + kb4*APLp - ka2*A*Lp - ka3*A*K - ka3*A*LK - ka4*A*LpP - ka4*A*P #* A
        dLp = kb7*APLp + kcat1*AKL + kcat1*LK + kb2*LpA + kb2*LpAK + kb2*LpAKL + kb2*LpAP + kcat1*LpAKL + kb2*LpAPLp + kb7*LpAPLp + kb7*LpP - ka2*A*Lp - ka2*AK*Lp - ka2*AP*Lp - ka7*AP*Lp - ka7*Lp*P - ka2*y*AKL*Lp - ka2*y*APLp*Lp - ka7*y*Lp*LpAP #* Lp
        dLpA = kb3*LpAK + kb3*LpAKL + kb4*LpAP + kb4*LpAPLp + ka2*A*Lp - kb2*LpA - ka3*K*LpA - ka4*LpA*P - ka3*y*LK*LpA - ka4*y*LpA*LpP #* LpA
        dLK = kb3*AKL + kb3*LpAKL + ka1*K*L - kb1*LK - kcat1*LK - ka3*A*LK - ka3*y*LK*LpA #* LK
        dLpP = kb4*APLp + kb4*LpAPLp + ka7*Lp*P - kb7*LpP - kcat7*LpP - ka4*A*LpP - ka4*y*LpA*LpP #* LpP
        dLpAK = kb1*LpAKL + kcat1*LpAKL + ka2*AK*Lp + ka3*K*LpA - kb2*LpAK - kb3*LpAK - ka1*y*L*LpAK #* LpAK
        dLpAP = kb7*LpAPLp + kcat7*LpAPLp + ka2*AP*Lp + ka4*LpA*P - kb2*LpAP - kb4*LpAP - ka7*y*Lp*LpAP #* LpAP
        dLpAKL = ka1*y*L*LpAK + ka2*y*AKL*Lp + ka3*y*LK*LpA - kb1*LpAKL - kb2*LpAKL - kb3*LpAKL - kcat1*LpAKL #* LpAKL
        dLpAPLp = ka2*y*APLp*Lp + ka7*y*Lp*LpAP + ka4*y*LpA*LpP - kb2*LpAPLp - kb4*LpAPLp - kb7*LpAPLp - kcat7*LpAPLp #* LpAPLp
        dAK = kb1*AKL + kb2*LpAK + kcat1*AKL + ka3*A*K - kb3*AK - ka1*AK*L - ka2*AK*Lp #* AK
        dAP = kb2*LpAP + kb7*APLp + kcat7*APLp + ka4*A*P - kb4*AP - ka2*AP*Lp - ka7*AP*Lp #* AP
        dAKL = kb2*LpAKL + ka3*A*LK + ka1*AK*L - kb1*AKL - kb3*AKL - kcat1*AKL - ka2*y*AKL*Lp #* AKL
        dAPLp = kb2*LpAPLp + ka7*AP*Lp + ka4*A*LpP - kb4*APLp - kb7*APLp - kcat7*APLp - ka2*y*APLp*Lp #* APLp
        SA[dL, dK, dP, dA, dLp, dLpA, dLK, dLpP, dLpAK, dLpAP, dLpAKL, dLpAPLp, dAK, dAP, dAKL, dAPLp]
end

function fullmodel_static(u, p, t)
        dL = p[2]*u[15] + p[12]*u[16] + p[2]*u[7] + p[2]*u[11] + p[12]*u[12] + p[12]*u[8] - p[1]*u[13]*u[1] - p[1]*u[2]*u[1] - p[1]*p[13]*u[1]*u[9] #* L
        dK = p[2]*u[7] + p[7]*u[13] + p[3]*u[7] + p[7]*u[9] - p[6]*u[4]*u[2] - p[1]*u[2]*u[1] - p[6]*u[2]*u[6] #* u[2]
        dP = p[9]*u[14] + p[9]*u[10] + p[11]*u[8] + p[12]*u[8] - p[8]*u[4]*u[3] - p[10]*u[5]*u[3] - p[8]*u[6]*u[3] #* u[3]
        dA = p[5]*u[6] + p[7]*u[13] + p[7]*u[15] + p[9]*u[14] + p[9]*u[16] - p[4]*u[4]*u[5] - p[6]*u[4]*u[2] - p[6]*u[4]*u[7] - p[8]*u[4]*u[8] - p[8]*u[4]*u[3] #* u[4]
        dLp = p[11]*u[16] + p[3]*u[15] + p[3]*u[7] + p[5]*u[6] + p[5]*u[9] + p[5]*u[11] + p[5]*u[10] + p[3]*u[11] + p[5]*u[12] + p[11]*u[12] + p[11]*u[8] - p[4]*u[4]*u[5] - p[4]*u[13]*u[5] - p[4]*u[14]*u[5] - p[10]*u[14]*u[5] - p[10]*u[5]*u[3] - p[4]*p[13]*u[15]*u[5] - p[4]*p[13]*u[16]*u[5] - p[10]*p[13]*u[5]*u[10] #* u[5]
        dLpA = p[7]*u[9] + p[7]*u[11] + p[9]*u[10] + p[9]*u[12] + p[4]*u[4]*u[5] - p[5]*u[6] - p[6]*u[2]*u[6] - p[8]*u[6]*u[3] - p[6]*p[13]*u[7]*u[6] - p[8]*p[13]*u[6]*u[8] #* u[6]
        dLK = p[7]*u[15] + p[7]*u[11] + p[1]*u[2]*u[1] - p[2]*u[7] - p[3]*u[7] - p[6]*u[4]*u[7] - p[6]*p[13]*u[7]*u[6] #* u[7]
        dLpP = p[9]*u[16] + p[9]*u[12] + p[10]*u[5]*u[3] - p[11]*u[8] - p[12]*u[8] - p[8]*u[4]*u[8] - p[8]*p[13]*u[6]*u[8] #* u[8]
        dLpAK = p[2]*u[11] + p[3]*u[11] + p[4]*u[13]*u[5] + p[6]*u[2]*u[6] - p[5]*u[9] - p[7]*u[9] - p[1]*p[13]*u[1]*u[9] #* u[9]
        dLpAP = p[11]*u[12] + p[12]*u[12] + p[4]*u[14]*u[5] + p[8]*u[6]*u[3] - p[5]*u[10] - p[9]*u[10] - p[10]*p[13]*u[5]*u[10] #* u[10]
        dLpAKL = p[1]*p[13]*u[1]*u[9] + p[4]*p[13]*u[15]*u[5] + p[6]*p[13]*u[7]*u[6] - p[2]*u[11] - p[5]*u[11] - p[7]*u[11] - p[3]*u[11] #* u[11]
        dLpAPLp = p[4]*p[13]*u[16]*u[5] + p[10]*p[13]*u[5]*u[10] + p[8]*p[13]*u[6]*u[8] - p[5]*u[12] - p[9]*u[12] - p[11]*u[12] - p[12]*u[12] #* u[12]
        dAK = p[2]*u[15] + p[5]*u[9] + p[3]*u[15] + p[6]*u[4]*u[2] - p[7]*u[13] - p[1]*u[13]*u[1] - p[4]*u[13]*u[5] #* u[13]
        dAP = p[5]*u[10] + p[11]*u[16] + p[12]*u[16] + p[8]*u[4]*u[3] - p[9]*u[14] - p[4]*u[14]*u[5] - p[10]*u[14]*u[5] #* u[14]
        dAKL = p[5]*u[11] + p[6]*u[4]*u[7] + p[1]*u[13]*u[1] - p[2]*u[15] - p[7]*u[15] - p[3]*u[15] - p[4]*p[13]*u[15]*u[5] #* u[15]
        dAPLp = p[5]*u[12] + p[10]*u[14]*u[5] + p[8]*u[4]*u[8] - p[9]*u[16] - p[11]*u[16] - p[12]*u[16] - p[4]*p[13]*u[16]*u[5] #* APLp
        SA[dL, dK, dP, dA, dLp, dLpA, dLK, dLpP, dLpAK, dLpAP, dLpAKL, dLpAPLp, dAK, dAP, dAKL, dAPLp]
end

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

pstatic = SA[0.009433439939827041, 2.3550169939427845, 832.7213093872278, 12.993995997539924, 6.150972501791291, 
                1.3481451097940793, 0.006201726090609513, 0.006277294665474662, 0.9250191811994848, 57.36471615394549, 0.04411989797898752, 42.288085868394326, 3631.050539219606]
u0static = SA[0.0, 0.5, 0.3, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

#? Timespan for integration
tspan = (0., 100.)

inplaceprob = ODEProblem(fullmodel_ode!, u0, tspan, p)
fullrn = make_fullrn()
ogprob = ODEProblem(fullrn, u0, tspan, p)
oopprob = ODEProblem(fullmodel_static, u0static, tspan, pstatic)

@btime solve($inplaceprob, saveat=0.1, save_idxs=1)
@btime solve($ogprob, saveat=0.1, save_idxs=1)
@btime solve($oopprob, saveat=0.1, save_idxs=1)




# fullmodel(u,p,t) = fullmodel!(similar(u),u,p,t)

# @code_warntype ODEProblem(fullmodel_ode!, u0, tspan, p)
odeprob = ODEProblem(fullmodel_ode!, u0, tspan, p)
# @code_warntype solve(odeprob, saveat=0.1, save_idxs=1)
# odeprobstatic = ODEProblem(fullmodel, u0static, tspan, p)
# fullprob = ODEProblem(fullrn, u0, tspan, p)

# @benchmark solve($odeprob, saveat=0.1, save_idxs=1)
# @benchmark solve($odeprobstatic, saveat=0.1, save_idxs=1)
# @benchmark solve($fullprob, saveat=0.1, save_idxs=1)