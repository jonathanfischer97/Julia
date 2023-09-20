using Plots 
using DifferentialEquations
using BifurcationKit, Setfield, LinearAlgebra, ForwardDiff, Parameters; const BK = BifurcationKit
using Peaks
using Statistics
default(lw = 2, size = (1000, 600))

norminf(x) = norm(x, Inf)

"""Calculates the period and amplitude of each individual in the population"""
function getPerAmp(sol::ODESolution, idx = 1)
    # Find peaks and calculate amplitudes and periods
    indx_max, vals_max = findmaxima(sol[idx,:], 1)
    indx_min, vals_min = findminima(sol[idx,:], 1)

    if length(indx_max) < 2 || length(indx_min) < 2
        return 0., 0.
    else
        # Calculate amplitudes and periods
        @inbounds amps = [(vals_max[i] - vals_min[i])/2 for i in 1:min(length(indx_max), length(indx_min))]
        @inbounds pers = [sol.t[indx_max[i+1]] - sol.t[indx_max[i]] for i in 1:(length(indx_max)-1)]

        # Calculate means of amplitudes and periods
        amp = mean(amps)
        per = mean(pers)

        return per, amp
    end
end


function fullmodel!(du, u, p, t=0)
    L, Lp, K, P, A, LpA, LK, LpP, LpAK, LpAP, LpAKL, LpAPLp, AK, AP, AKL, APLp = u #initial conditions
    @unpack ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, y = p #parameters
    du[1] = kb1*AKL + kcat7*APLp + kb1*LK + kb1*LpAKL + kcat7*LpAPLp + kcat7*LpP - ka1*AK*L - ka1*K*L - ka1*y*L*LpAK #* L
    du[2] = kb7*APLp + kcat1*AKL + kcat1*LK + kb2*LpA + kb2*LpAK + kb2*LpAKL + kb2*LpAP + kcat1*LpAKL + kb2*LpAPLp + kb7*LpAPLp + kb7*LpP - ka2*A*Lp - ka2*AK*Lp - ka2*AP*Lp - ka7*AP*Lp - ka7*Lp*P - ka2*y*AKL*Lp - ka2*y*APLp*Lp - ka7*y*Lp*LpAP #* Lp
    du[3] = kb1*LK + kb3*AK + kcat1*LK + kb3*LpAK - ka3*A*K - ka1*K*L - ka3*K*LpA #* K
    du[4] = kb4*AP + kb4*LpAP + kb7*LpP + kcat7*LpP - ka4*A*P - ka7*Lp*P - ka4*LpA*P #* P
    du[5] = kb2*LpA + kb3*AK + kb3*AKL + kb4*AP + kb4*APLp - ka2*A*Lp - ka3*A*K - ka3*A*LK - ka4*A*LpP - ka4*A*P #* A
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

fullmodel(u,p) = fullmodel!(similar(u),u,p,0)

parms = (ka1 = 0.009433439939827041, kb1 = 2.3550169939427845, kcat1 = 832.7213093872278, ka2 = 12.993995997539924, kb2 = 6.150972501791291,
        ka3 = 0.3481451097940793, kb3 = 0.006201726090609513, ka4 = 0.006277294665474662, kb4 = 0.9250191811994848, ka7 = 57.36471615394549, 
        kb7 = 0.04411989797898752, kcat7 = 42.288085868394326, y = 3631.050539219606)

u0 = [1.3943946650740193, 0.38465690613621295, 0.0016145597945468088, 0.06946855142927788, 0.7482934787301279, 
    0.555317742680265, 2.725521420682084e-8, 0.032528082931354375, 0.21314167678268764, 0.00010445506753742689, 
    0.012191541654353318, 0.19770110727469617, 0.27304805518766095, 0.00013079640806368744, 4.139325538860142e-6, 6.700688906891654e-5]

prob = BifurcationProblem(fullmodel, u0, parms, (@lens _.kcat7);
	recordFromSolution = (x, p) -> (L = x[1], Lp = x[2], K = x[3], P = x[4], A = x[5], LpA = x[6], LK = x[7], LpP = x[8], LpAK = x[9], LpAP = x[10], LpAKL = x[11], LpAPLp = x[12], AK = x[13], AP = x[14], AKL = x[15], APLp = x[16]))


# Adjust the Newton solver parameters
newton_params = NewtonPar(tol = 1e-10, maxIter = 1000)  # Increase maxIter and adjust the tolerance

# Update the opts_br for the BifurcationProblem
opts_br = ContinuationPar(θ=0.1, a=0.1,  pMin = 10.0, pMax = 5000.0, ds = 0.001, dsmax = 0.01, nInversion = 6, detectBifurcation = 3, maxBisectionSteps = 25, nev = 16, maxSteps = 200000, newtonOptions = newton_params)


#! Simulate the model
prob_de = ODEProblem(fullmodel!, u0, (0,10.), parms)
sol = solve(prob_de, Rodas5())
per, amp = getPerAmp(sol)
prob_de = ODEProblem(fullmodel!, u0, (0,per), parms, reltol = 1e-10, abstol = 1e-12)
sol = solve(prob_de, Rodas5())
plot(sol)

#! Helper functions to record and plot the periodic orbits
argspo = (recordFromSolution = (x, p) -> begin
		xtt = getPeriodicOrbit(p.prob, x, set(getParams(p.prob), getLens(p.prob), p.p))
		return (max = maximum(xtt[1,:]),
				min = minimum(xtt[1,:]),
				period = getPeriod(p.prob, x, set(getParams(p.prob), getLens(p.prob), p.p)))
	end,
	plotSolution = (x, p; k...) -> begin
		xtt = getPeriodicOrbit(p.prob, x, set(getParams(p.prob), getLens(p.prob), p.p))
		plot!(xtt.t, xtt[1,:]; label = "x", k...)
		plot!(xtt.t, xtt[2,:]; label = "y", k...)
		# plot!(br; subplot = 1, putspecialptlegend = false)
	end)

#! Compute periodic orbits with trapezoid method
# this is the function which builds probtrap from sol
probtrap, ci = generateCIProblem(PeriodicOrbitTrapProblem(M = 151;
	jacobian = :DenseAD, updateSectionEveryStep = 0),
	prob, sol, per)

opts_po_cont = setproperties(opts_br, maxSteps = 500, tolStability = 1e-8)
brpo_fold = continuation(probtrap, ci, PALC(), opts_po_cont;
	verbosity = 3, plot = true, normC = norminf,
	argspo...
	)

scene = plot(brpo_fold)


#! Parallel standard shooting
record_sh = recordFromSolution = (x, p) -> begin
		xtt = getPeriodicOrbit(p.prob, x, set(parms, p.prob.lens, p.p))
		return (max = maximum(xtt[1,:]),
				min = minimum(xtt[1,:]),
				period = getPeriod(p.prob, x, set(parms, p.prob.lens, p.p)))
	end
plot_sh  = (x, p; k...) -> begin
	xtt = getPeriodicOrbit(p.prob, x, set(parms, p.prob.lens, p.p))
	plot!(xtt.t, xtt[1,:]; label = "x", k...)
	plot!(xtt.t, xtt[2,:]; label = "y", k...)
	# plot!(br; subplot = 1, putspecialptlegend = false)
	end

probsh, cish = generateCIProblem( ShootingProblem(M=3),
    prob, prob_de, sol, per; alg = Rodas5())

opts_po_cont = setproperties(opts_br, maxSteps = 500, saveEigenvectors = true, detectLoop = true, tolStability = 1e-3)
br_fold_sh = continuation(probsh, cish, PALC(tangent = Bordered()), opts_po_cont;
	verbosity = 3, plot = true,
	recordFromSolution = record_sh,
	plotSolution = plot_sh,)

scene = plot(br_fold_sh)


#! Orthogonal collocation
# this is the function which builds probcoll from sol
probcoll, ci = generateCIProblem(PeriodicOrbitOCollProblem(26, 3; updateSectionEveryStep = 0),
	prob, sol, 3)

opts_po_cont = setproperties(opts_br, maxSteps = 50, saveEigenvectors = true, tolStability = 1e-8)
brpo_fold = continuation(probcoll, ci, PALC(), opts_po_cont;
	verbosity = 3, plot = true,
	argspo...
	)
scene = plot(brpo_fold)

