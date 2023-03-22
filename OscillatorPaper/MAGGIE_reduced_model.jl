using Plots 
using Catalyst
using DifferentialEquations
using Symbolics
using Latexify



redrn = @reaction_network redrn begin
    (ka1, kb1), L + K ↔ LK
    kcat1, LK → Lp + K
    (ka2, kb2), Lp + A ↔ LpA
    (ka3, kb3), A + K ↔ LpAK
    (ka4, kb4), A + P ↔ LpAP
    (ka7, kb7), Lp + P ↔ LpP
    kcat7, LpP → LpAP + P
end 

osys = convert(ODESystem, redrn)


function format_equations(eqs::Symbolics.Arr{Equation,1})
    for ode in eqs
        ode_str = string(ode)
        # ode_str = replace(ode_str, "~" => "=")
        ode_str = replace(ode_str, "(t)" => "")
        ode_str = replace(ode_str, "Differential(" => "d")
        ode_str = replace(ode_str, ") ~" => " =")
        println(ode_str)
    end
end

eqs = format_equations(osys.eqs)




function reduced_oscillator_odes!(dy, y, p, t)
    L, Lp, LpA, LpAP = y
    ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, DF, ktot, ptot, atot, liptot = p

    LK, LpAK, P, LpP, LpAPLp, LpAKL, A, K = calc_other_vars(y, p)

    dy[1] = kb1*LK + kb1*LpAKL + kcat7*LpAPLp + kcat7*LpP - ka1*K*L - ka1*DF*L*LpAK # L
    dy[2] = kcat1*LK + kb2*LpA + kcat1*LpAKL + kb7*LpAPLp + kb7*LpP - ka2*A*Lp - ka7*Lp*P - ka7*DF*Lp*LpAP # Lp
    dy[3] = kb3*LpAK + kb4*LpAP + ka2*A*Lp - kb2*LpA - ka3*K*LpA - ka4*LpA*P # LpA
    dy[4] = kb7*LpAPLp + kcat7*LpAPLp + ka4*LpA*P - kb4*LpAP - ka7*DF*Lp*LpAP # LpAP

function calc_other_vars(y, p)
    L, Lp, LpA, LpAP = y
    ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, DF, ktot, ptot, atot, liptot = p

    P = ((kb7 + kcat7) * (kb1 * kb3 * (L - liptot + Lp + LpA - LpAP + 2 * ptot) +
        kb3 * kcat1 * (L - liptot + Lp + LpA - LpAP + 2 * ptot) +
        ka1 * kb3 * L * (ktot + L - liptot + Lp + LpA - LpAP + 2 * ptot) +
        ka3 * kb1 * LpA * (ktot + L - liptot + Lp + LpA - LpAP + 2 * ptot) +
        ka3 * kcat1 * LpA * (ktot + L - liptot + Lp + LpA - LpAP + 2 * ptot) +
        DF * ka1 * ka3 * L * LpA * (2 * ktot + L - liptot + Lp + LpA - LpAP + 2 * ptot))) /
        ((2 * kb7 + 2 * kcat7 + ka7 * Lp) * (kb3 * (kcat1 + ka1 * L) +
        ka3 * (kcat1 + DF * ka1 * L) * LpA + kb1 * (kb3 + ka3 * LpA)))

    LpP = (ka7 * Lp * (kb1 * kb3 * (L - liptot + Lp + LpA - LpAP + 2 * ptot) +
        kb3 * kcat1 * (L - liptot + Lp + LpA - LpAP + 2 * ptot) +
        ka1 * kb3 *
    