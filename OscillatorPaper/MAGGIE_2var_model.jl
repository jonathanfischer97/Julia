using Plots 
using Catalyst
using DifferentialEquations
using Peaks
using Statistics
using BenchmarkTools

# %So we eliminated A, LpAKL, LpAPLp and K using mass conservation--this is
# %exact.
# %We eliminated LK and LpAK and P and LpP by 4 steady-state approxes. 
# %%set dLK=0 and dLpAKL=0, solve for LK and LpAK. 
# %set dLpP=0 and dK=0, solve for P and LpP

# %eliminate LpAP using dLpAP=0

# %try eliminating Lp using dLpAPLp/dt=0

# %L+K is reaction 1: ka1, kb1, kcat1
# %Lp+A is reaction 2: ka2, kb2
# %A+K is reaction 3: ka3, kb3
# %A+P is reaction 4: ka4, kb4
# %Lp+P is reaction 7: ka7, kb7, kcat7

function getLp(L, LpA, p, tots)
    # L, LpA = y
    ka1, kb1, kcat1, _, _, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, DF = p
    ktot, ptot, _, liptot = tots
    ka34 = ka3*ka4
    kb34 = kb3*kb4
    LpAsqred = LpA^2

        (kb1*kb34*kb7^2 + kb34*kb7^2*kcat1 + 2*kb1*kb34*kb7*kcat7 + 
        2*kb34*kb7*kcat1*kcat7 + kb1*kb34*kcat7^2 + kb34*kcat1*kcat7^2 + 
        ka7*kb1*kb34*kb7*L + ka1*kb34*kb7^2*L + ka7*kb34*kb7*kcat1*L + 
        ka7*kb1*kb34*kcat7*L + 2*ka1*kb34*kb7*kcat7*L + 
        ka7*kb34*kcat1*kcat7*L + ka1*kb34*kcat7^2*L + 
        ka1*ka7*kb34*kb7*ktot*L + ka1*ka7*kb34*kcat7*ktot*L + 
        ka1*ka7*kb34*kb7*L^2 + ka1*ka7*kb34*kcat7*L^2 - 
        ka7*kb1*kb34*kb7*liptot - ka7*kb34*kb7*kcat1*liptot - 
        ka7*kb1*kb34*kcat7*liptot - ka7*kb34*kcat1*kcat7*liptot - 
        ka1*ka7*kb34*kb7*L*liptot - ka1*ka7*kb34*kcat7*L*liptot + 
        ka7*kb1*kb34*kb7*LpA + ka4*kb1*kb3*kb7^2*LpA + 2*ka3*kb1*kb4*kb7^2*LpA + 
        ka7*kb34*kb7*kcat1*LpA + ka4*kb3*kb7^2*kcat1*LpA + 
        2*ka3*kb4*kb7^2*kcat1*LpA + ka7*kb1*kb34*kcat7*LpA + 
        2*ka4*kb1*kb3*kb7*kcat7*LpA + 4*ka3*kb1*kb4*kb7*kcat7*LpA + 
        ka7*kb34*kcat1*kcat7*LpA + 2*ka4*kb3*kb7*kcat1*kcat7*LpA + 
        4*ka3*kb4*kb7*kcat1*kcat7*LpA + ka4*kb1*kb3*kcat7^2*LpA + 
        2*ka3*kb1*kb4*kcat7^2*LpA + ka4*kb3*kcat1*kcat7^2*LpA + 
        2*ka3*kb4*kcat1*kcat7^2*LpA + ka3*ka7*kb1*kb4*kb7*ktot*LpA + 
        ka3*ka7*kb4*kb7*kcat1*ktot*LpA + ka3*ka7*kb1*kb4*kcat7*ktot*LpA + 
        ka3*ka7*kb4*kcat1*kcat7*ktot*LpA + DF*ka4*ka7*kb1*kb3*kb7*L*LpA + 
        2*ka3*ka7*kb1*kb4*kb7*L*LpA + ka1*ka7*kb34*kb7*L*LpA + 
        ka1*ka4*kb3*kb7^2*L*LpA + DF*ka4*ka7*kb3*kb7*kcat1*L*LpA + 
        2*ka3*ka7*kb4*kb7*kcat1*L*LpA + DF*ka4*ka7*kb1*kb3*kcat7*L*LpA + 
        2*ka3*ka7*kb1*kb4*kcat7*L*LpA + ka1*ka7*kb34*kcat7*L*LpA + 
        2*ka1*ka4*kb3*kb7*kcat7*L*LpA + DF*ka4*ka7*kb3*kcat1*kcat7*L*LpA + 
        2*ka3*ka7*kb4*kcat1*kcat7*L*LpA + ka1*ka4*kb3*kcat7^2*L*LpA + 
        DF*ka1*ka4*ka7*kb3*kb7*ktot*L*LpA + 2*DF*ka1*ka3*ka7*kb4*kb7*ktot*L*LpA + 
        DF*ka1*ka4*ka7*kb3*kcat7*ktot*L*LpA + 2*DF*ka1*ka3*ka7*kb4*kcat7*ktot*L*LpA + 
        DF*ka1*ka4*ka7*kb3*kb7*L^2*LpA + DF*ka1*ka4*ka7*kb3*kcat7*L^2*LpA - 
        DF*ka4*ka7*kb1*kb3*kb7*liptot*LpA - 2*ka3*ka7*kb1*kb4*kb7*liptot*LpA - 
        DF*ka4*ka7*kb3*kb7*kcat1*liptot*LpA - 2*ka3*ka7*kb4*kb7*kcat1*liptot*LpA - 
        DF*ka4*ka7*kb1*kb3*kcat7*liptot*LpA - 2*ka3*ka7*kb1*kb4*kcat7*liptot*LpA - 
        DF*ka4*ka7*kb3*kcat1*kcat7*liptot*LpA - 2*ka3*ka7*kb4*kcat1*kcat7*liptot*LpA - 
        DF*ka1*ka4*ka7*kb3*kb7*L*liptot*LpA - DF*ka1*ka4*ka7*kb3*kcat7*L*liptot*LpA + 
        DF*ka4*ka7*kb1*kb3*kb7*LpAsqred + 2*ka3*ka7*kb1*kb4*kb7*LpAsqred + 
        2*ka34*kb1*kb7^2*LpAsqred + DF*ka4*ka7*kb3*kb7*kcat1*LpAsqred + 
        2*ka3*ka7*kb4*kb7*kcat1*LpAsqred + 2*ka34*kb7^2*kcat1*LpAsqred + 
        DF*ka4*ka7*kb1*kb3*kcat7*LpAsqred + 2*ka3*ka7*kb1*kb4*kcat7*LpAsqred + 
        4*ka34*kb1*kb7*kcat7*LpAsqred + DF*ka4*ka7*kb3*kcat1*kcat7*LpAsqred + 
        2*ka3*ka7*kb4*kcat1*kcat7*LpAsqred + 4*ka34*kb7*kcat1*kcat7*LpAsqred + 
        2*ka34*kb1*kcat7^2*LpAsqred + 2*ka34*kcat1*kcat7^2*LpAsqred + 
        DF*ka34*ka7*kb1*kb7*ktot*LpAsqred + DF*ka34*ka7*kb7*kcat1*ktot*LpAsqred + 
        DF*ka34*ka7*kb1*kcat7*ktot*LpAsqred + 
        DF*ka34*ka7*kcat1*kcat7*ktot*LpAsqred + 2*DF*ka34*ka7*kb1*kb7*L*LpAsqred + 
        DF*ka1*ka4*ka7*kb3*kb7*L*LpAsqred + 2*DF*ka34*ka7*kb7*kcat1*L*LpAsqred + 
        2*DF*ka34*ka7*kb1*kcat7*L*LpAsqred + DF*ka1*ka4*ka7*kb3*kcat7*L*LpAsqred + 
        2*DF*ka34*ka7*kcat1*kcat7*L*LpAsqred + 
        2*DF^2*ka1*ka34*ka7*kb7*ktot*L*LpAsqred + 
        2*DF^2*ka1*ka34*ka7*kcat7*ktot*L*LpAsqred - 
        2*DF*ka34*ka7*kb1*kb7*liptot*LpAsqred - 
        2*DF*ka34*ka7*kb7*kcat1*liptot*LpAsqred - 
        2*DF*ka34*ka7*kb1*kcat7*liptot*LpAsqred - 
        2*DF*ka34*ka7*kcat1*kcat7*liptot*LpAsqred + 2*DF*ka34*ka7*kb1*kb7*LpA^3 + 
        2*DF*ka34*ka7*kb7*kcat1*LpA^3 + 2*DF*ka34*ka7*kb1*kcat7*LpA^3 + 
        2*DF*ka34*ka7*kcat1*kcat7*LpA^3 + ka7*kb1*kb34*kb7*ptot + 
        ka7*kb34*kb7*kcat1*ptot + ka7*kb1*kb34*kcat7*ptot + 
        ka7*kb34*kcat1*kcat7*ptot + ka1*ka7*kb34*kb7*L*ptot + 
        ka1*ka7*kb34*kcat7*L*ptot + 2*DF*ka4*ka7*kb1*kb3*kb7*LpA*ptot + 
        2*ka3*ka7*kb1*kb4*kb7*LpA*ptot + 2*DF*ka4*ka7*kb3*kb7*kcat1*LpA*ptot + 
        2*ka3*ka7*kb4*kb7*kcat1*LpA*ptot + 2*DF*ka4*ka7*kb1*kb3*kcat7*LpA*ptot + 
        2*ka3*ka7*kb1*kb4*kcat7*LpA*ptot + 2*DF*ka4*ka7*kb3*kcat1*kcat7*LpA*ptot + 
        2*ka3*ka7*kb4*kcat1*kcat7*LpA*ptot + 2*DF*ka1*ka4*ka7*kb3*kb7*L*LpA*ptot + 
        2*DF*ka1*ka4*ka7*kb3*kcat7*L*LpA*ptot + 4*DF*ka34*ka7*kb1*kb7*LpAsqred*ptot + 
        4*DF*ka34*ka7*kb7*kcat1*LpAsqred*ptot + 
        4*DF*ka34*ka7*kb1*kcat7*LpAsqred*ptot + 
        4*DF*ka34*ka7*kcat1*kcat7*LpAsqred*ptot - 
        sqrt((-(kb1*kb34*kb7^2) - kb34*kb7^2*kcat1 - 2*kb1*kb34*kb7*kcat7 - 
              2*kb34*kb7*kcat1*kcat7 - kb1*kb34*kcat7^2 - 
              kb34*kcat1*kcat7^2 - ka7*kb1*kb34*kb7*L - ka1*kb34*kb7^2*L - 
              ka7*kb34*kb7*kcat1*L - ka7*kb1*kb34*kcat7*L - 
              2*ka1*kb34*kb7*kcat7*L - ka7*kb34*kcat1*kcat7*L - 
              ka1*kb34*kcat7^2*L - ka1*ka7*kb34*kb7*ktot*L - 
              ka1*ka7*kb34*kcat7*ktot*L - ka1*ka7*kb34*kb7*L^2 - 
              ka1*ka7*kb34*kcat7*L^2 + ka7*kb1*kb34*kb7*liptot + 
              ka7*kb34*kb7*kcat1*liptot + ka7*kb1*kb34*kcat7*liptot + 
              ka7*kb34*kcat1*kcat7*liptot + ka1*ka7*kb34*kb7*L*liptot + 
              ka1*ka7*kb34*kcat7*L*liptot - ka7*kb1*kb34*kb7*LpA - 
              ka4*kb1*kb3*kb7^2*LpA - 2*ka3*kb1*kb4*kb7^2*LpA - 
              ka7*kb34*kb7*kcat1*LpA - ka4*kb3*kb7^2*kcat1*LpA - 
              2*ka3*kb4*kb7^2*kcat1*LpA - ka7*kb1*kb34*kcat7*LpA - 
              2*ka4*kb1*kb3*kb7*kcat7*LpA - 4*ka3*kb1*kb4*kb7*kcat7*LpA - 
              ka7*kb34*kcat1*kcat7*LpA - 2*ka4*kb3*kb7*kcat1*kcat7*LpA - 
              4*ka3*kb4*kb7*kcat1*kcat7*LpA - ka4*kb1*kb3*kcat7^2*LpA - 
              2*ka3*kb1*kb4*kcat7^2*LpA - ka4*kb3*kcat1*kcat7^2*LpA - 
              2*ka3*kb4*kcat1*kcat7^2*LpA - ka3*ka7*kb1*kb4*kb7*ktot*LpA - 
              ka3*ka7*kb4*kb7*kcat1*ktot*LpA - ka3*ka7*kb1*kb4*kcat7*ktot*LpA - 
              ka3*ka7*kb4*kcat1*kcat7*ktot*LpA - DF*ka4*ka7*kb1*kb3*kb7*L*LpA - 
              2*ka3*ka7*kb1*kb4*kb7*L*LpA - ka1*ka7*kb34*kb7*L*LpA - 
              ka1*ka4*kb3*kb7^2*L*LpA - DF*ka4*ka7*kb3*kb7*kcat1*L*LpA - 
              2*ka3*ka7*kb4*kb7*kcat1*L*LpA - DF*ka4*ka7*kb1*kb3*kcat7*L*LpA - 
              2*ka3*ka7*kb1*kb4*kcat7*L*LpA - ka1*ka7*kb34*kcat7*L*LpA - 
              2*ka1*ka4*kb3*kb7*kcat7*L*LpA - DF*ka4*ka7*kb3*kcat1*kcat7*L*LpA - 
              2*ka3*ka7*kb4*kcat1*kcat7*L*LpA - ka1*ka4*kb3*kcat7^2*L*LpA - 
              DF*ka1*ka4*ka7*kb3*kb7*ktot*L*LpA - 2*DF*ka1*ka3*ka7*kb4*kb7*ktot*L*LpA - 
              DF*ka1*ka4*ka7*kb3*kcat7*ktot*L*LpA - 
              2*DF*ka1*ka3*ka7*kb4*kcat7*ktot*L*LpA - DF*ka1*ka4*ka7*kb3*kb7*L^2*LpA - 
              DF*ka1*ka4*ka7*kb3*kcat7*L^2*LpA + DF*ka4*ka7*kb1*kb3*kb7*liptot*LpA + 
              2*ka3*ka7*kb1*kb4*kb7*liptot*LpA + DF*ka4*ka7*kb3*kb7*kcat1*liptot*LpA + 
              2*ka3*ka7*kb4*kb7*kcat1*liptot*LpA + 
              DF*ka4*ka7*kb1*kb3*kcat7*liptot*LpA + 
              2*ka3*ka7*kb1*kb4*kcat7*liptot*LpA + 
              DF*ka4*ka7*kb3*kcat1*kcat7*liptot*LpA + 
              2*ka3*ka7*kb4*kcat1*kcat7*liptot*LpA + 
              DF*ka1*ka4*ka7*kb3*kb7*L*liptot*LpA + 
              DF*ka1*ka4*ka7*kb3*kcat7*L*liptot*LpA - DF*ka4*ka7*kb1*kb3*kb7*LpAsqred - 
              2*ka3*ka7*kb1*kb4*kb7*LpAsqred - 2*ka34*kb1*kb7^2*LpAsqred - 
              DF*ka4*ka7*kb3*kb7*kcat1*LpAsqred - 2*ka3*ka7*kb4*kb7*kcat1*LpAsqred - 
              2*ka34*kb7^2*kcat1*LpAsqred - DF*ka4*ka7*kb1*kb3*kcat7*LpAsqred - 
              2*ka3*ka7*kb1*kb4*kcat7*LpAsqred - 4*ka34*kb1*kb7*kcat7*LpAsqred - 
              DF*ka4*ka7*kb3*kcat1*kcat7*LpAsqred - 2*ka3*ka7*kb4*kcat1*kcat7*LpAsqred - 
              4*ka34*kb7*kcat1*kcat7*LpAsqred - 2*ka34*kb1*kcat7^2*LpAsqred - 
              2*ka34*kcat1*kcat7^2*LpAsqred - DF*ka34*ka7*kb1*kb7*ktot*LpAsqred - 
              DF*ka34*ka7*kb7*kcat1*ktot*LpAsqred - 
              DF*ka34*ka7*kb1*kcat7*ktot*LpAsqred - 
              DF*ka34*ka7*kcat1*kcat7*ktot*LpAsqred - 
              2*DF*ka34*ka7*kb1*kb7*L*LpAsqred - DF*ka1*ka4*ka7*kb3*kb7*L*LpAsqred - 
              2*DF*ka34*ka7*kb7*kcat1*L*LpAsqred - 
              2*DF*ka34*ka7*kb1*kcat7*L*LpAsqred - DF*ka1*ka4*ka7*kb3*kcat7*L*LpAsqred - 
              2*DF*ka34*ka7*kcat1*kcat7*L*LpAsqred - 
              2*DF^2*ka1*ka34*ka7*kb7*ktot*L*LpAsqred - 
              2*DF^2*ka1*ka34*ka7*kcat7*ktot*L*LpAsqred + 
              2*DF*ka34*ka7*kb1*kb7*liptot*LpAsqred + 
              2*DF*ka34*ka7*kb7*kcat1*liptot*LpAsqred + 
              2*DF*ka34*ka7*kb1*kcat7*liptot*LpAsqred + 
              2*DF*ka34*ka7*kcat1*kcat7*liptot*LpAsqred - 
              2*DF*ka34*ka7*kb1*kb7*LpA^3 - 2*DF*ka34*ka7*kb7*kcat1*LpA^3 - 
              2*DF*ka34*ka7*kb1*kcat7*LpA^3 - 2*DF*ka34*ka7*kcat1*kcat7*LpA^3 - 
              ka7*kb1*kb34*kb7*ptot - ka7*kb34*kb7*kcat1*ptot - 
              ka7*kb1*kb34*kcat7*ptot - ka7*kb34*kcat1*kcat7*ptot - 
              ka1*ka7*kb34*kb7*L*ptot - ka1*ka7*kb34*kcat7*L*ptot - 
              2*DF*ka4*ka7*kb1*kb3*kb7*LpA*ptot - 2*ka3*ka7*kb1*kb4*kb7*LpA*ptot - 
              2*DF*ka4*ka7*kb3*kb7*kcat1*LpA*ptot - 2*ka3*ka7*kb4*kb7*kcat1*LpA*ptot - 
              2*DF*ka4*ka7*kb1*kb3*kcat7*LpA*ptot - 2*ka3*ka7*kb1*kb4*kcat7*LpA*ptot - 
              2*DF*ka4*ka7*kb3*kcat1*kcat7*LpA*ptot - 
              2*ka3*ka7*kb4*kcat1*kcat7*LpA*ptot - 
              2*DF*ka1*ka4*ka7*kb3*kb7*L*LpA*ptot - 
              2*DF*ka1*ka4*ka7*kb3*kcat7*L*LpA*ptot - 
              4*DF*ka34*ka7*kb1*kb7*LpAsqred*ptot - 
              4*DF*ka34*ka7*kb7*kcat1*LpAsqred*ptot - 
              4*DF*ka34*ka7*kb1*kcat7*LpAsqred*ptot - 
              4*DF*ka34*ka7*kcat1*kcat7*LpAsqred*ptot)^2 - 
           4*(-(ka7*kb1*kb34*kb7) - ka7*kb34*kb7*kcat1 - ka7*kb1*kb34*kcat7 - 
              ka7*kb34*kcat1*kcat7 - ka1*ka7*kb34*kb7*L - 
              ka1*ka7*kb34*kcat7*L - DF*ka4*ka7*kb1*kb3*kb7*LpA - 
              2*ka3*ka7*kb1*kb4*kb7*LpA - DF*ka4*ka7*kb3*kb7*kcat1*LpA - 
              2*ka3*ka7*kb4*kb7*kcat1*LpA - DF*ka4*ka7*kb1*kb3*kcat7*LpA - 
              2*ka3*ka7*kb1*kb4*kcat7*LpA - DF*ka4*ka7*kb3*kcat1*kcat7*LpA - 
              2*ka3*ka7*kb4*kcat1*kcat7*LpA - DF*ka1*ka4*ka7*kb3*kb7*L*LpA - 
              DF*ka1*ka4*ka7*kb3*kcat7*L*LpA - 2*DF*ka34*ka7*kb1*kb7*LpAsqred - 
              2*DF*ka34*ka7*kb7*kcat1*LpAsqred - 2*DF*ka34*ka7*kb1*kcat7*LpAsqred - 
              2*DF*ka34*ka7*kcat1*kcat7*LpAsqred)*
            (-(kb1*kb34*kb7^2*L) - kb34*kb7^2*kcat1*L - 
              2*kb1*kb34*kb7*kcat7*L - 2*kb34*kb7*kcat1*kcat7*L - 
              kb1*kb34*kcat7^2*L - kb34*kcat1*kcat7^2*L - 
              ka1*kb34*kb7^2*ktot*L - 2*ka1*kb34*kb7*kcat7*ktot*L - 
              ka1*kb34*kcat7^2*ktot*L - ka1*kb34*kb7^2*L^2 - 
              2*ka1*kb34*kb7*kcat7*L^2 - ka1*kb34*kcat7^2*L^2 + 
              kb1*kb34*kb7^2*liptot + kb34*kb7^2*kcat1*liptot + 
              2*kb1*kb34*kb7*kcat7*liptot + 2*kb34*kb7*kcat1*kcat7*liptot + 
              kb1*kb34*kcat7^2*liptot + kb34*kcat1*kcat7^2*liptot + 
              ka1*kb34*kb7^2*L*liptot + 2*ka1*kb34*kb7*kcat7*L*liptot + 
              ka1*kb34*kcat7^2*L*liptot - kb1*kb34*kb7^2*LpA - 
              kb34*kb7^2*kcat1*LpA - 2*kb1*kb34*kb7*kcat7*LpA - 
              2*kb34*kb7*kcat1*kcat7*LpA - kb1*kb34*kcat7^2*LpA - 
              kb34*kcat1*kcat7^2*LpA - ka3*kb1*kb4*kb7^2*ktot*LpA - 
              ka3*kb4*kb7^2*kcat1*ktot*LpA - 2*ka3*kb1*kb4*kb7*kcat7*ktot*LpA - 
              2*ka3*kb4*kb7*kcat1*kcat7*ktot*LpA - ka3*kb1*kb4*kcat7^2*ktot*LpA - 
              ka3*kb4*kcat1*kcat7^2*ktot*LpA - ka4*kb1*kb3*kb7^2*L*LpA - 
              2*ka3*kb1*kb4*kb7^2*L*LpA - ka1*kb34*kb7^2*L*LpA - 
              ka4*kb3*kb7^2*kcat1*L*LpA - 2*ka3*kb4*kb7^2*kcat1*L*LpA - 
              2*ka4*kb1*kb3*kb7*kcat7*L*LpA - 4*ka3*kb1*kb4*kb7*kcat7*L*LpA - 
              2*ka1*kb34*kb7*kcat7*L*LpA - 2*ka4*kb3*kb7*kcat1*kcat7*L*LpA - 
              4*ka3*kb4*kb7*kcat1*kcat7*L*LpA - ka4*kb1*kb3*kcat7^2*L*LpA - 
              2*ka3*kb1*kb4*kcat7^2*L*LpA - ka1*kb34*kcat7^2*L*LpA - 
              ka4*kb3*kcat1*kcat7^2*L*LpA - 2*ka3*kb4*kcat1*kcat7^2*L*LpA - 
              ka1*ka4*kb3*kb7^2*ktot*L*LpA - 2*DF*ka1*ka3*kb4*kb7^2*ktot*L*LpA - 
              2*ka1*ka4*kb3*kb7*kcat7*ktot*L*LpA - 
              4*DF*ka1*ka3*kb4*kb7*kcat7*ktot*L*LpA - ka1*ka4*kb3*kcat7^2*ktot*L*LpA - 
              2*DF*ka1*ka3*kb4*kcat7^2*ktot*L*LpA - ka1*ka4*kb3*kb7^2*L^2*LpA - 
              2*ka1*ka4*kb3*kb7*kcat7*L^2*LpA - ka1*ka4*kb3*kcat7^2*L^2*LpA + 
              ka4*kb1*kb3*kb7^2*liptot*LpA + 2*ka3*kb1*kb4*kb7^2*liptot*LpA + 
              ka4*kb3*kb7^2*kcat1*liptot*LpA + 2*ka3*kb4*kb7^2*kcat1*liptot*LpA + 
              2*ka4*kb1*kb3*kb7*kcat7*liptot*LpA + 4*ka3*kb1*kb4*kb7*kcat7*liptot*LpA + 
              2*ka4*kb3*kb7*kcat1*kcat7*liptot*LpA + 
              4*ka3*kb4*kb7*kcat1*kcat7*liptot*LpA + ka4*kb1*kb3*kcat7^2*liptot*LpA + 
              2*ka3*kb1*kb4*kcat7^2*liptot*LpA + ka4*kb3*kcat1*kcat7^2*liptot*LpA + 
              2*ka3*kb4*kcat1*kcat7^2*liptot*LpA + ka1*ka4*kb3*kb7^2*L*liptot*LpA + 
              2*ka1*ka4*kb3*kb7*kcat7*L*liptot*LpA + 
              ka1*ka4*kb3*kcat7^2*L*liptot*LpA - ka4*kb1*kb3*kb7^2*LpAsqred - 
              2*ka3*kb1*kb4*kb7^2*LpAsqred - ka4*kb3*kb7^2*kcat1*LpAsqred - 
              2*ka3*kb4*kb7^2*kcat1*LpAsqred - 2*ka4*kb1*kb3*kb7*kcat7*LpAsqred - 
              4*ka3*kb1*kb4*kb7*kcat7*LpAsqred - 2*ka4*kb3*kb7*kcat1*kcat7*LpAsqred - 
              4*ka3*kb4*kb7*kcat1*kcat7*LpAsqred - ka4*kb1*kb3*kcat7^2*LpAsqred - 
              2*ka3*kb1*kb4*kcat7^2*LpAsqred - ka4*kb3*kcat1*kcat7^2*LpAsqred - 
              2*ka3*kb4*kcat1*kcat7^2*LpAsqred - ka34*kb1*kb7^2*ktot*LpAsqred - 
              ka34*kb7^2*kcat1*ktot*LpAsqred - 2*ka34*kb1*kb7*kcat7*ktot*LpAsqred - 
              2*ka34*kb7*kcat1*kcat7*ktot*LpAsqred - 
              ka34*kb1*kcat7^2*ktot*LpAsqred - ka34*kcat1*kcat7^2*ktot*LpAsqred - 
              2*ka34*kb1*kb7^2*L*LpAsqred - ka1*ka4*kb3*kb7^2*L*LpAsqred - 
              2*ka34*kb7^2*kcat1*L*LpAsqred - 4*ka34*kb1*kb7*kcat7*L*LpAsqred - 
              2*ka1*ka4*kb3*kb7*kcat7*L*LpAsqred - 4*ka34*kb7*kcat1*kcat7*L*LpAsqred - 
              2*ka34*kb1*kcat7^2*L*LpAsqred - ka1*ka4*kb3*kcat7^2*L*LpAsqred - 
              2*ka34*kcat1*kcat7^2*L*LpAsqred - 
              2*DF*ka1*ka34*kb7^2*ktot*L*LpAsqred - 
              4*DF*ka1*ka34*kb7*kcat7*ktot*L*LpAsqred - 
              2*DF*ka1*ka34*kcat7^2*ktot*L*LpAsqred + 
              2*ka34*kb1*kb7^2*liptot*LpAsqred + 
              2*ka34*kb7^2*kcat1*liptot*LpAsqred + 
              4*ka34*kb1*kb7*kcat7*liptot*LpAsqred + 
              4*ka34*kb7*kcat1*kcat7*liptot*LpAsqred + 
              2*ka34*kb1*kcat7^2*liptot*LpAsqred + 
              2*ka34*kcat1*kcat7^2*liptot*LpAsqred - 2*ka34*kb1*kb7^2*LpA^3 - 
              2*ka34*kb7^2*kcat1*LpA^3 - 4*ka34*kb1*kb7*kcat7*LpA^3 - 
              4*ka34*kb7*kcat1*kcat7*LpA^3 - 2*ka34*kb1*kcat7^2*LpA^3 - 
              2*ka34*kcat1*kcat7^2*LpA^3 - ka4*kb1*kb3*kb7^2*LpA*ptot - 
              ka4*kb3*kb7^2*kcat1*LpA*ptot - 2*ka4*kb1*kb3*kb7*kcat7*LpA*ptot - 
              2*ka4*kb3*kb7*kcat1*kcat7*LpA*ptot - ka4*kb1*kb3*kcat7^2*LpA*ptot - 
              ka4*kb3*kcat1*kcat7^2*LpA*ptot - ka1*ka4*kb3*kb7^2*L*LpA*ptot - 
              2*ka1*ka4*kb3*kb7*kcat7*L*LpA*ptot - ka1*ka4*kb3*kcat7^2*L*LpA*ptot - 
              2*ka34*kb1*kb7^2*LpAsqred*ptot - 2*ka34*kb7^2*kcat1*LpAsqred*ptot - 
              4*ka34*kb1*kb7*kcat7*LpAsqred*ptot - 
              4*ka34*kb7*kcat1*kcat7*LpAsqred*ptot - 
              2*ka34*kb1*kcat7^2*LpAsqred*ptot - 2*ka34*kcat1*kcat7^2*LpAsqred*ptot)
            ))/(2*(-(ka7*kb1*kb34*kb7) - ka7*kb34*kb7*kcat1 - 
            ka7*kb1*kb34*kcat7 - ka7*kb34*kcat1*kcat7 - ka1*ka7*kb34*kb7*L - 
            ka1*ka7*kb34*kcat7*L - DF*ka4*ka7*kb1*kb3*kb7*LpA - 
            2*ka3*ka7*kb1*kb4*kb7*LpA - DF*ka4*ka7*kb3*kb7*kcat1*LpA - 
            2*ka3*ka7*kb4*kb7*kcat1*LpA - DF*ka4*ka7*kb1*kb3*kcat7*LpA - 
            2*ka3*ka7*kb1*kb4*kcat7*LpA - DF*ka4*ka7*kb3*kcat1*kcat7*LpA - 
            2*ka3*ka7*kb4*kcat1*kcat7*LpA - DF*ka1*ka4*ka7*kb3*kb7*L*LpA - 
            DF*ka1*ka4*ka7*kb3*kcat7*L*LpA - 2*DF*ka34*ka7*kb1*kb7*LpAsqred - 
            2*DF*ka34*ka7*kb7*kcat1*LpAsqred - 2*DF*ka34*ka7*kb1*kcat7*LpAsqred - 
		    2*DF*ka34*ka7*kcat1*kcat7*LpAsqred))
end

"""Given the unknowns, calculate the dependent variables"""
function calc_other_vars(y, p, tots)
    L, LpA = y
    ka1, kb1, kcat1, _, _, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, DF = p
    ktot, ptot, atot, liptot = tots

    Lp = getLp(L, LpA, p, tots)

    LpAP= -(((kb7 + kcat7)*(kb3*(kb7 + kcat7 + ka7*Lp - ka4*LpA)*
            (kcat1*(L - liptot + Lp + LpA) + ka1*L*(ktot + L - liptot + Lp + LpA)) + 
            ka3*LpA*(kb7 + kcat7 + ka7*Lp - ka4*LpA)*
            (2*DF*ka1*ktot*L + kcat1*(ktot + 2*(L - liptot + Lp + LpA))) + 
            kb1*(kb7 + kcat7 + ka7*Lp - ka4*LpA)*
            (kb3*(L - liptot + Lp + LpA) + ka3*LpA*(ktot + 2*(L - liptot + Lp + LpA)))
            + kb3*(kcat1 + ka1*L)*(ka7*Lp - 2*ka4*LpA)*ptot + 
            2*ka3*kcat1*LpA*(ka7*Lp - 2*ka4*LpA)*ptot + 
            kb1*(kb3 + 2*ka3*LpA)*(ka7*Lp - 2*ka4*LpA)*ptot))/
            ((kb3*(kb1 + kcat1 + ka1*L) + 2*ka3*(kb1 + kcat1)*LpA)*
            (kb7^2 + 2*kb7*kcat7 + kcat7^2 + 2*kb4*(kb7 + kcat7) + ka7*kb4*Lp + 
            2*DF*ka7*kb7*Lp + 2*DF*ka7*kcat7*Lp + DF*ka7^2*Lp^2 + 
            ka4*(kb7 + kcat7)*LpA)))


    common_sum = kb1 * kb3 * (L - liptot + Lp + LpA - LpAP + 2 * ptot) +
                 kb3 * kcat1 * (L - liptot + Lp + LpA - LpAP + 2 * ptot) +
                 ka1 * kb3 * L * (ktot + L - liptot + Lp + LpA - LpAP + 2 * ptot) +
                 ka3 * kb1 * LpA * (ktot + L - liptot + Lp + LpA - LpAP + 2 * ptot) +
                 ka3 * kcat1 * LpA * (ktot + L - liptot + Lp + LpA - LpAP + 2 * ptot) +
                 DF * ka1 * ka3 * L * LpA * (2 * ktot + L - liptot + Lp + LpA - LpAP + 2 * ptot)

    denominator = (2 * kb7 + 2 * kcat7 + ka7 * Lp) * (kb3 * (kcat1 + ka1 * L) + ka3 * (kcat1 + DF * ka1 * L) * LpA + kb1 * (kb3 + ka3 * LpA))


    P = ((kb7 + kcat7) * common_sum) / denominator

    LpP = ka7 * Lp * common_sum / denominator

    LK_denom = (kb1 + kcat1)^2 + 2 * DF * ka1 * (kb1 + kcat1) * L + DF * ka1^2 * L^2
    LK = ka1 * L * (kb1 * (ktot + L - liptot + Lp + LpA - LpAP - LpP - 2 * P + 2 * ptot) + kcat1 * (ktot + L - liptot + Lp + LpA - LpAP - LpP - 2 * P + 2 * ptot) + DF * ka1 * L * 
            (2 * ktot + L - liptot + Lp + LpA - LpAP - LpP - 2 * P + 2 * ptot)) / LK_denom

    LpAK = -((kb1 + kcat1) * (kb1 * (L - liptot + Lp + LpA - LpAP - LpP - 2 * P + 2 * ptot) + kcat1 * (L - liptot + Lp + LpA - LpAP - LpP - 2 * P + 2 * ptot) + ka1 * L * 
            (ktot + L - liptot + Lp + LpA - LpAP - LpP - 2 * P + 2 * ptot))) / LK_denom


    #mass constraints 
    LpAPLp= -LpAP - LpP - P + ptot
    LpAKL = 0.5*(-L + liptot - LK - Lp - LpA - LpAK + LpAP + LpP + 2*P - 2*ptot)
    A= 0.5*(2*atot + L - liptot + LK + Lp - LpA - LpAK - LpAP + LpP)
    K = 0.5* (2*ktot + L - liptot - LK + Lp + LpA - LpAK - LpAP - LpP - 2*P + 2*ptot)

    return Lp, LpAP, P, LpP, LK, LpAK, LpAPLp, LpAKL, A, K
end

"""ODE function for the reduced 2nd order oscillator model"""
function reduced_oscillator_odes!(dy, y, p, t)
    L, LpA = y #variables in the model
    ka1, kb1, _, ka2, kb2, ka3, kb3, ka4, kb4, _, _, kcat7, DF = p #parameters in the model
    tots = p[end-3:end] #total concentrations of the species as constraints

    Lp, LpAP, P, LpP, LK, LpAK, LpAPLp, LpAKL, A, K = calc_other_vars(y, p, tots) #calculate other variables as "constants and apply constraints"

    dy[1] = kb1*LK + kb1*LpAKL + kcat7*LpAPLp + kcat7*LpP - ka1*K*L - ka1*DF*L*LpAK # L
    dy[2] = kb3*LpAK + kb4*LpAP + ka2*A*Lp - kb2*LpA - ka3*K*LpA - ka4*LpA*P # LpA
end




"""Full oscillator model for comparison"""
osc_rn = @reaction_network osc_rn begin
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
end

    
#parameter list
"""ka1, kb1, kcat1, ka2, kb2, ka3, kb3, ka4, kb4, ka7, kb7, kcat7, DF"""
p = [0.125, 411.43070022437774, 42.44753785470635, 57.56861328754449, 1.0, 0.75, 0.001, 0.126, 
    102.45983604574394, 123.39827909372974, 25.805378439197266, 470.69414436040074, 3718.0197650684563]

#initial condition list
# L, Lp, K, P, A, LpA, LK, LpAK, LpAKL, LpP, LpAP, LpAPLp = u0

# u0 = [0., 3.0, 0.2, 0.3, 0.2, 0.0, 0., 0., 0., 0., 0., 0.] #working 
u0 = [1.01, 3.0, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 1.6, 0.3, 0.2]
L, Lp, LpA, LpAP, LK, LpAK, LpAKL, LpP, LpAPLp, A, K, P = u0

umap = [:L => L, :Lp => Lp, :LpA => LpA, :LpAP => LpAP, :LK => LK, :LpAK => LpAK, :LpAKL => LpAKL, :LpP => LpP, :LpAPLp => LpAPLp, :A => A, :K => K, :P => P]

#conserved quantities
ktot=K+LK+LpAK+LpAKL;
ptot=P+LpP+LpAP+LpAPLp;
atot=A+LpA+LpAK+LpAKL+LpAP+LpAPLp;
ltot=L+LK+LpAKL;
lptot=Lp+LpA+LpAK+LpAKL+LpP+LpAP+2*LpAPLp;
liptot=ltot+lptot;
tots = [ktot, ptot, atot, liptot]

#timespan for integration
tspan = (0., 100.);

#solve the reduced ODEs
prob = ODEProblem(reduced_oscillator_odes!, u0[1:2:3], tspan, vcat(p, tots))
sol = solve(prob) #solve adaptively
sol_fixed = solve(prob, Tsit5(), saveat=0.1) #solve with fixed time steps
sol_fixed2 = solve(prob, Tsit5(), saveat=0.01) #solve with smaller fixed time steps
sol_stiff = solve(prob, AutoTsit5(Rosenbrock23()), saveat=0.01) #solve with stiff solver
sol_interp = [sol(t) for t in range(0, stop=100, length=1000)] #interpolate solution to a finer grid

# output = [sol(t) for t in range(0, stop=100, length=1000)]
# sol = solve(remake(prob, p=vcat(p, tots)))


#solve the full model
prob2 = ODEProblem(osc_rn, umap, tspan, p)
sol2 = solve(prob2, save_idxs=[1,6])

#plot the results
p1 = plot(sol, label = ["L" "LpA"] ,  lw=2, title="Reduced Oscillator Model", xlabel="Time (s)", ylabel="Concentration");
p2 = plot(sol2, label = ["L" "LpA"], lw=2, title="Full Oscillator Model", xlabel="Time (s)", ylabel="Concentration");
plot(p1, p2, layout=(2,1), size=(800,800))





## Genetic algorithm to find the best parameters for the reduced model ## 
using Evolutionary
using FFTW

## COST FUNCTIONS and dependencies 
function getDif(indexes::Vector{Int}, arrayData::Vector{Float64}) #get summed difference between fft peak indexes
    arrLen = length(indexes)
    if arrLen < 2
        return 0.0
    end
    sum_diff = @inbounds sum(arrayData[indexes[i]] - arrayData[indexes[i+1]] for i in 1:(arrLen-1))
    sum_diff += arrayData[indexes[end]] #add the last element
    return sum_diff
end

function getSTD(peakindxs::Vector{Int}, arrayData::Vector{Float64}, window_ratio::Float64) #get average standard deviation of fft peak indexes
    window = max(1, round(Int, window_ratio * length(arrayData)))
    sum_std = @inbounds sum(std(arrayData[max(1, ind - window):min(length(arrayData), ind + window)]) for ind in peakindxs)
    return sum_std / length(peakindxs)
end


function getFrequencies(y::Vector{Float64})
    res = abs.(rfft(y))
    return res ./ cld(length(y), 2) #normalize amplitudes
end



function eval_fitness(p::Vector{Float64}, tots::Vector{Float64},  prob::ODEProblem, idxs)
    Ysol = nothing
    try 
        Ysol = solve(remake(prob, p=vcat(p,tots)), save_idxs=idxs)
    catch e
        if isa(e, DomainError)
            return 0.0
        else
            rethrow(e)
        end
    end

    if any(isnan.(Ysol.u)) || any(isless.(Ysol.u, 0.0))
        return 0.0
    end

    #interpolate evenly spaced time steps for the solution
    Y = [Ysol(t) for t in range(0, stop=100, length=1000)]

    #get the fft of the solution
    fftData = getFrequencies(Y)
    indexes = findmaxima(fftData,1)[1]
    if isempty(indexes)
        return 0.0
    end
    std = getSTD(indexes, fftData, 0.1)
    diff = getDif(indexes, fftData)
    return -std - diff
end


function make_fitness_function(prob::ODEProblem, tots; idxs = 1) # Create a fitness function that includes your ODE problem as a constant
    function fitness_function(p::Vector{Float64})
        return eval_fitness(p, tots, prob, idxs)  
    end
    return fitness_function
end

fitness_function = make_fitness_function(prob, tots) # Create a fitness function that includes your ODE problem as a constant



## parameter constraint ranges ##
ka_min, ka_max = 0.001, 1.
kb_min, kb_max = 0.001, 500.0
kcat_min, kcat_max = 0.001, 500.0


param_values = Dict(
    "ka1" => Dict("min" => ka_min, "max" => ka_max),
    "kb1" => Dict("min" => kb_min, "max" => kb_max),
    "kcat1" => Dict("min" => kcat_min, "max" => kcat_max),
    "ka2" => Dict("min" => ka_min, "max" => ka_max),
    "kb2" => Dict("min" => kb_min, "max" => kb_max),
    "ka3" => Dict("min" => ka_min, "max" => ka_max),
    "kb3" => Dict("min" => kb_min, "max" => kb_max),
    "ka4" => Dict("min" => ka_min, "max" => ka_max),
    "kb4" => Dict("min" => kb_min, "max" => kb_max),
    "ka7" => Dict("min" => ka_min, "max" => ka_max),
    "kb7" => Dict("min" => kb_min, "max" => kb_max),
    "kcat7" => Dict("min" => kcat_min, "max" => kcat_max),
    "y" => Dict("min" => 100., "max" => 5000.)
);

random_p = [rand(param_values[p]["min"]:0.01:param_values[p]["max"]) for p in keys(param_values)]
eval_fitness(random_p, tots, prob, 1)
# eval_fitness(random_p, tots, prob2, 1)
testprob = ODEProblem(reduced_oscillator_odes!, u0[1:2:3], tspan, vcat(random_p, tots))
testsol = solve(testprob,Tsit5(), p=vcat(p, tots), dense=true)
findmaxima(testsol[1,:], 10)
plot(testsol)

#Optimization parameters
constraints = BoxConstraints([param_values[p]["min"] for p in keys(param_values)], [param_values[p]["max"] for p in keys(param_values)])
opts = Evolutionary.Options(show_trace=true,show_every=1, store_trace=true, iterations=10, parallelization=:thread, abstol=1e-6, reltol=1e-6)

common_range = 0.5
valrange = fill(common_range, 13)
mthd = GA(populationSize = 5000, selection = tournament(500),
          crossover = TPX, crossoverRate = 0.5,
          mutation  = BGA(valrange, 2), mutationRate = 0.9)


# function callback(trace)

#Optimization
result = Evolutionary.optimize(fitness_function, constraints, mthd, opts)

newp = result.minimizer
newsol = solve(remake(prob, p=vcat(newp, tots)))
newsol2 = solve(remake(prob2, p = vcat(newp,tots)), save_idxs=[1,6])

#plot the results
p1 = plot(newsol, label = ["L" "LpA"] ,  lw=2, title="Reduced Oscillator Model", xlabel="Time (s)", ylabel="Concentration");
p2 = plot(newsol2, label = ["L" "LpA"], lw=2, title="Full Oscillator Model", xlabel="Time (s)", ylabel="Concentration");
plot(p1, p2, layout=(2,1), size=(800,800))


eval_fitness(newp, tots, prob, 1)
eval_fitness(newp, tots, prob2, 1)