using ModelingToolkit, Symbolics

# Define the variables
# @variables L Lp K P LpA LK LpP LpAK LpAP LpAKL LpAPLp
diffvars = @variables t L(t) LpA(t) 
algvars = @variables Lp(t) K(t) P(t) A(t) LK(t) LpP(t) LpAK(t) LpAP(t) LpAKL(t) LpAPLp(t)
parms = @parameters ka1 kb1 kcat1 ka2 kb2 ka3 kb3 ka4 kb4 ka7 kb7 kcat7 DF 
masstots = @parameters ktot ptot atot liptot 

function getLp(L, LpA, p, tots)
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



@constants begin common_sum = kb1 * kb3 * (L - liptot + Lp + LpA - LpAP + 2 * ptot) +
        kb3 * kcat1 * (L - liptot + Lp + LpA - LpAP + 2 * ptot) +
        ka1 * kb3 * L * (ktot + L - liptot + Lp + LpA - LpAP + 2 * ptot) +
        ka3 * kb1 * LpA * (ktot + L - liptot + Lp + LpA - LpAP + 2 * ptot) +
        ka3 * kcat1 * LpA * (ktot + L - liptot + Lp + LpA - LpAP + 2 * ptot) +
        DF * ka1 * ka3 * L * LpA * (2 * ktot + L - liptot + Lp + LpA - LpAP + 2 * ptot)

        denominator = (2 * kb7 + 2 * kcat7 + ka7 * Lp) * (kb3 * (kcat1 + ka1 * L) + ka3 * (kcat1 + DF * ka1 * L) * LpA + kb1 * (kb3 + ka3 * LpA))

        LK_denom = (kb1 + kcat1)^2 + 2 * DF * ka1 * (kb1 + kcat1) * L + DF * ka1^2 * L^2
end

# Define the algebraic equations
algeqs = [
        Lp ~ getLp(L, LpA, parms, masstots),

        LpAP ~ -(((kb7 + kcat7)*(kb3*(kb7 + kcat7 + ka7*Lp - ka4*LpA)*
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
                ka4*(kb7 + kcat7)*LpA))),

        P ~ ((kb7 + kcat7) * common_sum) / denominator,

        LpP ~ ka7 * Lp * common_sum / denominator,

        LK ~ ka1 * L * (kb1 * (ktot + L - liptot + Lp + LpA - LpAP - LpP - 2 * P + 2 * ptot) + kcat1 * (ktot + L - liptot + Lp + LpA - LpAP - LpP - 2 * P + 2 * ptot) + DF * ka1 * L * 
                (2 * ktot + L - liptot + Lp + LpA - LpAP - LpP - 2 * P + 2 * ptot)) / LK_denom,
    
        LpAK ~ -((kb1 + kcat1) * (kb1 * (L - liptot + Lp + LpA - LpAP - LpP - 2 * P + 2 * ptot) + kcat1 * (L - liptot + Lp + LpA - LpAP - LpP - 2 * P + 2 * ptot) + ka1 * L * 
                (ktot + L - liptot + Lp + LpA - LpAP - LpP - 2 * P + 2 * ptot))) / LK_denom,
    
    
        #mass constraints 
        LpAPLp~ -LpAP - LpP - P + ptot,

        LpAKL ~ 0.5*(-L + liptot - LK - Lp - LpA - LpAK + LpAP + LpP + 2*P - 2*ptot),

        A ~ 0.5*(2*atot + L - liptot + LK + Lp - LpA - LpAK - LpAP + LpP),

        K ~ 0.5* (2*ktot + L - liptot - LK + Lp + LpA - LpAK - LpAP - LpP - 2*P + 2*ptot),
]

D = Differential(t)
# Define the differential equations
diffeqs = [D(L) ~ kb1*LK + kb1*LpAKL + kcat7*LpAPLp + kcat7*LpP - ka1*K*L - ka1*DF*L*LpAK,
            D(LpA) ~ kb3*LpAK + kb4*LpAP + ka2*A*Lp - kb2*LpA - ka3*K*LpA - ka4*LpA*P]

# Define the ODE system
@named reduced_oscillator_odes = ODESystem(vcat(algeqs,diffeqs), t)

simplified_oscillator_odes = structural_simplify(reduced_oscillator_odes)