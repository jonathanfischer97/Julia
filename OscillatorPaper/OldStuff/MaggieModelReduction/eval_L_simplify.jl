using Symbolics

@variables ka1 kb1 ka3 kb3 ka4 kb4 ka7 kb7 kcat1 kcat7 ktot ptot liptot L LpA DF



expr=(kb1*kb3*kb4*kb7^2 + kb3*kb4*kb7^2*kcat1 + 2*kb1*kb3*kb4*kb7*kcat7 + 
         2*kb3*kb4*kb7*kcat1*kcat7 + kb1*kb3*kb4*kcat7^2 + kb3*kb4*kcat1*kcat7^2 + 
         ka7*kb1*kb3*kb4*kb7*L + ka1*kb3*kb4*kb7^2*L + ka7*kb3*kb4*kb7*kcat1*L + 
         ka7*kb1*kb3*kb4*kcat7*L + 2*ka1*kb3*kb4*kb7*kcat7*L + 
         ka7*kb3*kb4*kcat1*kcat7*L + ka1*kb3*kb4*kcat7^2*L + 
         ka1*ka7*kb3*kb4*kb7*ktot*L + ka1*ka7*kb3*kb4*kcat7*ktot*L + 
         ka1*ka7*kb3*kb4*kb7*L^2 + ka1*ka7*kb3*kb4*kcat7*L^2 - 
         ka7*kb1*kb3*kb4*kb7*liptot - ka7*kb3*kb4*kb7*kcat1*liptot - 
         ka7*kb1*kb3*kb4*kcat7*liptot - ka7*kb3*kb4*kcat1*kcat7*liptot - 
         ka1*ka7*kb3*kb4*kb7*L*liptot - ka1*ka7*kb3*kb4*kcat7*L*liptot + 
         ka7*kb1*kb3*kb4*kb7*LpA + ka4*kb1*kb3*kb7^2*LpA + 2*ka3*kb1*kb4*kb7^2*LpA + 
         ka7*kb3*kb4*kb7*kcat1*LpA + ka4*kb3*kb7^2*kcat1*LpA + 
         2*ka3*kb4*kb7^2*kcat1*LpA + ka7*kb1*kb3*kb4*kcat7*LpA + 
         2*ka4*kb1*kb3*kb7*kcat7*LpA + 4*ka3*kb1*kb4*kb7*kcat7*LpA + 
         ka7*kb3*kb4*kcat1*kcat7*LpA + 2*ka4*kb3*kb7*kcat1*kcat7*LpA + 
         4*ka3*kb4*kb7*kcat1*kcat7*LpA + ka4*kb1*kb3*kcat7^2*LpA + 
         2*ka3*kb1*kb4*kcat7^2*LpA + ka4*kb3*kcat1*kcat7^2*LpA + 
         2*ka3*kb4*kcat1*kcat7^2*LpA + ka3*ka7*kb1*kb4*kb7*ktot*LpA + 
         ka3*ka7*kb4*kb7*kcat1*ktot*LpA + ka3*ka7*kb1*kb4*kcat7*ktot*LpA + 
         ka3*ka7*kb4*kcat1*kcat7*ktot*LpA + DF*ka4*ka7*kb1*kb3*kb7*L*LpA + 
         2*ka3*ka7*kb1*kb4*kb7*L*LpA + ka1*ka7*kb3*kb4*kb7*L*LpA + 
         ka1*ka4*kb3*kb7^2*L*LpA + DF*ka4*ka7*kb3*kb7*kcat1*L*LpA + 
         2*ka3*ka7*kb4*kb7*kcat1*L*LpA + DF*ka4*ka7*kb1*kb3*kcat7*L*LpA + 
         2*ka3*ka7*kb1*kb4*kcat7*L*LpA + ka1*ka7*kb3*kb4*kcat7*L*LpA + 
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
         DF*ka4*ka7*kb1*kb3*kb7*LpA^2 + 2*ka3*ka7*kb1*kb4*kb7*LpA^2 + 
         2*ka3*ka4*kb1*kb7^2*LpA^2 + DF*ka4*ka7*kb3*kb7*kcat1*LpA^2 + 
         2*ka3*ka7*kb4*kb7*kcat1*LpA^2 + 2*ka3*ka4*kb7^2*kcat1*LpA^2 + 
         DF*ka4*ka7*kb1*kb3*kcat7*LpA^2 + 2*ka3*ka7*kb1*kb4*kcat7*LpA^2 + 
         4*ka3*ka4*kb1*kb7*kcat7*LpA^2 + DF*ka4*ka7*kb3*kcat1*kcat7*LpA^2 + 
         2*ka3*ka7*kb4*kcat1*kcat7*LpA^2 + 4*ka3*ka4*kb7*kcat1*kcat7*LpA^2 + 
         2*ka3*ka4*kb1*kcat7^2*LpA^2 + 2*ka3*ka4*kcat1*kcat7^2*LpA^2 + 
         DF*ka3*ka4*ka7*kb1*kb7*ktot*LpA^2 + DF*ka3*ka4*ka7*kb7*kcat1*ktot*LpA^2 + 
         DF*ka3*ka4*ka7*kb1*kcat7*ktot*LpA^2 + 
         DF*ka3*ka4*ka7*kcat1*kcat7*ktot*LpA^2 + 2*DF*ka3*ka4*ka7*kb1*kb7*L*LpA^2 + 
         DF*ka1*ka4*ka7*kb3*kb7*L*LpA^2 + 2*DF*ka3*ka4*ka7*kb7*kcat1*L*LpA^2 + 
         2*DF*ka3*ka4*ka7*kb1*kcat7*L*LpA^2 + DF*ka1*ka4*ka7*kb3*kcat7*L*LpA^2 + 
         2*DF*ka3*ka4*ka7*kcat1*kcat7*L*LpA^2 + 
         2*DF^2*ka1*ka3*ka4*ka7*kb7*ktot*L*LpA^2 + 
         2*DF^2*ka1*ka3*ka4*ka7*kcat7*ktot*L*LpA^2 - 
         2*DF*ka3*ka4*ka7*kb1*kb7*liptot*LpA^2 - 
         2*DF*ka3*ka4*ka7*kb7*kcat1*liptot*LpA^2 - 
         2*DF*ka3*ka4*ka7*kb1*kcat7*liptot*LpA^2 - 
         2*DF*ka3*ka4*ka7*kcat1*kcat7*liptot*LpA^2 + 2*DF*ka3*ka4*ka7*kb1*kb7*LpA^3 + 
         2*DF*ka3*ka4*ka7*kb7*kcat1*LpA^3 + 2*DF*ka3*ka4*ka7*kb1*kcat7*LpA^3 + 
         2*DF*ka3*ka4*ka7*kcat1*kcat7*LpA^3 + ka7*kb1*kb3*kb4*kb7*ptot + 
         ka7*kb3*kb4*kb7*kcat1*ptot + ka7*kb1*kb3*kb4*kcat7*ptot + 
         ka7*kb3*kb4*kcat1*kcat7*ptot + ka1*ka7*kb3*kb4*kb7*L*ptot + 
         ka1*ka7*kb3*kb4*kcat7*L*ptot + 2*DF*ka4*ka7*kb1*kb3*kb7*LpA*ptot + 
         2*ka3*ka7*kb1*kb4*kb7*LpA*ptot + 2*DF*ka4*ka7*kb3*kb7*kcat1*LpA*ptot + 
         2*ka3*ka7*kb4*kb7*kcat1*LpA*ptot + 2*DF*ka4*ka7*kb1*kb3*kcat7*LpA*ptot + 
         2*ka3*ka7*kb1*kb4*kcat7*LpA*ptot + 2*DF*ka4*ka7*kb3*kcat1*kcat7*LpA*ptot + 
         2*ka3*ka7*kb4*kcat1*kcat7*LpA*ptot + 2*DF*ka1*ka4*ka7*kb3*kb7*L*LpA*ptot + 
         2*DF*ka1*ka4*ka7*kb3*kcat7*L*LpA*ptot + 4*DF*ka3*ka4*ka7*kb1*kb7*LpA^2*ptot + 
         4*DF*ka3*ka4*ka7*kb7*kcat1*LpA^2*ptot + 
         4*DF*ka3*ka4*ka7*kb1*kcat7*LpA^2*ptot + 
         4*DF*ka3*ka4*ka7*kcat1*kcat7*LpA^2*ptot - 
         sqrt((-(kb1*kb3*kb4*kb7^2) - kb3*kb4*kb7^2*kcat1 - 2*kb1*kb3*kb4*kb7*kcat7 - 
              2*kb3*kb4*kb7*kcat1*kcat7 - kb1*kb3*kb4*kcat7^2 - 
              kb3*kb4*kcat1*kcat7^2 - ka7*kb1*kb3*kb4*kb7*L - ka1*kb3*kb4*kb7^2*L - 
              ka7*kb3*kb4*kb7*kcat1*L - ka7*kb1*kb3*kb4*kcat7*L - 
              2*ka1*kb3*kb4*kb7*kcat7*L - ka7*kb3*kb4*kcat1*kcat7*L - 
              ka1*kb3*kb4*kcat7^2*L - ka1*ka7*kb3*kb4*kb7*ktot*L - 
              ka1*ka7*kb3*kb4*kcat7*ktot*L - ka1*ka7*kb3*kb4*kb7*L^2 - 
              ka1*ka7*kb3*kb4*kcat7*L^2 + ka7*kb1*kb3*kb4*kb7*liptot + 
              ka7*kb3*kb4*kb7*kcat1*liptot + ka7*kb1*kb3*kb4*kcat7*liptot + 
              ka7*kb3*kb4*kcat1*kcat7*liptot + ka1*ka7*kb3*kb4*kb7*L*liptot + 
              ka1*ka7*kb3*kb4*kcat7*L*liptot - ka7*kb1*kb3*kb4*kb7*LpA - 
              ka4*kb1*kb3*kb7^2*LpA - 2*ka3*kb1*kb4*kb7^2*LpA - 
              ka7*kb3*kb4*kb7*kcat1*LpA - ka4*kb3*kb7^2*kcat1*LpA - 
              2*ka3*kb4*kb7^2*kcat1*LpA - ka7*kb1*kb3*kb4*kcat7*LpA - 
              2*ka4*kb1*kb3*kb7*kcat7*LpA - 4*ka3*kb1*kb4*kb7*kcat7*LpA - 
              ka7*kb3*kb4*kcat1*kcat7*LpA - 2*ka4*kb3*kb7*kcat1*kcat7*LpA - 
              4*ka3*kb4*kb7*kcat1*kcat7*LpA - ka4*kb1*kb3*kcat7^2*LpA - 
              2*ka3*kb1*kb4*kcat7^2*LpA - ka4*kb3*kcat1*kcat7^2*LpA - 
              2*ka3*kb4*kcat1*kcat7^2*LpA - ka3*ka7*kb1*kb4*kb7*ktot*LpA - 
              ka3*ka7*kb4*kb7*kcat1*ktot*LpA - ka3*ka7*kb1*kb4*kcat7*ktot*LpA - 
              ka3*ka7*kb4*kcat1*kcat7*ktot*LpA - DF*ka4*ka7*kb1*kb3*kb7*L*LpA - 
              2*ka3*ka7*kb1*kb4*kb7*L*LpA - ka1*ka7*kb3*kb4*kb7*L*LpA - 
              ka1*ka4*kb3*kb7^2*L*LpA - DF*ka4*ka7*kb3*kb7*kcat1*L*LpA - 
              2*ka3*ka7*kb4*kb7*kcat1*L*LpA - DF*ka4*ka7*kb1*kb3*kcat7*L*LpA - 
              2*ka3*ka7*kb1*kb4*kcat7*L*LpA - ka1*ka7*kb3*kb4*kcat7*L*LpA - 
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
              DF*ka1*ka4*ka7*kb3*kcat7*L*liptot*LpA - DF*ka4*ka7*kb1*kb3*kb7*LpA^2 - 
              2*ka3*ka7*kb1*kb4*kb7*LpA^2 - 2*ka3*ka4*kb1*kb7^2*LpA^2 - 
              DF*ka4*ka7*kb3*kb7*kcat1*LpA^2 - 2*ka3*ka7*kb4*kb7*kcat1*LpA^2 - 
              2*ka3*ka4*kb7^2*kcat1*LpA^2 - DF*ka4*ka7*kb1*kb3*kcat7*LpA^2 - 
              2*ka3*ka7*kb1*kb4*kcat7*LpA^2 - 4*ka3*ka4*kb1*kb7*kcat7*LpA^2 - 
              DF*ka4*ka7*kb3*kcat1*kcat7*LpA^2 - 2*ka3*ka7*kb4*kcat1*kcat7*LpA^2 - 
              4*ka3*ka4*kb7*kcat1*kcat7*LpA^2 - 2*ka3*ka4*kb1*kcat7^2*LpA^2 - 
              2*ka3*ka4*kcat1*kcat7^2*LpA^2 - DF*ka3*ka4*ka7*kb1*kb7*ktot*LpA^2 - 
              DF*ka3*ka4*ka7*kb7*kcat1*ktot*LpA^2 - 
              DF*ka3*ka4*ka7*kb1*kcat7*ktot*LpA^2 - 
              DF*ka3*ka4*ka7*kcat1*kcat7*ktot*LpA^2 - 
              2*DF*ka3*ka4*ka7*kb1*kb7*L*LpA^2 - DF*ka1*ka4*ka7*kb3*kb7*L*LpA^2 - 
              2*DF*ka3*ka4*ka7*kb7*kcat1*L*LpA^2 - 
              2*DF*ka3*ka4*ka7*kb1*kcat7*L*LpA^2 - DF*ka1*ka4*ka7*kb3*kcat7*L*LpA^2 - 
              2*DF*ka3*ka4*ka7*kcat1*kcat7*L*LpA^2 - 
              2*DF^2*ka1*ka3*ka4*ka7*kb7*ktot*L*LpA^2 - 
              2*DF^2*ka1*ka3*ka4*ka7*kcat7*ktot*L*LpA^2 + 
              2*DF*ka3*ka4*ka7*kb1*kb7*liptot*LpA^2 + 
              2*DF*ka3*ka4*ka7*kb7*kcat1*liptot*LpA^2 + 
              2*DF*ka3*ka4*ka7*kb1*kcat7*liptot*LpA^2 + 
              2*DF*ka3*ka4*ka7*kcat1*kcat7*liptot*LpA^2 - 
              2*DF*ka3*ka4*ka7*kb1*kb7*LpA^3 - 2*DF*ka3*ka4*ka7*kb7*kcat1*LpA^3 - 
              2*DF*ka3*ka4*ka7*kb1*kcat7*LpA^3 - 2*DF*ka3*ka4*ka7*kcat1*kcat7*LpA^3 - 
              ka7*kb1*kb3*kb4*kb7*ptot - ka7*kb3*kb4*kb7*kcat1*ptot - 
              ka7*kb1*kb3*kb4*kcat7*ptot - ka7*kb3*kb4*kcat1*kcat7*ptot - 
              ka1*ka7*kb3*kb4*kb7*L*ptot - ka1*ka7*kb3*kb4*kcat7*L*ptot - 
              2*DF*ka4*ka7*kb1*kb3*kb7*LpA*ptot - 2*ka3*ka7*kb1*kb4*kb7*LpA*ptot - 
              2*DF*ka4*ka7*kb3*kb7*kcat1*LpA*ptot - 2*ka3*ka7*kb4*kb7*kcat1*LpA*ptot - 
              2*DF*ka4*ka7*kb1*kb3*kcat7*LpA*ptot - 2*ka3*ka7*kb1*kb4*kcat7*LpA*ptot - 
              2*DF*ka4*ka7*kb3*kcat1*kcat7*LpA*ptot - 
              2*ka3*ka7*kb4*kcat1*kcat7*LpA*ptot - 
              2*DF*ka1*ka4*ka7*kb3*kb7*L*LpA*ptot - 
              2*DF*ka1*ka4*ka7*kb3*kcat7*L*LpA*ptot - 
              4*DF*ka3*ka4*ka7*kb1*kb7*LpA^2*ptot - 
              4*DF*ka3*ka4*ka7*kb7*kcat1*LpA^2*ptot - 
              4*DF*ka3*ka4*ka7*kb1*kcat7*LpA^2*ptot - 
              4*DF*ka3*ka4*ka7*kcat1*kcat7*LpA^2*ptot)^2 - 
           4*(-(ka7*kb1*kb3*kb4*kb7) - ka7*kb3*kb4*kb7*kcat1 - ka7*kb1*kb3*kb4*kcat7 - 
              ka7*kb3*kb4*kcat1*kcat7 - ka1*ka7*kb3*kb4*kb7*L - 
              ka1*ka7*kb3*kb4*kcat7*L - DF*ka4*ka7*kb1*kb3*kb7*LpA - 
              2*ka3*ka7*kb1*kb4*kb7*LpA - DF*ka4*ka7*kb3*kb7*kcat1*LpA - 
              2*ka3*ka7*kb4*kb7*kcat1*LpA - DF*ka4*ka7*kb1*kb3*kcat7*LpA - 
              2*ka3*ka7*kb1*kb4*kcat7*LpA - DF*ka4*ka7*kb3*kcat1*kcat7*LpA - 
              2*ka3*ka7*kb4*kcat1*kcat7*LpA - DF*ka1*ka4*ka7*kb3*kb7*L*LpA - 
              DF*ka1*ka4*ka7*kb3*kcat7*L*LpA - 2*DF*ka3*ka4*ka7*kb1*kb7*LpA^2 - 
              2*DF*ka3*ka4*ka7*kb7*kcat1*LpA^2 - 2*DF*ka3*ka4*ka7*kb1*kcat7*LpA^2 - 
              2*DF*ka3*ka4*ka7*kcat1*kcat7*LpA^2)*
            (-(kb1*kb3*kb4*kb7^2*L) - kb3*kb4*kb7^2*kcat1*L - 
              2*kb1*kb3*kb4*kb7*kcat7*L - 2*kb3*kb4*kb7*kcat1*kcat7*L - 
              kb1*kb3*kb4*kcat7^2*L - kb3*kb4*kcat1*kcat7^2*L - 
              ka1*kb3*kb4*kb7^2*ktot*L - 2*ka1*kb3*kb4*kb7*kcat7*ktot*L - 
              ka1*kb3*kb4*kcat7^2*ktot*L - ka1*kb3*kb4*kb7^2*L^2 - 
              2*ka1*kb3*kb4*kb7*kcat7*L^2 - ka1*kb3*kb4*kcat7^2*L^2 + 
              kb1*kb3*kb4*kb7^2*liptot + kb3*kb4*kb7^2*kcat1*liptot + 
              2*kb1*kb3*kb4*kb7*kcat7*liptot + 2*kb3*kb4*kb7*kcat1*kcat7*liptot + 
              kb1*kb3*kb4*kcat7^2*liptot + kb3*kb4*kcat1*kcat7^2*liptot + 
              ka1*kb3*kb4*kb7^2*L*liptot + 2*ka1*kb3*kb4*kb7*kcat7*L*liptot + 
              ka1*kb3*kb4*kcat7^2*L*liptot - kb1*kb3*kb4*kb7^2*LpA - 
              kb3*kb4*kb7^2*kcat1*LpA - 2*kb1*kb3*kb4*kb7*kcat7*LpA - 
              2*kb3*kb4*kb7*kcat1*kcat7*LpA - kb1*kb3*kb4*kcat7^2*LpA - 
              kb3*kb4*kcat1*kcat7^2*LpA - ka3*kb1*kb4*kb7^2*ktot*LpA - 
              ka3*kb4*kb7^2*kcat1*ktot*LpA - 2*ka3*kb1*kb4*kb7*kcat7*ktot*LpA - 
              2*ka3*kb4*kb7*kcat1*kcat7*ktot*LpA - ka3*kb1*kb4*kcat7^2*ktot*LpA - 
              ka3*kb4*kcat1*kcat7^2*ktot*LpA - ka4*kb1*kb3*kb7^2*L*LpA - 
              2*ka3*kb1*kb4*kb7^2*L*LpA - ka1*kb3*kb4*kb7^2*L*LpA - 
              ka4*kb3*kb7^2*kcat1*L*LpA - 2*ka3*kb4*kb7^2*kcat1*L*LpA - 
              2*ka4*kb1*kb3*kb7*kcat7*L*LpA - 4*ka3*kb1*kb4*kb7*kcat7*L*LpA - 
              2*ka1*kb3*kb4*kb7*kcat7*L*LpA - 2*ka4*kb3*kb7*kcat1*kcat7*L*LpA - 
              4*ka3*kb4*kb7*kcat1*kcat7*L*LpA - ka4*kb1*kb3*kcat7^2*L*LpA - 
              2*ka3*kb1*kb4*kcat7^2*L*LpA - ka1*kb3*kb4*kcat7^2*L*LpA - 
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
              ka1*ka4*kb3*kcat7^2*L*liptot*LpA - ka4*kb1*kb3*kb7^2*LpA^2 - 
              2*ka3*kb1*kb4*kb7^2*LpA^2 - ka4*kb3*kb7^2*kcat1*LpA^2 - 
              2*ka3*kb4*kb7^2*kcat1*LpA^2 - 2*ka4*kb1*kb3*kb7*kcat7*LpA^2 - 
              4*ka3*kb1*kb4*kb7*kcat7*LpA^2 - 2*ka4*kb3*kb7*kcat1*kcat7*LpA^2 - 
              4*ka3*kb4*kb7*kcat1*kcat7*LpA^2 - ka4*kb1*kb3*kcat7^2*LpA^2 - 
              2*ka3*kb1*kb4*kcat7^2*LpA^2 - ka4*kb3*kcat1*kcat7^2*LpA^2 - 
              2*ka3*kb4*kcat1*kcat7^2*LpA^2 - ka3*ka4*kb1*kb7^2*ktot*LpA^2 - 
              ka3*ka4*kb7^2*kcat1*ktot*LpA^2 - 2*ka3*ka4*kb1*kb7*kcat7*ktot*LpA^2 - 
              2*ka3*ka4*kb7*kcat1*kcat7*ktot*LpA^2 - 
              ka3*ka4*kb1*kcat7^2*ktot*LpA^2 - ka3*ka4*kcat1*kcat7^2*ktot*LpA^2 - 
              2*ka3*ka4*kb1*kb7^2*L*LpA^2 - ka1*ka4*kb3*kb7^2*L*LpA^2 - 
              2*ka3*ka4*kb7^2*kcat1*L*LpA^2 - 4*ka3*ka4*kb1*kb7*kcat7*L*LpA^2 - 
              2*ka1*ka4*kb3*kb7*kcat7*L*LpA^2 - 4*ka3*ka4*kb7*kcat1*kcat7*L*LpA^2 - 
              2*ka3*ka4*kb1*kcat7^2*L*LpA^2 - ka1*ka4*kb3*kcat7^2*L*LpA^2 - 
              2*ka3*ka4*kcat1*kcat7^2*L*LpA^2 - 
              2*DF*ka1*ka3*ka4*kb7^2*ktot*L*LpA^2 - 
              4*DF*ka1*ka3*ka4*kb7*kcat7*ktot*L*LpA^2 - 
              2*DF*ka1*ka3*ka4*kcat7^2*ktot*L*LpA^2 + 
              2*ka3*ka4*kb1*kb7^2*liptot*LpA^2 + 
              2*ka3*ka4*kb7^2*kcat1*liptot*LpA^2 + 
              4*ka3*ka4*kb1*kb7*kcat7*liptot*LpA^2 + 
              4*ka3*ka4*kb7*kcat1*kcat7*liptot*LpA^2 + 
              2*ka3*ka4*kb1*kcat7^2*liptot*LpA^2 + 
              2*ka3*ka4*kcat1*kcat7^2*liptot*LpA^2 - 2*ka3*ka4*kb1*kb7^2*LpA^3 - 
              2*ka3*ka4*kb7^2*kcat1*LpA^3 - 4*ka3*ka4*kb1*kb7*kcat7*LpA^3 - 
              4*ka3*ka4*kb7*kcat1*kcat7*LpA^3 - 2*ka3*ka4*kb1*kcat7^2*LpA^3 - 
              2*ka3*ka4*kcat1*kcat7^2*LpA^3 - ka4*kb1*kb3*kb7^2*LpA*ptot - 
              ka4*kb3*kb7^2*kcat1*LpA*ptot - 2*ka4*kb1*kb3*kb7*kcat7*LpA*ptot - 
              2*ka4*kb3*kb7*kcat1*kcat7*LpA*ptot - ka4*kb1*kb3*kcat7^2*LpA*ptot - 
              ka4*kb3*kcat1*kcat7^2*LpA*ptot - ka1*ka4*kb3*kb7^2*L*LpA*ptot - 
              2*ka1*ka4*kb3*kb7*kcat7*L*LpA*ptot - ka1*ka4*kb3*kcat7^2*L*LpA*ptot - 
              2*ka3*ka4*kb1*kb7^2*LpA^2*ptot - 2*ka3*ka4*kb7^2*kcat1*LpA^2*ptot - 
              4*ka3*ka4*kb1*kb7*kcat7*LpA^2*ptot - 
              4*ka3*ka4*kb7*kcat1*kcat7*LpA^2*ptot - 
              2*ka3*ka4*kb1*kcat7^2*LpA^2*ptot - 2*ka3*ka4*kcat1*kcat7^2*LpA^2*ptot)
           ))/(2*(-(ka7*kb1*kb3*kb4*kb7) - ka7*kb3*kb4*kb7*kcat1 - 
           ka7*kb1*kb3*kb4*kcat7 - ka7*kb3*kb4*kcat1*kcat7 - ka1*ka7*kb3*kb4*kb7*L - 
           ka1*ka7*kb3*kb4*kcat7*L - DF*ka4*ka7*kb1*kb3*kb7*LpA - 
           2*ka3*ka7*kb1*kb4*kb7*LpA - DF*ka4*ka7*kb3*kb7*kcat1*LpA - 
           2*ka3*ka7*kb4*kb7*kcat1*LpA - DF*ka4*ka7*kb1*kb3*kcat7*LpA - 
           2*ka3*ka7*kb1*kb4*kcat7*LpA - DF*ka4*ka7*kb3*kcat1*kcat7*LpA - 
           2*ka3*ka7*kb4*kcat1*kcat7*LpA - DF*ka1*ka4*ka7*kb3*kb7*L*LpA - 
           DF*ka1*ka4*ka7*kb3*kcat7*L*LpA - 2*DF*ka3*ka4*ka7*kb1*kb7*LpA^2 - 
           2*DF*ka3*ka4*ka7*kb7*kcat1*LpA^2 - 2*DF*ka3*ka4*ka7*kb1*kcat7*LpA^2 - 
				   2*DF*ka3*ka4*ka7*kcat1*kcat7*LpA^2))


simplified_expr = simplify(expr, threaded=true)
println(simplified_expr)

