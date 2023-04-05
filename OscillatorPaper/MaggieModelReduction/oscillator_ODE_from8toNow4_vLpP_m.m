%The oscillator model has 12 variables, or 12 unknowns, and 12 differential
%equations plus 4 mass conservation equations. 

%THIS ONE WORKS WITH 4 VARIABLE UNKNOWNS:
%L, Lp, LpA, and LpAP.

%So we eliminated A, LpAKL, LpAPLp and K using mass conservation--this is
%exact.
%We eliminated LK and LpAK and P and LpP by 4 steady-state approxes. 
%%set dLK=0 and dLpAKL=0, solve for LK and LpAK. 
%set dLpP=0 and dK=0, solve for P and LpP

%L+K is reaction 1: ka1, kb1, kcat1
%Lp+A is reaction 2: ka2, kb2
%A+K is reaction 3: ka3, kb3
%A+P is reaction 4: ka4, kb4
%Lp+P is reaction 7: ka7, kb7, kcat7

%INPUT: 
%1) all rates constants in parmVec
%2) All initial concentrations in ICvec
%3) Time window or timePoint array
%Units of the forwards rates must match units of ICvec, i.e. k: uM-1s-1,
%and ICs are in uM. Or if IC are in copy numbers, the forward rates have to
%have the volume divided out, to be in /copy/s. 
%Time is in seconds
%
%Output
%1) timepoints of solution
%2) structure containing all concentrations/copies of species. 
%SAMPLE INPUTS:
% parm = 
% 
%   struct with fields:
% 
%       ka1: 0.055000000000000
%       kb1: 19.800000000000001
%     kcat1: 241
%       ka2: 1
%       kb2: 0.950000000000000
%       ka3: 41
%       kb3: 193
%       kb4: 0.130000000000000
%       ka7: 0.620000000000000
%       kb7: 3.390000000000000
%     kcat7: 4.600000000000000
%        DF: 750
%       ka4: 0.190000000000000
% ic = 
% 
%   struct with fields:
% 
%          L: 0.200000000000000
%         Lp: 3.000000000000000
%          K: 0.200000000000000
%          P: 0.300000000000000
%          A: 0.600000000000000
%         LK: 0.100000000000000
%        LpA: 0.010000000000000
%       LpAK: 0.100000000000000
%     LpAPLp: 0
%      LpAKL: 0
%        LpP: 0.100000000000000
%       LpAP: 0.010000000000000
% time=[0,100]
function[timePts, cs]=oscillator_ODE_from8toNow4_vLpP(parmVec, ICvec, time)

%initialize copy numbers for all 12 species. All units should be matched
%between species and on-rates
L=ICvec.L;
Lp=ICvec.Lp;
K=ICvec.K;
LK=ICvec.LK;
A=ICvec.A;
LpA=ICvec.LpA;
LpAK=ICvec.LpAK;
LpAKL=ICvec.LpAKL;
P=ICvec.P;
LpP=ICvec.LpP;
LpAP=ICvec.LpAP;
LpAPLp=ICvec.LpAPLp;

%INITIALIZE TOTAL MASS
ktot=K+LK+LpAK+LpAKL
ptot=P+LpP+LpAP+LpAPLp
atot=A+LpA+LpAK+LpAKL+LpAP+LpAPLp
ltot=L+LK+LpAKL
lptot=Lp+LpA+LpAK+LpAKL+LpP+LpAP+2*LpAPLp
liptot=ltot+lptot;
mass.ktot=ktot;
mass.ptot=ptot;
mass.atot=atot;
mass.ltot=ltot;
mass.lptot=lptot;
mass.liptot=liptot;
%initialize parameters
ka1=parmVec.ka1;
kb1=parmVec.kb1;
kcat1=parmVec.kcat1;
ka2=parmVec.ka2;
kb2=parmVec.kb2;
ka3=parmVec.ka3;
kb3=parmVec.kb3;
ka4=parmVec.ka4;
kb4=parmVec.kb4;
ka7=parmVec.ka7;
kb7=parmVec.kb7;
kcat7=parmVec.kcat7;
DF=parmVec.DF;%this is the DF

%initial conditions of the unknowns
    
 y0(1)=L;
y0(2)=Lp;
y0(3)=LpA;
y0(4)=LpAP;


 display('initial VALUES')
 [cs.A(1), cs.K(1), cs.LpAKL(1), cs.LpAPLp(1),cs.LpAK(1), cs.LK(1), cs.P(1), cs.LpP(1)]=calc_other_vars(y0,parmVec, mass);
  display('LpAK')
  cs.LpAK(1)
   display('LpP')
   cs.LpP(1)
   display('A')
   cs.A(1)
   display('K')
   cs.K(1)
   display('P')
   cs.P(1)
   display('LpAKL')
   cs.LpAKL(1)
   display('LpAPLp')
   cs.LpAPLp(1)
 
 %display begin iterating!
 display('START ITERATING THE ODE')
 
%%%%%%%%%%%%SOLVE THE ODE, USE A STIFF SOLVER %%%%%%%%%%%%%%%% 
opt=odeset('RelTol',1E-7,'AbsTol',1E-9);
 [timePts,conc] = ode15s(@(t,y) odes_4(t,y,parmVec, mass),time,y0,opt);
  
%map back to vars
cs.time=timePts;    
cs.L=conc(:,1);

cs.Lp=conc(:,2);

cs.LpA=conc(:,3);
cs.LpAP=conc(:,4);
%cs.LpAPLp=conc(:,12);
N=length(conc);
cs.A=zeros(N,1);
cs.K=zeros(N,1);
cs.LpAKL=zeros(N,1);
cs.LpAPLp=zeros(N,1);
cs.LK=zeros(N,1);
cs.LpAK=zeros(N,1);
cs.P=zeros(N,1);
cs.LpP=zeros(N,1);
for i=1:1:N   
    [cs.A(i), cs.K(i), cs.LpAKL(i), cs.LpAPLp(i),cs.LpAK(i), cs.LK(i), cs.P(i), cs.LpP(i)]=calc_other_vars(conc(i,:),parmVec, mass);
end
 %PLOT SECONDARY VARIABLES
  fignum=30
 try close(fignum+2)
 end
 f=figure(fignum+2)
 ax=axes('Parent',f,'FontSize',20,'LineWidth',1,'XScale','linear','YScale','linear');
 hold(ax)
 plot(timePts, cs.LpAK,'c--','LineWidth',2)%LpAK
 plot(timePts, cs.LpP,'b--','LineWidth',2)%LpP
 plot(timePts, cs.LK,'k--','LineWidth',2)%LK
 plot(timePts, cs.LpAP,'m--','LineWidth',2)%LpA
  plot(timePts, cs.LpAKL,'g--','LineWidth',2)%LpAKL
   plot(timePts, cs.LpAPLp,'y--','LineWidth',2)%LpAPLp
   plot(timePts, cs.K,'k-','LineWidth',2)%LpAKL
   plot(timePts, cs.P,'b-','LineWidth',2)%LpAPLp
 
 legend('LpAK','LpP','LK','LpAP','LpAKL','LpAPLp','K','P')
 
 %CALCULATE THE TOTAL MASS VS TIME
fktot=cs.K+cs.LK+cs.LpAK+cs.LpAKL;
fptot=cs.P+cs.LpP+cs.LpAP+cs.LpAPLp;
fatot=cs.A+cs.LpA+cs.LpAK+cs.LpAKL+cs.LpAP+cs.LpAPLp;
fltot=cs.L+cs.LK+cs.LpAKL;
flptot=cs.Lp+cs.LpA+cs.LpAK+cs.LpAKL+cs.LpP+cs.LpAP+2*cs.LpAPLp;
lipidTot=fltot+flptot;
try
     close(fignum+1)
end
 %PLOT TOTAL MASS
f=figure(fignum+1);
 ax=axes('Parent',f,'FontSize',20,'LineWidth',1,'XScale','linear','YScale','linear');
 hold(ax)
 plot(timePts, fktot,'r-','LineWidth',2)%K
 plot(timePts, fptot,'b-','LineWidth',2)%P
 plot(timePts, fatot,'k--','LineWidth',2)%A
 plot(timePts, fltot,'m--','LineWidth',2)%L
 plot(timePts, flptot,'c--','LineWidth',2)%Lp
 plot(timePts, lipidTot,'g--','LineWidth',2)%Lp
 legend('Ktot','Ptot','Atot','Ltot','Lptot','LipidTotal')
 
 
 
  %PLOT THE TIME VARIANCE OF THE UNKNOWNS
 try
     close(fignum)
 end
 
 f=figure(fignum)
 ax=axes('Parent',f,'FontSize',20,'LineWidth',1,'XScale','linear','YScale','linear');
 hold(ax)
 plot(timePts, cs.L,'r-','LineWidth',2)%L
 plot(timePts, cs.Lp,'b-','LineWidth',2)%Lp
 plot(timePts, cs.A,'k--','LineWidth',2)%A
 plot(timePts, cs.LpA,'m--','LineWidth',2)%LpA
 plot(timePts, cs.LpAK,'c--','LineWidth',2)%LpAK
 legend('L','Lp','A','LpA','LpAK')
 
end %END OF MAIN


%ODE DEFINITIONS
function dy = odes_4(t,y,parmVec, mass)

dy = zeros(2,1);    % a column vector

%initialize parameters
ka1=parmVec.ka1;
kb1=parmVec.kb1;
kcat1=parmVec.kcat1;
ka2=parmVec.ka2;
kb2=parmVec.kb2;
ka3=parmVec.ka3;
kb3=parmVec.kb3;
ka4=parmVec.ka4;
kb4=parmVec.kb4;
ka7=parmVec.ka7;
kb7=parmVec.kb7;
kcat7=parmVec.kcat7;
DF=parmVec.DF;%this is the DF
%mass conservation
ktot=mass.ktot;
ptot=mass.ptot;
atot=mass.atot;
liptot=mass.liptot;

%DEFINE CONSTANTS:


%L y1
%Lp y2
%LpA y3
%LpAP y4
   
%set dLpP=0 and dK=0, solve for P and LpP
P =((kb7 + kcat7)* (kb1 *kb3* (y(1) - liptot + y(2) + y(3) - y(4) + 2* ptot) +...
       kb3* kcat1* (y(1) - liptot + y(2) + y(3) - y(4) + 2 *ptot) + ...
      ka1 *kb3* y(1)* (ktot + y(1) - liptot + y(2) + y(3) - y(4) + 2* ptot) + ...
      ka3 *kb1* y(3)* (ktot + y(1) - liptot + y(2) + y(3) - y(4) + 2* ptot) + ...
      ka3* kcat1 *y(3) *(ktot + y(1) - liptot + y(2) + y(3) - y(4) + 2* ptot) + ...
      DF* ka1* ka3 *y(1)* y(3)* (2 *ktot + y(1) - liptot + y(2) + y(3) - y(4) + ...
         2* ptot)))/((2* kb7 + 2 *kcat7 + ka7* y(2)) *(kb3 *(kcat1 + ka1* y(1)) + ...
      ka3 *(kcat1 + DF* ka1* y(1)) *y(3) + kb1* (kb3 + ka3* y(3))));
      
LpP = (ka7 *y(2)* (kb1* kb3 *(y(1) - liptot + y(2) + y(3) - y(4) + 2 *ptot) + ...
     kb3 *kcat1 *(y(1) - liptot + y(2) + y(3) - y(4) + 2* ptot) + ...
     ka1 *kb3* y(1) *(ktot + y(1) - liptot + y(2) + y(3) - y(4) + 2* ptot) + ...
     ka3 *kb1 *y(3) *(ktot + y(1) - liptot + y(2) + y(3) - y(4) + 2* ptot) + ...
     ka3* kcat1 *y(3)* (ktot + y(1) - liptot + y(2) + y(3) - y(4) + 2 *ptot) + ...
     DF *ka1* ka3* y(1) *y(3)* (2* ktot + y(1) - liptot + y(2) + y(3) - y(4) + ...
        2 *ptot)))/((2 *kb7 + 2 *kcat7 + ka7* y(2))* (kb3* (kcat1 + ka1* y(1)) + ...
     ka3 *(kcat1 + DF* ka1* y(1))* y(3) + kb1* (kb3 + ka3* y(3))));

 %set dLK=0 and dLpAKL=0, solve for LK and LpAK. 
LK = (ka1 *y(1)* (kb1 *(ktot + y(1) - liptot + y(2) + y(3) - y(4) - LpP - 2 *P + ...
         2 *ptot) + ...
      kcat1 *(ktot + y(1) - liptot + y(2) + y(3) - y(4) - LpP - 2* P + ...
         2 *ptot) + ...
      DF *ka1* y(1) *(2 *ktot + y(1) - liptot + y(2) + y(3) - y(4) - LpP - 2 *P + ...
         2* ptot)))/((kb1 + kcat1)^2 + 2 *DF *ka1 *(kb1 + kcat1)* y(1) + ...
    DF *ka1^2* y(1)^2);

LpAK = -(((kb1 + ...
        kcat1)* (kb1 *(y(1) - liptot + y(2) + y(3) - y(4) - LpP - 2 *P + ...
           2* ptot) + ...
        kcat1* (y(1) - liptot + y(2) + y(3) - y(4) - LpP - 2 *P + 2 *ptot) + ...
        ka1* y(1)* (ktot + y(1) - liptot + y(2) + y(3) - y(4) - LpP - 2* P + ...
           2 *ptot)))/((kb1 + kcat1)^2 + 2* DF *ka1* (kb1 + kcat1) *y(1) + ...
      DF *ka1^2* y(1)^2));

%These 4 variable definitions come from conservation of mass.
LpAPLp= -y(4) - LpP - P + ptot; 
LpAKL = 0.5*(-y(1) + liptot - LK - y(2) - y(3) - LpAK + y(4) + LpP + 2*P - 2*ptot);
A= 0.5*(2*atot + y(1) - liptot + LK + y(2) - y(3) - LpAK - y(4) + LpP);
K = 0.5* (2*ktot + y(1) - liptot - LK + y(2) + y(3) - LpAK - y(4) - LpP - 2*P + 2*ptot);

%THE 8 VARIABLE VERSION IS EXACT, ELIMINATE 4 VARS USING MASS CONSERVATION
% dy(1) = kb1*y(2) + kb1*LpAKL + kcat7*LpAPLp + kcat7*y(7) - ka1*K*y(1) - ka1*DF*y(1)*y(5);%L
% 
% dy(2) = ka1*K*y(1) - kb1*y(2) - kcat1*y(2);%LK
% dy(3) = kcat1*y(2) + kb2*y(4) + kcat1*LpAKL + kb7*LpAPLp + kb7*y(7) - ka2*A*y(3) - ka7*y(3)*y(6) - ka7*DF*y(3)*y(8);%Lp
% 
% dy(4) = kb3*y(5) + kb4*y(8) + ka2*A*y(3) - kb2*y(4) - ka3*K*y(4) - ka4*y(4)*y(6);%LpA
% dy(5) = kb1*LpAKL + kcat1*LpAKL + ka3*K*y(4) - kb3*y(5) - ka1*DF*y(1)*y(5);%LpAK
% 
% dy(6) = kb4*y(8) + kb7*y(7) + kcat7*y(7) - ka7*y(3)*y(6) - ka4*y(4)*y(6);%P
% dy(7) = ka7*y(3)*y(6) - kb7*y(7) - kcat7*y(7);%LpP
% dy(8) = kb7*LpAPLp + kcat7*LpAPLp + ka4*y(4)*y(6) - kb4*y(8) - ka7*DF*y(3)*y(8);%LpAP

dy(1) = kb1*LK + kb1*LpAKL + kcat7*LpAPLp + kcat7*LpP - ka1*K*y(1) - ka1*DF*y(1)*LpAK;%L

dy(2) = kcat1*LK + kb2*y(3) + kcat1*LpAKL + kb7*LpAPLp + kb7*LpP - ka2*A*y(2) - ka7*y(2)*P - ka7*DF*y(2)*y(4);%Lp

dy(3) = kb3*LpAK + kb4*y(4) + ka2*A*y(2) - kb2*y(3) - ka3*K*y(3) - ka4*y(3)*P;%LpA
dy(4) = kb7*LpAPLp + kcat7*LpAPLp + ka4*y(3)*P - kb4*y(4) - ka7*DF*y(2)*y(4);%LpAP



end



%Given the unknowns, calculate all the dependent variables.
function[A, K, LpAKL, LpAPLp, LpAK, LK, P, LpP]=calc_other_vars(y, parmVec,mass)
ktot=mass.ktot;
ptot=mass.ptot;
atot=mass.atot;
liptot=mass.liptot;
ka1=parmVec.ka1;
kb1=parmVec.kb1;
kcat1=parmVec.kcat1;
ka2=parmVec.ka2;
kb2=parmVec.kb2;
ka3=parmVec.ka3;
kb3=parmVec.kb3;
ka4=parmVec.ka4;
kb4=parmVec.kb4;
ka7=parmVec.ka7;
kb7=parmVec.kb7;
kcat7=parmVec.kcat7;
DF=parmVec.DF;%this is the DF

%set dLpP=0 and dK=0, solve for P and LpP
P =((kb7 + kcat7)* (kb1 *kb3* (y(1) - liptot + y(2) + y(3) - y(4) + 2* ptot) +...
       kb3* kcat1* (y(1) - liptot + y(2) + y(3) - y(4) + 2 *ptot) + ...
      ka1 *kb3* y(1)* (ktot + y(1) - liptot + y(2) + y(3) - y(4) + 2* ptot) + ...
      ka3 *kb1* y(3)* (ktot + y(1) - liptot + y(2) + y(3) - y(4) + 2* ptot) + ...
      ka3* kcat1 *y(3) *(ktot + y(1) - liptot + y(2) + y(3) - y(4) + 2* ptot) + ...
      DF* ka1* ka3 *y(1)* y(3)* (2 *ktot + y(1) - liptot + y(2) + y(3) - y(4) + ...
         2* ptot)))/((2* kb7 + 2 *kcat7 + ka7* y(2)) *(kb3 *(kcat1 + ka1* y(1)) + ...
      ka3 *(kcat1 + DF* ka1* y(1)) *y(3) + kb1* (kb3 + ka3* y(3))));
      
LpP = (ka7 *y(2)* (kb1* kb3 *(y(1) - liptot + y(2) + y(3) - y(4) + 2 *ptot) + ...
     kb3 *kcat1 *(y(1) - liptot + y(2) + y(3) - y(4) + 2* ptot) + ...
     ka1 *kb3* y(1) *(ktot + y(1) - liptot + y(2) + y(3) - y(4) + 2* ptot) + ...
     ka3 *kb1 *y(3) *(ktot + y(1) - liptot + y(2) + y(3) - y(4) + 2* ptot) + ...
     ka3* kcat1 *y(3)* (ktot + y(1) - liptot + y(2) + y(3) - y(4) + 2 *ptot) + ...
     DF *ka1* ka3* y(1) *y(3)* (2* ktot + y(1) - liptot + y(2) + y(3) - y(4) + ...
        2 *ptot)))/((2 *kb7 + 2 *kcat7 + ka7* y(2))* (kb3* (kcat1 + ka1* y(1)) + ...
     ka3 *(kcat1 + DF* ka1* y(1))* y(3) + kb1* (kb3 + ka3* y(3))));

  %set dLK=0 and dLpAKL=0
LK = (ka1 *y(1)* (kb1 *(ktot + y(1) - liptot + y(2) + y(3) - y(4) - LpP - 2 *P + ...
         2 *ptot) + ...
      kcat1 *(ktot + y(1) - liptot + y(2) + y(3) - y(4) - LpP - 2* P + ...
         2 *ptot) + ...
      DF *ka1* y(1) *(2 *ktot + y(1) - liptot + y(2) + y(3) - y(4) - LpP - 2 *P + ...
         2* ptot)))/((kb1 + kcat1)^2 + 2 *DF *ka1 *(kb1 + kcat1)* y(1) + ...
    DF *ka1^2* y(1)^2);
LpAK = -(((kb1 + ...
        kcat1)* (kb1 *(y(1) - liptot + y(2) + y(3) - y(4) - LpP - 2 *P + ...
           2* ptot) + ...
        kcat1* (y(1) - liptot + y(2) + y(3) - y(4) - LpP - 2 *P + 2 *ptot) + ...
        ka1* y(1)* (ktot + y(1) - liptot + y(2) + y(3) - y(4) - LpP - 2* P + ...
           2 *ptot)))/((kb1 + kcat1)^2 + 2* DF *ka1* (kb1 + kcat1) *y(1) + ...
      DF *ka1^2* y(1)^2));



LpAPLp= -y(4) - LpP - P + ptot; 
LpAKL = 0.5*(-y(1) + liptot - LK - y(2) - y(3) - LpAK + y(4) + LpP + 2*P - 2*ptot);
A= 0.5*(2*atot + y(1) - liptot + LK + y(2) - y(3) - LpAK - y(4) + LpP);
K = 0.5* (2*ktot + y(1) - liptot - LK + y(2) + y(3) - LpAK - y(4) - LpP - 2*P + 2*ptot);

end

