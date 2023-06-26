Wspan = [0 1000]; % Range for the independent variable W
To = 450;% Initial Temperature (in K)
Po = 10 ;% Initial Pressure (in atm)
y0 = [0;To;Po]; % Initial values for the dependent variables i.e X
vo = 20;% Initial Volumetric Flow rate (in L/s)
R = 8.314;% Gas constant
Cao = (Po*101325 / (R *To*1000));% Initial Concentration (in mol/L)
FAo=Cao * vo;
Ea = 31.4 * 1000;%Activation energy for forward reaction(in J/mol)
alpha = 0.019;
Er = 51.4 * 1000;%Activation energy for backward reaction(in J/mol)
[w,y]=ode45(@(w,y)ODEfun(y,FAo,Cao,Ea,Er,R,alpha,To,Po),Wspan,y0);
tiledlayout(1,2)
set (gcf,'Position',[0,0.1,800,300])
nexttile
plot(w,y(:,1))
ylabel('X');
xlabel('W(kg)');
axis([0 100 0 0.02])
title('Conversion Profile')
grid on
nexttile
plot(w,y(:,2));
ylabel('Temperature(K)');
xlabel('W(kg)');
axis([0 100 450 455])
title('Temperature Profile')
grid on

function dYfuncvecdt = ODEfun(Yfuncvec,FAo,Cao,Ea,Er,R,alpha,To,Po)
X = Yfuncvec(1); 
T = Yfuncvec(2);
P = Yfuncvec(3);
% Explicit equations
epsilon = (2-1/1);
kf=0.133*exp((Ea/R)*((1/To)-(1/T)));
kr=exp((Er/R)*((1/To)-(1/T)));
T_factor = To/T;
p_in_term = (-alpha*(Po^2))/T_factor;
p_term_ergun = p_in_term/2*P;
P_factor = P/Po;
Ca = (Cao*(1-X)/(1+epsilon*X))*T_factor*P_factor;
Cb = (Cao*(X)/(1+epsilon*X))*T_factor*P_factor;
Cc = (Cao*(X)/(1+epsilon*X))*T_factor*P_factor;
rt = -kf*Ca+kr*Cb*Cc;  
jacket_term = 0.08;%Ua/Pb for reactor jacket = 0.08 (in J/s-kgcat-K)
Cpa = 40 ; %Molar heat capacity of reactant A (in J/mol-K)
DeltaHoRx = -20000 ;%Heat of Reaction ( in J/mol)
heat_gen = rt*DeltaHoRx;
heat_removal = jacket_term*(T-323);
heat_term = heat_gen - heat_removal;
heat_constant = FAo*Cpa;

% Differential equations
dXdW= -rt/FAo; 
dTdW = heat_term/heat_constant;
dPdW = p_term_ergun*(1+epsilon*X);
dYfuncvecdt = [dXdW;dTdW;dPdW]; 

end
