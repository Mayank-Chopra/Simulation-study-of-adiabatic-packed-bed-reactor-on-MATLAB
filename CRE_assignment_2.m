Wspan = [0 100]; % Range for the independent variable W
Po = 10 ;% Initial Pressure (in atm)
y0 = [0;Po]; % Initial values for the dependent variables i.e X
vo = 20;% Initial Volumetric Flow rate (in L/s)
R = 8.314;% Gas constant 
To = 450;% Initial Temperature (in K)
Cao = (Po*101325 / (R *To*1000));% Initial Concentration (in mol/L)
FAo=Cao * vo;
Ea = 31.4 * 1000;%Activation energy for reaction(in J/mol)
alpha = 0.019;
[w,y]=ode45(@(w,y)ODEfun(y,FAo,Cao,Ea,R,alpha,To,Po),Wspan,y0);
tiledlayout(1,2)
set (gcf,'Position',[0,0.1,800,300])
z=size(y);
Temp_prof = zeros(z,"double");
for i=1:z(1,1)
    Temp_prof(i) = 450 + 500*y(i,1);
end
nexttile
plot(w,y(:,1));
ylabel('X');
xlabel('W(kg)');
axis([0 100 0 0.02])
title('Conversion Profile')
grid on
nexttile
plot(w,Temp_prof(:,1))
ylabel('Temperature(K)');
xlabel('W(kg)');
axis([0 100 450 455])
title('Temperature Profile')
grid on

function dYfuncvecdt = ODEfun(Yfuncvec,FAo,Cao,Ea,R,alpha,To,Po)
X = Yfuncvec(1);
P = Yfuncvec(2);
% Explicit equations
epsilon = (2-1/1);
T_prof = 450+500*X;%Got after energy balance on reactor
k=0.133*exp((Ea/R)*((1/450)-(1/T_prof)));
T_factor = To/T_prof;
p_in_term = (-alpha*(Po^2))/T_factor;
p_term_ergun = p_in_term/2*P;
P_factor = P/Po;
Ca = (Cao*(1-X)/(1+epsilon*X))*T_factor*P_factor; 
rt = -k*Ca;  
% Differential equations
dXdW = -rt/FAo; 
dPdW = p_term_ergun*(1+epsilon*X); 
dYfuncvecdt = [dXdW;dPdW]; 

end
