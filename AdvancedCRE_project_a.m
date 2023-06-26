Wspan = [0 100]; % Range for the independent variable W
y0 = 0; % Initial values for the dependent variables i.e X

vo = 20;% Initial Volumetric Flow rate (in L/s)
Po = 10 * 101325;% Initial Pressure (in Pa)
R = 8.314;% Gas constant 
To = 450;% Initial Temperature (in K)
Cao = (Po / (R * To * 1000));% Initial Concentration (in mol/L)
FAo=Cao * vo;
Ea = 31.4 * 1000;%Activation energy for reaction(in J/mol)
[w,y]=ode45(@(w,y)ODEfun(y,FAo,Cao,Ea,R),Wspan,y0);
tiledlayout(1,2)
set (gcf,'Position',[0,0.1,800,300])
z=size(w);
Temp_prof = zeros(z,"double");
for i=1:length(y)
    Temp_prof(i) = 450 + 500*y(i);
end
nexttile
plot(w,y)
ylabel('X');
xlabel('W(kg)');
axis([0 100 0 1])
title('Conversion Profile')
grid on
nexttile
plot(w,Temp_prof)
ylabel('Temperature(K)');
xlabel('W(kg)');
axis([0 100 0 1000])
title('Temperature Profile')
grid on

function dYfuncvecdt = ODEfun(Y,FAo,Cao,Ea,R) 
X = Y(1); 
% Explicit equations
epsilon = (2-1/1);
T_prof = 450+500*X;%Got after energy balance on reactor
k=0.133*exp((Ea/R)*((1/450)-(1/T_prof)));
T_factor = 450/T_prof;
Ca = (Cao*(1-X)/(1+epsilon*X))*T_factor; 
rt = -k*Ca;  
% Differential equations
dXdW= -rt/FAo; 
dYfuncvecdt = dXdW; 

end

