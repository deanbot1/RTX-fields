%ODE model for rtx

%% Parameters before fitting
% kon20 = 0.003880961;
% kon16 = 64.71922;
% koff_RA = 0.0468;
% koff16 = 16.956;

%% Parameters after fitting
kon_RA = 7.82554959564354e-07;
kon_FR = 38.0110113102053;%V158
% kon16 = 44.862740158308;%F158
koff_RA = 0.0468;
koff_FR = 16.956;

%% Initial conditions for CD20 and CD16

CD20_IC = 28732.5519292874;% CD20_Z138 molecules per target cell
% CD20_IC = 285738.642756014;% CD20_SUDHL4 molecules per target cell
CD16_IC = 11480.3061574918;% V158 molecules per effector cell
% CD16_IC = 23012.3854351209;% F158 molecules per effector cell

Atot = CD20_IC;
Ftot = CD16_IC;

%% Computation by ode solver and analytical sol
R0=10.^(-4:0.1:8);%to check
% options = odeset('RelTol',1e-3,'AbsTol',1e-15);
options = odeset('RelTol',1e-3,'AbsTol',1e-6);%default values
tf = 20;
tspan = [0,tf];
a=zeros(1,length(R0));
C_f=zeros(1,length(R0));
for idx=1:length(R0)
    idx
    y0=[CD20_IC,CD16_IC,0, 0, 0];
    a(idx) = CPX(R0(idx),kon_RA,kon_FR,koff_RA,koff_FR,tspan,y0,options); % solver
    % analytical sol:
    C_f(idx) = Atot + ((koff_RA + R0(idx)*kon_RA)*(koff_RA*koff_FR - (Atot^2*R0(idx)^2*kon_RA^2*kon_FR^2 - 2*Atot*Ftot*R0(idx)^2*kon_RA^2*kon_FR^2 + 2*Atot*R0(idx)^3*kon_RA^2*kon_FR^2 + 2*Atot*R0(idx)^2*koff_RA*kon_RA*kon_FR^2 + 2*Atot*R0(idx)^2*koff_FR*kon_RA^2*kon_FR + 2*Atot*R0(idx)*koff_RA*koff_FR*kon_RA*kon_FR + Ftot^2*R0(idx)^2*kon_RA^2*kon_FR^2 + 2*Ftot*R0(idx)^3*kon_RA^2*kon_FR^2 + 2*Ftot*R0(idx)^2*koff_RA*kon_RA*kon_FR^2 + 2*Ftot*R0(idx)^2*koff_FR*kon_RA^2*kon_FR + 2*Ftot*R0(idx)*koff_RA*koff_FR*kon_RA*kon_FR + R0(idx)^4*kon_RA^2*kon_FR^2 + 2*R0(idx)^3*koff_RA*kon_RA*kon_FR^2 + 2*R0(idx)^3*koff_FR*kon_RA^2*kon_FR + R0(idx)^2*koff_RA^2*kon_FR^2 + 4*R0(idx)^2*koff_RA*koff_FR*kon_RA*kon_FR + R0(idx)^2*koff_FR^2*kon_RA^2 + 2*R0(idx)*koff_RA^2*koff_FR*kon_FR + 2*R0(idx)*koff_RA*koff_FR^2*kon_RA + koff_RA^2*koff_FR^2)^(1/2) + R0(idx)*koff_RA*kon_FR + R0(idx)*koff_FR*kon_RA + R0(idx)^2*kon_RA*kon_FR - Atot*R0(idx)*kon_RA*kon_FR + Ftot*R0(idx)*kon_RA*kon_FR))/(2*R0(idx)*kon_RA*(koff_RA*kon_FR + R0(idx)*kon_RA*kon_FR));

end

%% basic error
delta_curves = abs(a-C_f);
percentageDifference = delta_curves ./ a; % Percent by element
meanDifference = mean(percentageDifference); % Average percentage over all elements

%% Figure to compare    
figure(6)
plot(R0,a,'linewidth',4)
set(gca,'XScale','log')
xlabel('Rituximab')
ylabel('Complex')
hold on;
plot(R0,C_f,'linewidth',1)
hold off;
legend('solver','analytical');
text(1e-3,500,sprintf('Final time %d',tf))
text(1e-3,1000,sprintf('Mean percentage diff %d',meanDifference))




function ABC = CPX(R0,kon_RA,kon_FR,koff_RA,koff_FR,tspan,y0,options)
    [T,Y]=ode23s(@(t,y) model1(t,y,R0,kon_RA,kon_FR,koff_RA,koff_FR),tspan,y0,options);
    ABC = Y(end,5);
end

function dydt = model1(t,y,R0,kon_RA,kon_FR,koff_RA,koff_FR)

A=y(1);
% B=y(2);
F=y(2);
RA=y(3);
FR=y(4);
FRA=y(5);

dydt(1) = koff_RA*(RA + FRA) - kon_RA*(R0*A + FR*A);
dydt(2) = koff_FR*(FR + FRA) - kon_FR*(F*R0 + F*RA);
dydt(3) = kon_RA*R0*A - koff_RA*RA - kon_FR*F*RA + koff_FR*FRA;
dydt(4) = kon_FR*F*R0 - koff_FR*FR - kon_RA*FR*A + koff_RA*FRA;
dydt(5) = kon_FR*F*RA + kon_RA*FR*A - (koff_FR+koff_RA)*FRA;

% dydt(1)=-kon20*R0.*A+koff20*FR-kon20*A.*RA+koff20*FRA;
% dydt(2)=-kon16*R0.*F+koff16*R0.*F-kon16*F.*FR + koff16*FRA;
% dydt(3)=kon20*R0.*A-kon16*F.*FR-koff20*FR+koff16*FRA;
% dydt(4)=kon16*R0.*F - kon20*A.*RA-koff16*RA+koff20*FRA;
% dydt(5)=kon16*F.*FR + kon20*A.*RA-koff16*FRA-koff20*FRA;
dydt=[dydt(1);dydt(2);dydt(3);dydt(4);dydt(5)];

end
