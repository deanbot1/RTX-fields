function [T,E,Estar,LDH,perf,CPX] = adcx(tf_mol,tf_et,nr_t_mol,nr_t_et,T0,E0toT0,Estar0,...
    g,r,kexp,gamma,...
    CD20,CD16,RTX,kon20,koff20,kon16,koff16,gamma_perf)

A0 = CD20;
C0 = CD16;
R0 = RTX;

% Run molecular rate equations to steady state
%f_rate = reactions(A0,R0,C0,kon20,kon16,koff20,koff16,tf_mol,nr_t_mol);
f_rate = SSadcc(RTX,true,CD16,CD20,kon16,koff16,kon20,koff20);

%% Effector-target code

% Initial condition
y0(1) = T0;  % T target cell
y0(2) = E0toT0*T0;  % E effector cell CD16 based on E:T ratio
y0(3) = 0;   % D dead cell
NE    = y0(2) + Estar0; % only E, no E* at first

pars = [g,r,kexp,NE,gamma,f_rate];

%% Simulate the Effector-Target ODEs
tspan = linspace(0,tf_et,nr_t_et);
[t,y] = ode23s(@et_code_rhs,tspan,y0,[],pars);

% Store variables
T = y(:,1);
E = y(:,2);
Estar = NE-E;
LDH   = y(:,3);
perf  = gamma_perf.*y(:,3);
CPX   = f_rate;

% %% Plotting
% figure(1)
% plot(t,y(:,1),'LineWidth',2)
% hold on
% plot(t,y(:,2),'LineWidth',2)
% plot(t,y(:,3),'LineWidth',2)
% ylim([0 5.5])
% legend('Target cells','Effector cells','Dead cells/lysis')
% xlabel('Time')
% ylabel('Cells/volume')
% set(gca,'FontSize',16)


% % Parameters f(R,CD16,CD20) function
% kon20  = 1;
% kon16  = 0.5;
% koff20 = 0.2;
% koff16 = 0.1;
% tf_mol = 200;
% nr_t_mol = 1000;

% % Parameters et_code
% g    = 0.01;
% r    = 5;
% kexp = 1;
% gamma = 0.5;
%
% T0 = 2;  % T target cell
% E0 = 5;  % E effector cell CD16
% CD20 = 20;
% CD16 = 30;
% RTX = 10^1;