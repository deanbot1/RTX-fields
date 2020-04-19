function [T,E,Estar,LDH,perf,CPX] = adcx(tspan,...
	T0,E0toT0,Estar0,g,r,kexp,gamma,...
    CD20,CD16,RTX,kon20,koff20,kon16,koff16,h,gammaPerf)

A0 = CD20;
C0 = CD16;
R0 = RTX;

% Run molecular rate equations to steady state
%f_rate = reactions(A0,R0,C0,kon20,kon16,koff20,koff16,tf_mol,nr_t_mol);
f_rate = SSadcc(RTX,CD16,CD20,kon16,koff16,kon20,koff20,h);

%% Effector-target code

% Initial condition
y0(1) = T0;  % T target cell
y0(2) = E0toT0*T0;  % E effector cell CD16 based on E:T ratio
y0(3) = 0;   % D dead cell
NE    = y0(2) + Estar0; % only E, no E* at first

pars = [g,r,kexp,NE,gamma,f_rate];

%% Simulate the Effector-Target ODEs
[t,y] = ode23s(@et_code_rhs,tspan,y0,[],pars);

% Store variables
T = y(:,1);
E = y(:,2);
Estar = NE-E;
LDH   = y(:,3);
perf  = gammaPerf.*y(:,3);
CPX   = f_rate;
