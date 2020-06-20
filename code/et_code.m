clear
close all

% Parameters et_code
g    = 0.01;
r    = 5;
kexp = 1;

% Initial conditionclear

y0(1) = 2;  % T target cell
y0(2) = 5;  % E effector cell CD16
y0(3) = 0;  % D dead cell
NE    = y0(2); % only E, no E* at first

pars = [g,r,kexp,NE];

% Parameters f(R,CD16,CD20) function
kdRCD20 = 1;
kdRCD16 = 1;
gamma   = 0.5;
R       = 0.1;
CD20    = 0.1;
CD16    = 0.1;

pars_f(1) = kdRCD20;
pars_f(2) = kdRCD16;
pars_f(3) = gamma;
pars_f(4) = CD20;
pars_f(5) = CD16;
pars_f(6) = R;

pars = [pars pars_f];

%% Simulate the ODEs
tspan = 0:.1:100;
[t,y] = ode45(@et_code_rhs,tspan,y0,[],pars);

%% Plotting
figure(1)
plot(t,y(:,1),'LineWidth',2)
hold on
plot(t,y(:,2),'LineWidth',2)
plot(t,y(:,3),'LineWidth',2)
legend('T','E','D')
set(gca,'FontSize',16)

