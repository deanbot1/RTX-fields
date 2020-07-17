%% Analyze parameter estimates from bootstrap runs

clear
close all

% % data from 100 runs (each) with 10% standard deviation and 1 start
% pbigbootall = [];
% for index = 1:5
%     load(['results_bootstrap_July16/bootstrap' num2str(index) '.mat'])
%     pbigbootall = [pbigbootall pbigboot];
% end
% save('bootstrap_500runs_July16.mat')

load('bootstrap_500runs_July16.mat')

%% Plot one example to make sure the fitting is working
figure()
plot(1:1:ndata, Ysim, 'r*')
hold on
plot(1:1:ndata, Ymodboot, 'r-')
plot(1:1:ndata, Ymod, 'b-')
plot(1:1:ndata, Yexpt, 'b*')
xlim([0, ndata])
xlabel('[RTX]')
ylabel('% ADCC')
legend('simulated data','fit to sim data', 'pbest model', 'real data', 'Location', 'NorthWest')
title('Example fit')
legend boxoff
set(gca,'FontSize',16,'LineWidth',1.5)

%% Use plotmatrix to visualize parameters
% Get back original expt structure with real data obs
[expt,pinit,pxform,cvs] = parse_par_expt(Tpar,Texp,'g_Z138',pbest.g_Z138,...
    'koff20',pbest.koff20,'g_SUDHL4',pbest.g_SUDHL4,'koff16_F158',pbest.koff16_F158);
paramnames = fieldnames(pinit); % used to be pbest
figure()
[S,AX,BigAx,H,HAx] = plotmatrix(pbigbootall');
title(BigAx,'A Comparison of Parameter Estimates ')
for i = 1:length(paramnames)
    ylabel(AX(i,1),paramnames{i})
    xlabel(AX(length(pbestactual),i),paramnames{i})
end
set(gca,'FontSize',16,'LineWidth',1.5)

% %% get back original expt structure with real data obs
% [expt,pinit,pxform,cvs] = parse_par_expt(Tpar,Texp);

%% Confidence intervals

% 99% CI from bootstrapping
phatranges = prctile(pbigboot,[0.5,99.5], 2);

yMean = mean(pbigboot,2);              % mean of each parameter value
ySEM  = std(pbigboot,0,2)/sqrt(nruns); % standard error of the mean for each parameter
CI95  = tinv([0.05 0.95],nruns-1);     % 95% probability intervals of t-distribution
yCI95 = bsxfun(@times,ySEM',CI95(:))'; % 95% confidence intervals for each parameter

%%
figure()
plot(yMean,'Linewidth',2)                                      % Plot Mean Of All Experiments
hold on
plot(yMean+yCI95,'Linewidth',2)                                % Plot 95% Confidence Intervals Of All Experiments
hold off
grid
set(gca,'FontSize',16,'LineWidth',1.5)

%%
figure()
errorbar(yMean,yCI95(:,2),'*','Linewidth',1.5)
xlabel('Parameter index')
ylabel('log(parameter)')
set(gca,'FontSize',16,'LineWidth',1.5)
