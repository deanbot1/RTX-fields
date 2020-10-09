%identifiability.m for ADCX model parameters

% Load in parameter 
clear; close all; 
rng('shuffle')
parameter_filename = 'adcx_parameters_experiments.csv'; % special parameter and experimental settings filename
results_fileroot = 'adcx_parameter_results';
pbest = load('pbest.mat');
pbest = pbest.pbest;

index = 1;
[Tpar,Texp] = read_par_expt(parameter_filename);
[expt,pinit,pxform,cvs] = parse_par_expt(Tpar,Texp,'g_Z138',pbest.g_Z138,...
    'koff20',pbest.koff20,'g_SUDHL4',pbest.g_SUDHL4,'koff16_F158',pbest.koff16_F158);

%% Bootstrapping
% Find the residuals and model from the current best fitting parameters and 
% store them
Ne = length(expt);
Yexpt = [];
ttnames = fieldnames(expt(1).pmap); % names of target parameter names
pbig = pstruct2vec(pbest,pxform); % gather full model inputs as vector
% For pbest, generate model  
Ymod = modelfunskinny(pbig,expt, pxform);
for j = 1:Ne
    Yexpt = vertcat(Yexpt, expt(j).obs);
end
% residuals as single vector
residuals = Ymod-Yexpt;
ndata = length(residuals);
sigma = std(residuals);


%% Run a loop to creat and fit synthetic data
nruns   = 100;
nstarts = 1;
bweight = 1;
rand_index_matrix  = randi([1 ndata], [ndata,nruns]); % store
%%
tic
pbestlog    = pstruct2vec(pbest,pxform);
pbestactual = exp(pbestlog);
pvec0logall = zeros(length(pbestlog),nruns);
pbigboot    = zeros(length(pbestlog),nruns);
for i = 1:nruns
    disp(i)
    % create synthetic data by:
    % randomly sample with replacement from residual list and add to the
    % model fit from original pbest
    rand_int  = rand_index_matrix(:,i);
    Ysim = Ymod + residuals(rand_int);
    
    % Put simulated data  temporarily into expt.obs
    m = 1; % counter for place along 
    for k = 1:length(expt)
        nobs = length(expt(k).obs);
        expt(k).obs = Ysim(m:m+nobs-1);
        m    = m + nobs;
    end   
      
    % fit the synthetic data: (perhaps make this a function)
    errfun = @(Ypred,Yobs)sum(sum((Ypred(~isnan(Yobs))-Yobs(~isnan(Yobs))).^2))./(sigma.^2); % sum squared error function ignoring NaNs
    pvec0  = pstruct2vec(pinit,pxform);
    cvec   = cvstruct2vec(cvs,pinit,pxform);
    ii     = cvec > eps & cvec < Inf;  
    ofun   = @(p)objfun2(p,expt,pxform,errfun) + bweight*sum(((p(ii)-pvec0(ii))./cvec(ii)).^2); % same objective function + bayesian penalty
    
    % run fminsearch local optimizer
    % use at least 1000 maxiter when running a single experiment... 
    options  = optimset('maxiter',10000,'maxfunevals',20000); % set the options for your favorite optimizer
    % Create matrix of initial conditions by drawing from normal distribution
    pvec0all = normrnd(pbestactual,0.1*pbestactual);
    pvec0log = log(pvec0all);
    [pbooti,fval]    = fminsearch(ofun,pvec0log,options);
    pvec0logall(:,i) = pvec0log;
    
    % Store bootstrapped parameter estimates (in log space)
    pbigboot(:,i) = pbooti;

    % For pbest, generate model  
    Ymodboot = modelfunskinny(pbooti,expt, pxform); 

end
toc

save(['bootstrap' num2str(index) '.mat'])

% %% Plot one example to make sure the fitting is working
% figure;
% plot(1:1:ndata, Ysim, 'r*')
% hold on
% plot(1:1:ndata, Ymodboot, 'r-')
% plot(1:1:ndata, Ymod, 'b-')
% plot(1:1:ndata, Yexpt, 'b*')
% xlim([0, ndata]) 
% xlabel('[RTX]')
% ylabel('% ADCC')
% legend('simulated data','fit to sim data', 'pbest model', 'real data', 'Location', 'NorthWest')
% legend boxoff
% set(gca,'FontSize',20,'LineWidth',1.5)

% 
% %% Use plotmatrix to visualize parameters
% paramnames = fieldnames(pbest);
% figure;
% [S,AX,BigAx,H,HAx] = plotmatrix(pbigboot');
% title(BigAx,'A Comparison of Parameter Estimates ')
% for i = 1:length(paramnames)
% ylabel(AX(i,1),paramnames{i})
% xlabel(AX(9,i),paramnames{i})
% end
% 
% %% get back original expt structure with real data obs
% [expt,pinit,pxform,cvs] = parse_par_expt(Tpar,Texp);
% 
