%% estimates ADCC parameters
% This file looks for an ../io/<mfilename>.csv file for parameter settings
% And writes a file ../io/<mfilename>_results.csv for parameter estimates
% also writes some figures to ../io directory based on <mfilename>
%
% 2019-20 D. Bottino, Takeda Pharmaceuticals
% with open source improvements by V. Ciocanel (Duke) and K. Johnson (UT
% austin): specifically now nlsqnonlin instead of fminsearch

clear all; close all
mname = mfilename;
if strcmp(mname(1:8),'LiveEdit'), error('dont run this in chunks, bad things will happen'); end
parameter_filename = sprintf('../io/%s.csv',mname); % special parameter and experimental settings filename
results_fileroot   = sprintf('../io/%s_results',mname);
bweight = 1; % how much to weight bayesian penalty relative to chi squared error function, should be always 0 or 1, unless you have a very good reason to make it bigger!

%% clear the workspace and close all figure windows
[Tpar,Texp] = read_par_expt(parameter_filename);
[expt,pinit,pxform,cvs] = parse_par_expt(Tpar,Texp);

%% plot goodness of fit of initial guesses in pinit

figure('Position',get(0,'ScreenSize'))
rtxspan = 10.^[-2:0.25:6]';
plot_expt_grid(expt,pinit,1,rtxspan,'Xscale','log','Ylim',[-10 100],'Xtick',10.^[-2:2:6]);
print(sprintf('../io/%s_init_fit.png',mname),'-dpng');

	
%% ok now setup parameter estimation problem and run it
% errfun = @(Ypred,Yobs)sum(sum((Ypred(~isnan(Yobs))-Yobs(~isnan(Yobs))).^2)); % sum squared error function ignoring NaNs
% for lsqnonlin, we simply need to write an error fun that stacks the
% errors vertically.
errfunlong = @(Ypred,Yobs)(Ypred(~isnan(Yobs))-Yobs(~isnan(Yobs))); % vector of residuals
werrfunlong = @(Ypred,Yobs,Yerr)(Ypred(~isnan(Yobs))-Yobs(~isnan(Yobs)))./Yerr; % vector of residuals, inverslely weighted by standard deviation of measurement error

pvec0 = pstruct2vec(pinit,pxform);
cvec = cvstruct2vec(cvs,pinit,pxform);
ii = cvec > eps & cvec < Inf;

ofunlong = @(p)objfunlong(p,expt,pxform,werrfunlong, pvec0, cvec, bweight, ii); % returns vector of residuals // same objective function + bayesian penalty
finit = sum((ofunlong(pvec0)).^2);

% Run lsqnonlin optimizer
% *** NOTE *** the OptimalityTolerance and StepTolerance below may be too big, or could be
% deleted, to make lsqnonlin work a little harder...
options = optimoptions(@lsqnonlin,'MaxIterations',10000,'MaxFunctionEvaluations',20000,'Display','iter','OptimalityTolerance',0.001,'StepTolerance',0.001); % set the options for your favorite optimizer
% options = optimoptions(@lsqnonlin,'MaxIterations',10000,'MaxFunctionEvaluations',20000,'Display','iter',...
% 	'OptimalityTolerance',0.001,'StepTolerance',0.001,'UseParallel',true);
[pbestlong,resnorm] = lsqnonlin(ofunlong,pvec0,[],[],options);
pbest = pvec2struct(pbestlong,pxform);

%% Generate outputs for identifiability analysis
save(sprintf('../io/%s_pbest.mat',mname));

%% make a fat results table
Tout = save_fat_results(Tpar,pbest,expt,[results_fileroot,'_fat.csv']);

%% save skinny results table to disk too
Tskinny = save_skinny_results(pxform,cvs,pinit,pbest,[results_fileroot,'_skinny.csv']);

%% plot the final results using the best-fit parameters
% that's it!

figure('Position',get(0,'ScreenSize'))
rtxspan = 10.^[-2:0.25:6]';
plot_expt_grid(expt,pbest,1,rtxspan,'Xscale','log','Ylim',[-10 100],'Xtick',10.^[-2:2:6]);
print(sprintf('../io/%s_final_fit.png',mname),'-dpng');

% UNCOMMENT THIS WHEN YOURE SPANNING MULTIPLE E:T RATIOS
% who knows, it might work
% figure('Position',get(0,'ScreenSize'))
% etspan = 2.^[-1:.25:4]';
% plot_expt_grid(expt,pbest,2,etspan,'Xscale','log','Ylim',[-10 100],'Xtick',2.^[-1:1:4]);
% print(sprintf('../io/%s_final_fit2.png',mname),'-dpng');


%%  make a plot of Tskinny 
figure;
param_plot(Tskinny,sprintf('Estimation results, bayes weight = %d',bweight))
grid on;
print(sprintf('../io/%s_final_params.png',mname),'-dpng');



