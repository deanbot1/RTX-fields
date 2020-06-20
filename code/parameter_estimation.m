%% parameter_estimation.m estimates ADCC parameters...
%  
%
% D. Bottino, Takeda Pharmaceuticals

clear all; close all
parameter_filename = 'adcx_parameters_experiments.csv'; % special parameter and experimental settings filename
results_fileroot = 'adcx_parameter_results';
bweight = 1; % how much to weight bayesian penalty relative to chi squared error function, should be always 0 or 1, unless you have a very good reason to make it bigger!


%% clear the workspace and close all figure windows

[Tpar,Texp] = read_par_expt(parameter_filename);
[expt,pinit,pxform,cvs] = parse_par_expt(Tpar,Texp);

%% plot goodness of fit of initial guesses in pinit

figure('Position',[239   558   990   420]);
xspan = 10.^[-2:0.25:6]';
plot_expt2(expt,pinit,xspan,'Xscale','log','Ylim',[-10 100],'Xtick',10.^[-2:2:6]);
print('../out/init_fit.png','-dpng');
	
%% ok now setup parameter estimation problem and run it

errfun = @(Ypred,Yobs)sum(sum((Ypred(~isnan(Yobs))-Yobs(~isnan(Yobs))).^2)); % sum squared error function ignoring NaNs
pvec0 = pstruct2vec(pinit,pxform);
cvec = cvstruct2vec(cvs,pinit,pxform);
ii = cvec > eps & cvec < Inf;

%ofun = @(p)(objfun2(p,expt,pxform,errfun)); % single parameter vector objective function in transformed space
ofun = @(p)objfun2(p,expt,pxform,errfun) + bweight*sum(((p(ii)-pvec0(ii))./cvec(ii)).^2); % same objective function + bayesian penalty

% run fminsearch local optimizer
% use at least 1000 maxiter when running a single experiment...

options = optimset('maxiter',5000,'maxfunevals',20000,'Display','iter'); % set the options for your favorite optimizer
pbigbest = fminsearch(ofun,pvec0,options); % run your favorite optimizer
pbest = pvec2struct(pbigbest,pxform)


%% plot the final results using the best-fit parameters
% that's it!
figure('Position',[239   558   990   420]);
xspan = 10.^[-2:0.25:6]';
plot_expt2(expt,pbest,xspan,'Xscale','log','Ylim',[-10 100],'Xtick',10.^[-2:2:6]);
print('../out/final_fit.png','-dpng');

%% Generate outputs for identifiability analysis
save('pbest', 'pbest')

%% make a fat results table

Tout = save_fat_results(Tpar,pbest,expt,[results_fileroot,'_fat.csv']);

%% save skinny results table to disk too

Tskinny = save_skinny_results(pxform,cvs,pinit,pbest,[results_fileroot,'_skinny.csv']);

%%  make a plot of Tskinny 

figure;
param_plot(Tskinny,sprintf('Estimation results, bayes weight = %d',bweight))
print('../out/final_params.png','-dpng');
