% Profile likelihoods calculation for ADCX model parameters

% Load in parameter
clear
close all
rng('shuffle')
parameter_filename = 'adcx_parameters_experiments.csv';   % special parameter and experimental settings filename
results_fileroot   = 'adcx_prof_likelihoods_center';
load('bootstrap_500runs_July25_pbestnew0_bw0.mat','pbigbootall','pxform')
pbigboot           = pbigbootall;

par         = 1; % parameter kept fixed
bweight     = 0; % how much to weight bayesian penalty relative to chi squared error function, should be always 0 or 1, unless you have a very good reason to make it bigger!
nr_runs     = 20;

% Grab the 99th percentiles from bootstrapping parameter ranges to sample from
bootbounds = prctile(pbigboot',[0.5,99.5], 1);
% Use the spread in the bootstrap parameter estiamtes to sample symmetrically
diamvec = bootbounds(2,:) - bootbounds(1,:);
diamvec = diamvec';

% Load pbest parameters
pbestload  = load('pbest.mat'); % 11 parameters optimized
pbest      = pbestload.pbest; % structure
% Load pbest from 7-parameter search
pbest0load = load('pbest_bw0.mat'); % 7 parameters optimized, bw 0 here 
pbest0     = pbest0load.pbestnew;

pbestlog   = pstruct2vec(pbest0,pxform); % vector, in log form
fval       = pbest0load.fval;
factor     = 0.4;
phatbounds = horzcat((pbestlog-factor*diamvec), (pbestlog+factor*diamvec));
phatbounds = phatbounds';

[Tpar,Texp] = read_par_expt(parameter_filename); % usual parameter file
list_names  = fields(pxform);
obj         = zeros(length(list_names),nr_runs); % objective function value/fval
pj          = zeros(length(list_names),nr_runs); % objective function value/fval

%% Profile likelihoods
for k = par  % length(list_names)
    
    disp(['Parameter ' list_names(k)])
    list_pars   = setdiff(1:length(list_names),k);
    values_boot = pbigbootall(k,:);
    pj(k,:)     = linspace(phatbounds(1,k),phatbounds(2,k),nr_runs);
%     pj(k,:) = linspace(min(values_boot),max(values_boot),nr_runs);
    
    for l = 1:nr_runs
        l
        [expt,pinit,pxform,cvs] = parse_par_expt(Tpar,Texp,'g_Z138',pbest.g_Z138,...
        'koff20',pbest.koff20,'g_SUDHL4',pbest.g_SUDHL4,'koff16_F158',pbest.koff16_F158,...
        list_names{k},pj(k,l)); % list_pars(k)

        % Run parameter estimation
        errfun = @(Ypred,Yobs)sum(sum((Ypred(~isnan(Yobs))-Yobs(~isnan(Yobs))).^2)); % sum squared error function ignoring NaNs
        pvec0  = pstruct2vec(pinit,pxform); % pvec0 from pinit, Excel spreadsheet
%         % Correct pvec0 so that it takes values from pbestlog
%         pvec0  =  pbestlog(list_pars); % pvec0 from pbestlog
        cvec   = cvstruct2vec(cvs,pinit,pxform);
        ii     = cvec > eps & cvec < Inf;
        
        %ofun = @(p)(objfun2(p,expt,pxform,errfun)); % single parameter vector objective function in transformed space
        ofun   = @(p)objfun2(p,expt,pxform,errfun) + bweight*sum(((p(ii)-pvec0(ii))./cvec(ii)).^2); % same objective function + bayesian penalty
        
        % run fminsearch local optimizer
        % options = optimset('maxiter',10000,'maxfunevals',20000); %,'Display','iter' % set the options for your favorite optimizer
        options = optimset('maxiter',10000,'maxfunevals',20000,'Display','iter'); 
        pvec0log = pbestlog(list_pars); % this starts from pbest0; is that right?
        [pbigbest,obj_best] = fminsearch(ofun,pvec0log,options); % run your favorite optimizer
        pestimate = pvec2struct(pbigbest,pxform);
        obj(k,l)  = obj_best;       
    end
end

%% Save results
if ~exist(results_fileroot, 'dir')
   mkdir(results_fileroot)
end
save([results_fileroot '/prof_likelihoods_bw' num2str(bweight)...
    '_runs' num2str(nr_runs) '_par' num2str(par) '.mat'])

%% Plotting
threshold = chi2inv(0.95,length(list_names))/2 + fval;
for k = par
    
figure(k)
hold on
plot(pj(k,:),obj(k,:),'*-','LineWidth',2)
plot(pbestlog(k),fval,'r*','LineWidth',2)
plot([pj(k,1) pj(k,end)],[threshold threshold],'r--','LineWidth',2)
% ylim([800 1000])
xlabel(list_names{k})
ylabel('fval')
set(gca,'FontSize',16)

end




