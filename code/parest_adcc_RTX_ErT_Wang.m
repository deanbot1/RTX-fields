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
parameter_filename = sprintf('../ioWang/%s.csv',mname); % special parameter and experimental settings filename
results_fileroot   = sprintf('../ioWang/%s_results',mname);
bweight = 1; % how much to weight bayesian penalty relative to chi squared error function, should be always 0 or 1, unless you have a very good reason to make it bigger!

%% clear the workspace and close all figure windows
[Tpar,Texp] = read_par_expt(parameter_filename);
[expt,pinit,pxform,cvs] = parse_par_expt(Tpar,Texp);

%% plot goodness of fit of initial guesses in pinit

% figure('Position',get(0,'ScreenSize'))
figure
rtxspan = 10.^[-2:0.25:6]';
localplotfun_RTX(expt,pinit,1,rtxspan,'Xscale','log','Ylim',[-10 100],'Xtick',10.^[-2:2:6]);
print(sprintf('../ioWang/%s_init_fit.png',mname),'-dpng');

% ErTspan = 1:2:50;
etspan = 2.^[-1:.25:5]';
localplotfun_ErT(expt,pinit,1,etspan,'Xscale','log','Ylim',[-10 100],'Xtick',2.^[-1:1:5]);
print(sprintf('../ioWang/%s_init_fit_ErT.png',mname),'-dpng');
% plot_expt_grid(expt,pinit,1,rtxspan,'Xscale','log','Ylim',[-10 100],'Xtick',10.^[-2:2:6]);
	
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
save(sprintf('../ioWang/%s_pbest.mat',mname));

%% make a fat results table
Tout = save_fat_results(Tpar,pbest,expt,[results_fileroot,'_fat.csv']);

%% save skinny results table to disk too
Tskinny = save_skinny_results(pxform,cvs,pinit,pbest,[results_fileroot,'_skinny.csv']);

%% plot the final results using the best-fit parameters
% that's it!

figure
rtxspan = 10.^[-2:0.25:6]';
localplotfun_RTX(expt,pbest,1,rtxspan,'Xscale','log','Ylim',[-10 100],'Xtick',10.^[-2:2:6]);
print(sprintf('../ioWang/%s_final_fit.png',mname),'-dpng');

etspan = 2.^[-1:.25:5]';
localplotfun_ErT(expt,pbest,1,etspan,'Xscale','log','Ylim',[-10 100],'Xtick',2.^[-1:1:5]);
print(sprintf('../ioWang/%s_final_fit_ErT.png',mname),'-dpng');

% % figure('Position',get(0,'ScreenSize'))
% figure
% rtxspan = 10.^[-2:0.25:6]';
% plot_expt_grid(expt,pbest,1,rtxspan,'Xscale','log','Ylim',[-10 100],'Xtick',10.^[-2:2:6]);
% print(sprintf('../ioWang/%s_final_fit.png',mname),'-dpng');

% %%
% % % UNCOMMENT THIS WHEN YOURE SPANNING MULTIPLE E:T RATIOS
% % % who knows, it might work - it does not work so far
% % % figure('Position',get(0,'ScreenSize'))
% etspan = 2.^[-1:.25:4]';
% plot_expt_grid(expt,pbest,2,etspan,'Xscale','log','Ylim',[-10 100],'Xtick',2.^[-1:1:4]);
% print(sprintf('../ioWang/%s_final_fit2.png',mname),'-dpng');


%%  make a plot of Tskinny 
figure;
param_plot(Tskinny,sprintf('Estimation results, bayes weight = %d',bweight))
grid on;
print(sprintf('../ioWang/%s_final_params.png',mname),'-dpng');



%% local plot function


function localplotfun_RTX(expt,pbest,xdim,xspan,varargin)

% plots experiments/model in a grid
% the columns are the experiments in expt
% the rows are different unique expt(j).xval(:,~xdim) values
% expt is a structure array indexed by experiment number
% pbest is a 'tall' parameter vector, maybe?
% xspan is requested simulation vector for first X variable (column xdim in
% expt(j).xval
% varargin is comma separated list of name,value pairs passed as graphics
% options, ie, 'linewidth','2', for the model curve.

Ne = length(expt);
ttnames = fieldnames(expt(1).pmap); % names of target paramter names

xvalbig = cat(1,expt.xval);
fdim = 3-xdim; % dimension along which to 'facet' plots
ufac = unique(xvalbig(:,fdim));
nfac = length(ufac); % number of vertical 'facets'; i.e, each E:T ratio
subplot2d = @(nrow,ncol,rowi,colj)subplot(nrow,ncol,ncol*(rowi-1)+mod(colj-1,ncol)+1); % 2d subplot

% First set of plots, varying RTX concentration
k = 5; % E:T ratio where we have multiple RTX concentrations data
for j = 1:Ne-1 % loop over experiments (only relevant ones)
        subplot2d(1,Ne-1,1,j);
        igood = find(expt(j).xval(:,fdim)==ufac(k));
        if ~isempty(igood)
            
            errorbar(expt(j).xval(igood,xdim),expt(j).obs(igood),expt(j).err(igood),'ro','markerfacecolor','r'); hold on
            
            for kk = 1:length(ttnames)
                tname = ttnames{kk};
                pfunc = expt(j).pmap.(tname);
                pstruct.(tname) = pfunc(pbest);
            end
        end
        xmat = zeros(length(xspan),2);
        xmat(:,xdim) = xspan;
        xmat(:,fdim) = ufac(k);
        Ypred = expt(j).model(pstruct,xmat);
        plot(xspan,Ypred,'b.-');
        set(gca,varargin{:});
        grid on
        
        if j==1
            ylabel(sprintf('%s,%s=%3.3g',char(expt(j).Ynames),expt(j).Xname{fdim},ufac(k)));
        end
        title(expt(j).name,'Interpreter','none');
    xlabel(expt(j).Xname{xdim});
end

end



function localplotfun_ErT(expt,pbest,xdim,xspan1,varargin)

% plots experiments/model in a grid
% the columns are the experiments in expt
% the rows are different unique expt(j).xval(:,~xdim) values
% expt is a structure array indexed by experiment number
% pbest is a 'tall' parameter vector, maybe?
% xspan is requested simulation vector for first X variable (column xdim in
% expt(j).xval
% varargin is comma separated list of name,value pairs passed as graphics
% options, ie, 'linewidth','2', for the model curve.

Ne = length(expt);
ttnames = fieldnames(expt(1).pmap); % names of target paramter names
xvalbig = cat(1,expt.xval);
% subplot2d = @(nrow,ncol,rowi,colj)subplot(nrow,ncol,ncol*(rowi-1)+mod(colj-1,ncol)+1); % 2d subplot

figure
% Second set of plots, varying E:T ratio at constant RTX concentration
j    = Ne;
xdim = 2;
fdim = 3-xdim; % dimension along which to 'facet' plots
ufac1 = expt(j).xval(1,fdim);
% subplot2d(1,Ne,1,j);
errorbar(expt(j).xval(:,2),expt(j).obs,expt(j).err,'ro','markerfacecolor','r'); hold on
for kk = 1:length(ttnames)
    tname = ttnames{kk};
    pfunc = expt(j).pmap.(tname);
    pstruct.(tname) = pfunc(pbest);
end
xmat1 = zeros(length(xspan1),2);
xmat1(:,xdim) = xspan1;
xmat1(:,fdim) = ufac1;
Ypred1 = expt(j).model(pstruct,xmat1);
plot(xspan1,Ypred1,'b.-');
set(gca,varargin{:});
grid on
ylabel(sprintf('%s,%s=%3.3g',char(expt(j).Ynames),expt(j).Xname{fdim},ufac1));
title(expt(j).name,'Interpreter','none');

end

%%
% function localplotfun(expt,pinit)
% 
% figure;
% subplot2d = @(nrow,ncol,rowi,colj)subplot(nrow,ncol,ncol*(rowi-1)+mod(colj-1,ncol)+1); % 2d subplot
% 
% Ne = length(expt);
% for j = 1:Ne
%     ttnames = fieldnames(expt(j).pmap);
%     for kk = 1:length(ttnames)
%         tname = ttnames{kk};
%         pfunc = expt(j).pmap.(tname);
%         pstruct.(tname) = pfunc(pinit);
%     end
%     
%     
%     for k = 1:length(expt(j).xhead)
%         x.(expt(j).xhead{k})=k;
%     end
%     
%     Uk = unique(expt(j).xval(:,x.EtoT));	Nk = length(Uk);
%     for k = 1:Nk
%         subplot2d(Nk,Ne,k,j);
%         igood = find(expt(j).xval(:,x.EtoT)==Uk(k));
%         Ui = unique(expt(j).xval(igood,x.RTXngml)); Ni = length(Ui); %was: x.RTXuM
%         colors = 'rbmgck';
%         for i = 1:Ni
%             igreat = intersect(igood,find(expt(j).xval(:,x.RTXngml)==Ui(i))); %was: x.RTXuM
%             h(i)=errorbar(expt(j).xval(igreat,x.TAK981uM),expt(j).obs(igreat),expt(j).err(igreat),'o','Color',colors(i)); hold on;
%             leg{i} = sprintf('RTX=%3.3g',Ui(i));
%             TAK981concs = 10.^[-2:.1:0]';
%             xsim = [TAK981concs,Uk(k)*ones(size(TAK981concs)),Ui(i)*ones(size(TAK981concs))];
%             ysim = expt(j).model(pstruct,xsim);
%             plot(xsim(:,x.TAK981uM),ysim,'-','LineWidth',2,'Color',colors(i));
%         end
%         title(sprintf('E:T = %3.3g',Uk(k)));
%         set(gca,'Ylim',[0 100]);
%         set(gca,'Xscale','log');
%         grid on
%     end
%     xlabel('[TAK-981]');
% end
% 
% end



