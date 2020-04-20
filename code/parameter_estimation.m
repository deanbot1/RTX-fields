%% parameter_estimation.m estimates ADCC parameters...
%  
%
% D. Bottino, Takeda Pharmaceuticals

clear all; close all
parameter_filename = 'adcx_parameters_experiments.csv'; % special parameter and experimental settings filename
results_fileroot = 'adcx_parameter_results';


%% clear the workspace and close all figure windows

Tbig = readtable(parameter_filename);

ires = find(strcmp(Tbig.name,'RESERVED'));
Tpar = Tbig(1:ires-1,:);
Tpar.default = arrayfun(@(s)str2num(s{1}),Tpar.default);

vnames = Tpar.Properties.VariableNames;
EE = arrayfun(@(s)sscanf(s{1},'E%d_*'),vnames,'UniformOutput',false);
EI = find(arrayfun(@(s)~isempty(s{1}),EE));
Nexp = length(EI);

Texp = Tbig(ires+1:end,:);
Texp.Properties.RowNames = Texp.name;

for j = 1:Nexp
	i = EI(j);
	ename = vnames{i};
	expt(j).name = ename;
	deffun = @(rname)ifelseval(isempty(Texp{rname,ename}{1}),Texp{rname,'default'},Texp{rname,ename});
	expt(j).Xname = deffun('XNAME');
	expt(j).Ynames = {deffun('YNAME')};
	datafile = fullfile(deffun('DATAPATH'),deffun('DATAFILE'));
	dat = readtable(datafile);
	expt(j).xval = dat.(deffun('XHEADER'));
	expt(j).obs = dat.(deffun('YHEADER'));
	expt(j).model = str2func(deffun('MODEL'));
	expt(j).weight = str2num(deffun('WEIGHT'));
end

% now build par in the older way and define defpmap? No, may want to bypass
% old way, which involves parsing strings...or maybe yes... ugh
clear pxform
clear pinit cvs % initial parameters and CV structure

for j = 1:Nexp
	clear pmap
	ej = expt(j);
	for i = 1:height(Tpar)
		pname = Tpar.name{i};
		if Tpar.cv(i) < eps
			% it's fixed, no fitting
			pmap.(pname) = @(pfit)ifelseval(isempty(Tpar.(ej.name){i}),Tpar.default(i),str2num(Tpar.(ej.name){i})); % it's just a numeric value, fixed
			pinit.(pname) = Tpar.default(i);
			pxform.(pname) = Tpar.xform{i};
			cvs.(pname) = Tpar.cv(i);
		elseif isempty(Tpar.(ej.name){i}) % then it's fit "globally' ie one value across all expts
			pmap.(pname) = @(pfit)pfit.(pname);
			pinit.(pname) = Tpar.default(i);
			pxform.(pname) = Tpar.xform{i};
			cvs.(pname) = Tpar.cv(i);
		else
			if isempty(str2num(Tpar.(ej.name){i}))
				% then it's a string and needs to be parsed out somehow
				tname = [pname '_' Tpar.(ej.name){i}];
				pmap.(pname) = @(pfit)pfit.(tname);
				pinit.(tname) = Tpar.default(i);
				pxform.(tname) = Tpar.xform{i};
				cvs.(tname) = Tpar.cv(i);
			else
				% then it's a value and it's fixed somehow
				pmap.(pname) = @(pfit)str2num(Tpar.(ej.name){i}); % fixed to # in expt
				pinit.(pname) = feval(pmap.(pname)); % hopefully this doesn't get used anywhere
				pxform.(pname) = Tpar.xform{i};
				cvs.(pname) = Tpar.cv(i);
			end
		end
	end
	expt(j).pmap = pmap;
end



%defpmap = default_pmap(par);



%% plot goodness of fit of initial guesses in pinit

figure('Position',[239   558   990   420]);
xspan = 10.^[-2:0.25:6]';
plot_expt2(expt,pinit,xspan,'Xscale','log','Ylim',[-10 100],'Xtick',10.^[-2:2:6]);

	
%% ok now setup parameter estimation problem

errfun = @(Ypred,Yobs)sum(sum((Ypred(~isnan(Yobs))-Yobs(~isnan(Yobs))).^2)); % sum squared error function ignoring NaNs
%[pbig0,prow,pcol] = psquash(fxf(par.value),par.fit,Ne); % for initial guesses we use value field in the par table
pvec0 = pstruct2vec(pinit,pxform);
cvec = cvstruct2vec(cvs,pinit,pxform);
ii = cvec > eps & cvec < Inf;

%ofun = @(p)(objfun2(p,expt,pxform,errfun)); % single parameter vector objective function in transformed space
ofun = @(p)objfun2(p,expt,pxform,errfun) + sum(((p(ii)-pvec0(ii))./cvec(ii)).^2); % same objective function + bayesian penalty


% run fminsearch local optimizer using best guess from monte
% use at least 1000 maxiter when running a single experiment...
options = optimset('maxiter',500,'maxfunevals',10000,'Display','iter'); % set the options for your favorite optimizer
pbigbest = fminsearch(ofun,pvec0,options); % run your favorite optimizer
%pbigbest = fminsearch(ofun,pbigbest,options); % run your favorite optimizer again...
% save pbest
% load pbest.mat (if you want to skip parameter estimation)
pbest = pvec2struct(pbigbest,pxform)


%% plot the final results using the best-fit parameters
% that's it!
figure('Position',[239   558   990   420]);
xspan = 10.^[-2:0.25:6]';
plot_expt2(expt,pbest,xspan,'Xscale','log','Ylim',[-10 100],'Xtick',10.^[-2:2:6]);


%% make a fat results table
Tout = Tpar;
Npar = length(Tout.name);
for j = 1:Nexp
	clear ppp
	for k = 1:Npar
		pname = Tout.name{k};
		ppp(k,1) = feval(expt(j).pmap.(pname),pbest);
	end
	Tout.(expt(j).name) = ppp;
end

disp(Tout)
writetable(Tout,[results_fileroot,'_fat.csv']); % save to disk
disp(sprintf('%s written to: %s','Tout',[results_fileroot,'_fat.csv']));

%% save skinny results table to disk too
clear name cv initial final
name = sort(fieldnames(pbest));
for k = 1:length(name)
	cv(k,1) = cvs.(name{k});
	initial(k,1) = pinit.(name{k});
	final(k,1) = pbest.(name{k});
end

Tskinny = table(name,cv,initial,final);
disp(Tskinny)
writetable(Tskinny,[results_fileroot,'_skinny.csv']); % save to disk
disp(sprintf('%s written to: %s','Tskinny',[results_fileroot,'_skinny.csv']));



