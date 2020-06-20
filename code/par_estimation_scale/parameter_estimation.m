%% parameter_estimation.m estimates ADCC parameters...
% rewriting this...
% Additionally this toolbox allows user to specify a transformation for
% each parameter:
%
% *For a strictly positive parameter, 'log' is recommended
% *For a fraction between 0 and 1, 'logit' is recommended
% *For a real number, just use 'identity' or any other string that doesn't
% match the two supported strings, and it will be ignored.
%
% lower and upper bounds are not explicitly handled by this package, maybe
% some other day, but they are enforced as a consequence of the
% transformation chosen above
%
% D. Bottino, Takeda Pharmaceuticals

%% clear the workspace and close all figure windows
clear all; close all
scale_par = 1;
if scale_par == 1
    parameter_filename = 'adcx_parameter_table_Kd_scaled.csv';
else
    parameter_filename = 'adcx_parameter_table_Kd.csv';
end
par = readtable(parameter_filename);
disp(par);

defpmap = default_pmap(par);

%% boh
% specify the model.
% a model is defined as a function taking in parameter and dependent
% variable values (say a time vector) and spitting out Y values for each
% dependent variable value. The model can either be an anonymous function
% or it can be its own standalone function. In this case I'm using a model
% that takes in a parameter STRUCTURE with fields that are parameter names.
% If your model already expects a parameter column vector and 'knows' the
% order of the parameters, you can do it that way too in which case you
% don't need the v2s function above.

% now we generate the experiment structure
% which in real life would be done by loading in experimental observations,
% not using the model to simulate the data, but this is a demo
data_path = '../../data/';
data_filenames = {'Herter_4A_V158onZ138.csv','Herter_4B_V158onSU-DHL4.csv'...
    'Herter_4C_F158onZ138.csv','Herter_4D_F158onSU-DHL4.csv', ...
    'Herter_4E_V158onZ138_High_Concentration.csv'};
Ne = length(data_filenames); % number of distinct experiments to fit to
%exptnames = {'exptA','exptC','exptZ'}; % labels for experiments
clear pdef
pj = 0; % p vector index thingy
for i = 1:Ne
    temp = readtable([data_path data_filenames{i}]);
    expt(i).name = strrep(data_filenames{i}(1:end-4),'-','_');
    iscore = strfind(expt(i).name,'_'); iscore = iscore(2);
    expt(i).name = {expt(i).name(1:iscore-1),expt(i).name(iscore+1:end)};
    expt(i).Ynames = {'%ADCC'}; % these can be different for each experiment but it's good practice to have all of them for all experiments, if possible. The names are matched for parameter estimation, so don't make any spelling variations!
    expt(i).Xname = '[RTX]';
    expt(i).model = @(pstruct,R_conc)adcx_wrapper(pstruct,R_conc);
    expt(i).xval = temp.Var1;
    expt(i).obs = temp.Var2;
    expt(i).pmap = defpmap;
    
    % this is where the magic happens == mapping of nonunique parameters to
    % parameters
    ion = strfind(expt(i).name{2},'on'); ion = ion(1); % first occurence
    switch expt(i).name{2}(ion+2:ion+3)
        case 'Z1'
            cellline = 'Z138';
        case 'SU'
            cellline = 'SUDHL4';
    end
    gname = ['g_' cellline];
    expt(i).pmap.g = @(pfit)pfit.(gname);
    cd20name = ['CD20_' cellline];
    expt(i).pmap.CD20 = @(pfit)pfit.(cd20name);
    CD16SNP = expt(i).name{2}(ion-4:ion-1);
    CD16name = ['kon16_' CD16SNP];
    expt(i).pmap.kon16 = @(pfit)pfit.(CD16name);
    % 	CD16name = ['CD16_' CD16SNP];
    % 	expt(i).pmap.CD16 = @(pfit)pfit.(CD16name);
end

% now we loop through the expt structure and build our pfit vector
% structure mapping

k = 0;
for j = 1:height(par)
    if par.fit(j) ~= 0
        k = k+1;
        pinit.(par.paramnames{j}) = par.value(j);
        pxform.(par.paramnames{j}) = par.xform{j};
    end
end

% figure('Position',[239   558   990   420]);
if scale_par == 1
    xspan = 10.^[-2:0.25:6]'.*50.*pi; % Span of RTX concentration
else
    xspan = 10.^[-2:0.25:6]'; % Span of RTX concentration
end
plot_expt2(expt,pinit,xspan,'Xscale','log','Ylim',[-10 100],'Xtick',10.^[-2:2:6]);


%% ok now setup parameter estimation problem

errfun = @(Ypred,Yobs)sum(sum((Ypred(~isnan(Yobs))-Yobs(~isnan(Yobs))).^2)); % sum squared error function ignoring NaNs
%[pbig0,prow,pcol] = psquash(fxf(par.value),par.fit,Ne); % for initial guesses we use value field in the par table
ofun  = @(p)(objfun2(p,expt,pxform,errfun)); % single parameter vector objective function in transformed space
pvec0 = pstruct2vec(pinit,pxform); % taking log (pxform) of the initial par values (pinit)

% run fminsearch local optimizer using best guess from monte
% use at least 1000 maxiter when running a single experiment...
options  = optimset('maxiter',2500,'maxfunevals',10000,'Display','iter'); % set the options for your favorite optimizer
pbigbest = fminsearch(ofun,pvec0,options); % run your favorite optimizer
%pbigbest = fminsearch(ofun,pbigbest,options); % run your favorite optimizer again...
% save pbest
% load pbest.mat (if you want to skip parameter estimation)
pbest = pvec2struct(pbigbest,pxform)


%% plot the final results using the best-fit parameters
% that's it!
figure('Position',[239   558   990   420]);
if scale_par == 1
    xspan = 10.^[-2:0.25:6]'.*50.*pi; % Span of RTX concentration
else
    xspan = 10.^[-2:0.25:6]'; % Span of RTX concentration
end
plot_expt2(expt,pbest,xspan,'Xscale','log','Ylim',[-10 100],'Xtick',10.^[-2:2:6]);

%% create a summary table of results
par.Properties.RowNames = par.paramnames;
par.bestest = par.value;
fnames = fieldnames(pbest);
for j = 1:length(fnames)
    par{fnames{j},'bestest'} = pbest.(fnames{j});
end
disp(par);

if scale_par == 1
    writetable(par,'parameter_estimation_results_Kd_scaled.csv')
else
    writetable(par,'parameter_estimation_results_Kd.csv')
end

%% if you want to read in the results table and skip param est just skip the last three blocks and jump to here

if scale_par == 1
    parfin = readtable('parameter_estimation_results_Kd_scaled.csv');
else
    parfin = readtable('parameter_estimation_results_Kd.csv');
end

k = 0;
for j = 1:height(parfin)
    if parfin.fit(j) ~= 0
        k = k+1;
        pfinal.(parfin.paramnames{j}) = parfin.bestest(j);
    end
end

figure('Position',[239   558   990   420]);
xspan = 10.^[-2:0.25:6]';
plot_expt2(expt,pfinal,xspan,'Xscale','log','Ylim',[-10 100],'Xtick',10.^[-2:2:6]);

