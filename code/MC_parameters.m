% MC scheme for bootstrapped parameter sets
close all; clear all; clc
paramstruct = load('bootstrap_500runs.mat');
pbest = paramstruct.pbest;
pbigboot = paramstruct.pbigbootall;
ofun = paramstruct.ofun; % This returns the fval for a parameter set
pinit = paramstruct.pinit;
% residuals = paramstruct.residuals;
% sigma = std(residuals);
%% Use plotmatrix to visualize bootstrapped parameters
nparams = size(pbigboot,1);
paramnames = fieldnames(pinit); 
figure;
[S,AX,BigAx,H,HAx] = plotmatrix(pbigboot');
title(BigAx,'A Comparison of Bootstrapped Parameter Estimates ')
for i = 1:nparams
ylabel(AX(i,1),paramnames{i})
xlabel(AX(7,i),paramnames{i})
end
%% Monte Carlo
nruns = 500;
% Grab the 99th percentiles from bootstrapping parameter ranges to sample
% from
phatbounds =prctile(pbigboot',[0.5,99.5], 1);
nparams = length(pbest);

% Uniform sampling from k-dimensional hyper-rectangle
pxform = paramstruct.pxform;
randmat = rand(nparams, nruns);
pbig = pstruct2vec(pbest,pxform);
fvalbest = ofun(pbig);
pbestMC = pbig;
%% Run forward model from uniform sampling from parameter bounds 
for i = 1:nruns
    % get your parameters
    phati = randmat(:,i).*(phatbounds(2,:)-phatbounds(1,:)) + phatbounds(1,:);
    % store parameters
    pbigMC(:,i) =phati';
    % compute function evaluation and save
    fval(i) =ofun(phati');
    if fval(i)< fvalbest
        fvalbest = fval(i);
        pbestMC = phati';
    end
end
% Not sure why but this returns "matrix is badly scaled, results may be
% inaccurate"
%% Find the parameters that fall within certain fval ranges 
DO = [0.68, 0.90, 0.95]; 

% From your distribution of fvals, find the 68, 90, and 95th lowest fvals
for i= 1:length(DO)
    ikeep = [];
    % fbest + deltaf from chi-squared table with nparams df and D0 % CI
    ikeep = fval<= fvalbest + chi2inv(DO(i), nparams);

    % save those parameter sets
    paramsMC{i}=pbigMC(:,ikeep);
end

%% Plot using plotmatrix
% Not sure how to plot on the same matrix
figure;
for i = 1:length(DO)
[S,AX,BigAx,H,HAx] = plotmatrix(paramsMC{i}');
hold on

title(BigAx,'A Comparison of Parameter Estimates ')
for j = 1:nparams
ylabel(AX(j,j),paramnames{j})
xlabel(AX(7,j),paramnames{j})
end
end
save('paramsMC', 'paramsMC')