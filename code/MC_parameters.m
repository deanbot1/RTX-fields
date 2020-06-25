% MC scheme for bootstrapped parameter sets
close all; clear all; clc
paramstruct = load('bootstrap_500runs.mat');
pbest = paramstruct.pbest;
pbigboot = paramstruct.pbigbootall;
ofun = paramstruct.ofun; % This returns the fval for a parameter set
pinit = paramstruct.pinit;
residuals = paramstruct.residuals;
sigma = std(residuals);
nparams = size(pbigboot,1);
paramnames = fieldnames(pinit); 
%% Use plotmatrix to visualize bootstrapped parameters

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
fvalbest = ofun(pbig)./(sigma.^2);
pbestMC = pbig;
%% Run forward model from uniform sampling from parameter bounds 
for i = 1:nruns
    % get your parameters
    phati = randmat(:,i).*(phatbounds(2,:)-phatbounds(1,:)) + phatbounds(1,:);
  
    % store parameters
    pbigMC(:,i) =phati';
    % compute function evaluation and save
    fval(i) =ofun(phati')/(sigma.^2);
    if fval(i)< fvalbest
        fvalbest = fval(i);
        pbestMC = phati';
    end
end
% Not sure why but this returns "matrix is badly scaled, results may be
% inaccurate"
%% Compare fvals to fvalbest
DO = [0.68, 0.90, 0.95]; 
figure;
histogram(fval, 100)
hold on
plot(fvalbest, 1, '*')
for i = 1:length(DO)
plot(fvalbest + chi2inv(DO(i), nparams), '*')
end
legend('fvals', 'fvalbest', '68th','90th', '95th')
xlabel('fval')
ylabel('frequency')
%% Find the parameters that fall within certain fval ranges 
DO = [0.68, 0.90, 0.95]; 

% From your distribution of fvals, find the 68, 90, and 95th lowest fvals
for i= 1:length(DO)
    ikeep = [];
    % fbest + deltaf from chi-squared table with nparams df and D0 % CI
    for j = 1:length(fval)
        if fval(j)<= fvalbest + chi2inv(DO(i), nparams)
            indkeep=j; % record the row
            ikeep = vertcat(ikeep, indkeep); % save the row
        end
    end

    % save those parameter sets
    paramsMC{i}=pbigMC(:,ikeep);
end
%% Plot overlapping histograms, also need to plot parameter relationships...
figure;
for i = 1:7
 subplot(1,7,i)
for j = 1:2length(DO)
    k = length(DO)+1-j;
    parammat = paramsMC{k};
    histogram(parammat(i,:),linspace(phatbounds(1,i), phatbounds(2,i), 20))
    hold on
    xlabel([paramnames{i}])
    legend('95th', '90th', '68th')
    set(gca,'FontSize',16,'LineWidth',1.5)
end
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