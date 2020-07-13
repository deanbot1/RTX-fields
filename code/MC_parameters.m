% Monte Carlo sampling of parameter space


close all; clear all; clc

paramstruct = load('bootstrap_500runs_July4.mat');
pbest = paramstruct.pbest;
pbigboot = paramstruct.pbigbootall;
ofun = paramstruct.ofun; % This returns the fval for a parameter set
pinit = paramstruct.pinit;
residuals = paramstruct.residuals;
sigma = std(residuals);
nparams = size(pbigboot,1);
paramnames = fieldnames(pinit); 
pbestnames = fieldnames(pbest);




% Convert pbest to a matrix 
pbesttab = struct2table(pbest);
pbestmat = pbesttab{:,:};






%% Define parameter range
nruns = 500;
% Grab the 99th percentiles from bootstrapping parameter ranges to sample
% from
bootbounds =prctile(pbigboot',[0.5,99.5], 1);
% Use the spread in the bootstrap parameter estiamtes to sample
% symmetrically
diamvec = bootbounds(2,:) - bootbounds(1,:);
diamvec = diamvec';

% Uniform sampling from k-dimensional hyper-rectangle
pxform = paramstruct.pxform;
randmat = rand(nparams, nruns);
pbig = pstruct2vec(pbest,pxform);
fvalbest = ofun(pbig)./(sigma.^2);
pbestMC = pbig;
factor = 0.4;
phatbounds = horzcat((pbestMC-factor*diamvec), (pbestMC+factor*diamvec));
phatbounds = phatbounds';

%% Use plotmatrix to visualize bootstrapped parameters

figure;
[S,AX,BigAx,H,HAx] = plotmatrix(pbigboot');
for k = 1:nparams
    hold(HAx(1,k),'on')
plot(HAx(1,k),pbestMC(k),1,'g*','LineWidth',2)
for m = 1:nparams
    hold(AX(m,k),'on')
    if m~=k
    plot(AX(m,k), pbestMC(k), pbestMC(m), 'r*', 'LineWidth',2)
    end
end
end
% S are the scatter plots, H are the histograms
title(BigAx,'A Comparison of Bootstrapped Parameter Estimates ')
for j = 1:nparams
   
ylabel(AX(j,1),paramnames{j})
xlabel(AX(7,j),paramnames{j})
end
%% Run forward model from uniform sampling from parameter bounds 
fval = [];

for i = 1:nruns
    w = warning('on', 'all');
    id = w.identifier;
    warning('off',id);
    % get your parameters
    phati = randmat(:,i)'.*(phatbounds(2,:)-phatbounds(1,:)) + phatbounds(1,:);
  
    % store parameters
    pbigMC(:,i) =phati';
    % compute function evaluation and save
    ofunval = ofun(phati')./(sigma.^2);
    if isreal(ofunval)
        fval = vertcat(fval,ofunval);
    end
   
    if ofunval< fvalbest
        fvalbest = ofunval;
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
plot([(fvalbest + chi2inv(DO(i), nparams)) (fvalbest + chi2inv(DO(i), nparams))],[0 300], '--','LineWidth',2)
end
xlim([0 50])
legend('fvals', 'fvalbest', '68th','90th', '95th')
legend boxoff
xlabel('fval')
ylabel('frequency')
set(gca,'FontSize',14,'LineWidth',1.5)
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

%% Plot using plotmatrix
% Not sure how to plot on the same matrix
figure;

% first plot the 95th percentile using plot matrix
[S,AX,BigAx,H,HAx] = plotmatrix(paramsMC{3}');
for j = 1:nparams
ylabel(AX(j,1),paramnames{j})
xlabel(AX(7,j),paramnames{j})
end


for k = 1:nparams
    for i = 1:2
    hold(HAx(1,k),'on')
    % Add the parameter matrices for each percentile
    pmat = paramsMC{i};
    histogram(HAx(1,k), pmat(k,:))
    
    % Add the scatter plot for each percentile
    for m = 1:nparams
    hold(AX(m,k),'on')
    if m~=k
    plot(AX(m,k), pmat(k,:), pmat(m,:), '.')
    end
    end
    end
end

% Show were the best parameter value falls
hold on
    for k = 1:nparams
    hold(HAx(1,k),'on')
    plot(HAx(1,k),pbestMC(k),1,'g*','LineWidth',2)
   for m = 1:nparams
    hold(AX(m,k),'on')
    if m~=k
    plot(AX(m,k), pbestMC(k), pbestMC(m), 'r*', 'LineWidth',2)
    end
    end
    end
title(BigAx,'Monte Carlo Parameter Estimates')
for j = 1:nparams
ylabel(AX(j,j),paramnames{j})
xlabel(AX(7,j),paramnames{j})
end
save('paramsMC', 'paramsMC')