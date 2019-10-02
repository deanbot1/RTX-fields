 % Combine the parameter estimates and the full distribution of CD20 and
 % CD16 levels to run the full model and observe the effects on the tumor
 % dynamics

%
close all; clear all; clc;

% Import data and 
tableCD20 = csvread('../data/Hiraga_2B1_CD20_BEFORE_RTX.csv');
valuesCD20 = tableCD20(:,1);
% valuesCD20 = rawvaluesCD20./max(rawvaluesCD20); % step to normalize values
freqCD20 = tableCD20(:,2);
pdCD20   = freqCD20./sum(freqCD20);

tableCD16 = csvread('../data/Srpan_2A_CD16_0.csv');
valuesCD16 = tableCD16(:,1); % values x axis
freqCD16 = tableCD16(:,2);   % remove negative frequencies
ireal = freqCD16>=0;
pdCD16 = freqCD16(ireal)./sum(freqCD16(ireal));
valuesCD16 = valuesCD16(ireal);

%% Plot the distributions of CD20 and CD16 expression levels from literature

figure;
plot(valuesCD20, pdCD20,'r-', 'LineWidth', 3)
set(gca,'FontSize',16, 'Xscale', 'log')
xlabel('number of CD20 receptors per tumor cell')
ylabel('pdf')
title('CD20 initial distribution on tumor cells (Srpan 2016 Fig2A)')


figure;
plot(valuesCD16, pdCD16,'b-', 'LineWidth', 3)
set(gca,'FontSize',16, 'Xscale', 'log')
xlabel('number of CD16 receptors per NK cell')
ylabel('pdf')
xlim([5e2 10e3])
title('CD16 initial distribution on NK cells (Hiraga Fig 2B)')



%% Set the parameters from the parameter estimation table

param_table =readtable('parameter_estimation_results.csv');

% Loop through the table and set all of the paramnames to the values as
% variables 
pnames = param_table.paramnames;
% Make a new structure 
pstruct = struct('cell', 0);
% Assign the fields to the parameter names, and the values to the best
% estimates
for i = 1:length(pnames)
    pstruct = setfield(pstruct,char(pnames(i)), param_table.bestest(i));

end



%% Sample from your distribution

% Scale CD20 distribution
nsamps = 100;
[CD16samps]= randsmpl(pdCD16, nsamps, valuesCD16);
[CD20samps] = randsmpl(pdCD20, nsamps, valuesCD20);
% test the histcounts function using the values form the data

NCD160 = histcounts(CD16samps, valuesCD16);
NCD200 = histcounts(CD20samps, valuesCD20);



%% Run the loop that iterates through the CD16 and CD20 samples and runs
% the forward adcx and reaction_ss model

lambda = max(CD16samps)*1e4;



% Call cellpopmodel to obtain the individual tumor and effector cell
% trajectories
expt = 'Z138';
[Tmat, Emat, CD20samps, CD16samps, delCD16, LDHmat, perfmat] = cellpopmodel(CD20samps,...
    CD16samps,nsamps, lambda, pstruct, expt);
LDHvec = sum(LDHmat,2);
perfvec = sum(perfmat,2);


%% Plotting
figure(5)
hold off
plot(sum(Tmat,2),'LineWidth',2)
hold on
plot(sum(Emat,2),'LineWidth',2)
plot(LDHvec,'LineWidth',2)
xlabel('Time (hours)')
%xlim([1 100])
ylabel('Number of cells/molecules')
title('Short term dynamics')
legend('Target cells','Effector cells (CD16)','LDH')
legend boxoff
set(gca,'FontSize',16)

figure(6)
hold off
plot(sum(Tmat,2),'LineWidth',2)
hold on
plot(sum(Emat,2),'LineWidth',2)
plot(LDHvec,'LineWidth',2)
xlabel('Time (hours)')
ylabel('Number of cells/molecules')
title('Long Range Dynamics')
legend('Target cells','Effector cells (CD16)','LDH')
legend boxoff
set(gca,'FontSize',16)
Tvec = sum(Tmat,2)
i120 = find(Tvec>1.2*Tvec(1),1, 'first')

