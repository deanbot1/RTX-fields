%Out:[Tmat, Emat,CD20samps, CD16samps, CPXvec, LDHvec, perfvec]
%In: (CD20dist0, CD16dist0,tvec, T0, E0,adcxpars,CPXpars)
% Initialize some fake values for your function inputs
close all; clear all; clc;

% Inport data and 
tableCD20 = csvread('../data/Hiraga_2B1_CD20_BEFORE_RTX.csv');
rawvaluesCD20 = tableCD20(:,1);
valuesCD20= rawvaluesCD20./max(rawvaluesCD20); % step to normalize values
freqCD20 = tableCD20(:,2);
pdCD20 = freqCD20./sum(freqCD20);

tableCD16 = csvread('../data/Srpan_2A_CD16_0.csv');
valuesCD16 = tableCD16(:,1); % normalized values from 0 to 1
freqCD16 = tableCD16(:,2);
% remove negative frequencies
ireal = freqCD16>=0;
pdCD16 = freqCD16(ireal)./sum(freqCD16(ireal));
valuesCD16 = valuesCD16(ireal);




%% This was only for generating simulated data from a distribution
CD20dist0 = makedist('Normal','mu',100,'sigma',5); % estimate 100 receptors per T 
CD16dist0 = makedist('Normal', 'mu',150', 'sigma', 5); % estimate 150 receptors per E 
% From your data you will have x as values, y as frequencies
valuesCD20 = 25:1:175;
valuesCD16 = 100:1:200;
pdCD20 = pdf(CD20dist0, valuesCD20);
pdCD16 = pdf(CD16dist0, valuesCD16);

% Check that data is truly normalized pdf
if sum(pdCD20)~=1
    pdCD20 = pdCD20./sum(pdCD20);
end
if sum(pdCD16)~=1
    pdCD16 = pdCD16./sum(pdCD16);
end
%% Set your initial parameters arbitrarily

tvec = 0:1:100; % simulate 100 hours

T0 = 1e6; 
E0 = 0.5e6;
%ADXCpars
g=0.01;
r=5;
kexp =1;
adcxpars = [g, r, kexp];

%CPXpars
R = 0.1;
kd20=1;
kd16=1;
CPXpars = [R, kd20, kd16];


%% Visualize your initial CD20 and CD16 distributions from literature
figure;
plot(valuesCD20, pdCD20, 'r', 'LineWidth',3)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('% Max MFI of T cell (\alpha CD20 receptors per T cell)')
title('Hiraga Fig 2B CD20 Histogram')

figure;
plot(valuesCD16, pdCD16, 'b','LineWidth',3)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('% Max MFI of E cell (\alpha CD16 receptors per E cell)')
ylabel('Frequency (PDF)')
title('Srpan 2018 Fig 2A CD16 Histogram')
%% Sample from your distribution

nsamps = 10000;
[CD16samps]= randsmpl(pdCD16, nsamps, valuesCD16);
[CD20samps] = randsmpl(pdCD20, nsamps, valuesCD20);
% Tabulate the randomly drawn samples
CD16samptbl = tabulate(CD16samps);
CD16sampvals = CD16samptbl(:,1);
CD16sampcts = CD16samptbl(:,2);
CD16samppdf = CD16samptbl(:,3)/100;

CD20samptbl = tabulate(CD20samps);
CD20sampvals = CD20samptbl(:,1);
CD20sampcts = CD20samptbl(:,2);
CD20samppdf = CD20samptbl(:,3)/100;


%% Check that your sampled distributions match your simulated data

figure;
plot(CD16sampvals, CD16samppdf, 'c', 'LineWidth',3)
hold on
plot(valuesCD16, pdCD16, 'b','LineWidth',3)
set(gca,'FontSize',20,'LineWidth',1.5)
legend('Normalized sampled data', 'Srpan 2018 Fig 2A data pdf')
legend boxoff
xlabel('% Max MFI of E cell (\alpha CD16 receptors per E cell)')
ylabel('Frequency (pdf)')
title('Sample from literature CD16 histogram')
%%
figure;
plot(CD20sampvals, CD20samppdf, 'm', 'LineWidth', 3)
hold on
plot(valuesCD20, pdCD20, 'r','LineWidth',3)
legend('Normalized sampled data', 'Hiraga Fig 2B data pdf')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('% Max MFI of T cell (\alpha CD20 receptors per T cell)')
ylabel('Number of T cells')
title('Sample from literature CD20 histogram')
%% Run the loop that iterates through the CD16 and CD20 samples and runs
% the forward adcx and reaction_ss model

for i = 1:nsamps
    % generate one column of your sampled data by calling the adcx model
    % which calles the reaction_ss function
    [Ti, Ei, Estari, LDHi, perfi, CPXi] = adcx(tvec, T0,E0, adcxpars,...
        CD20samps(i), CD16samps(i),CPXpars);
    Tmat(:,i) = Ti;
    Emat(:,i) = Estari;
    LDHmat(:,i) = LDHi;
    perfmat(:,i) =perfi;
    delCD16(i,1) = CPXi;
    
end
% For a given CD16, CD20 distribution, here are the expected LDH and perf
% trajectories over time
LDHvec = sum(LDHmat,2);
perfvec = sum(perfmat,2);



