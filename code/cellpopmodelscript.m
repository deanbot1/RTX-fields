%Out:[Tmat, Emat,CD20samps, CD16samps, CPXvec, LDHvec, perfvec]
%In: (CD20dist0, CD16dist0,tvec, T0, E0,adcxpars,CPXpars)
% Initialize some fake values for your function inputs
close all; clear all; clc;

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


%% Visualize your initial CD20 and CD16 distributions
figure;
plot(valuesCD20, pdCD20, 'r', 'LineWidth',3)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('CD20 receptors per T cell')
ylabel('PDF of CD20 levels')

figure;
plot(valuesCD16, pdCD16, 'b','LineWidth',3)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('CD16 receptors per E cell')
ylabel('PDF of CD16 levels')
%% Sample from your distribution

nsamps = 1000;
[CD16samps]= randsmpl(pdCD16, nsamps, valuesCD16);
[CD20samps] = randsmpl(pdCD20, nsamps, valuesCD20);

%% Check that your sampled distributions match your simulated data
figure;
histogram(CD16samps, 'Normalization', 'probability')
hold on
plot(valuesCD16, pdCD16, 'r', 'LineWidth',3)
set(gca,'FontSize',20,'LineWidth',1.5)
legend('Normalized sampled data', 'simulated data pdf')
legend boxoff
xlabel('CD16 receptors per E cell')
ylabel('Number of E cells (E or E*)')
title('Sample from literature CD16 histogram')

figure;
histogram(CD20samps, 'Normalization', 'probability')
hold on
plot(valuesCD20, pdCD20, 'b','LineWidth',3)
legend('Normalized sampled data', 'simulated data pdf')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('CD20 receptors per T cell')
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



