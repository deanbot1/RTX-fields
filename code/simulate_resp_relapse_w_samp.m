%Out:[Tmat, Emat,CD20samps, CD16samps, CPXvec, LDHvec, perfvec]
%In: (CD20dist0, CD16dist0,tvec, T0, E0,adcxpars,CPXpars)
% Initialize some fake values for your function inputs
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

%%

figure;
plot(valuesCD20, pdCD20,'r-', 'LineWidth', 3)
set(gca,'FontSize',16, 'Xscale', 'log')
xlabel('number of CD20 receptors per tumor cell')
ylabel('pdf')
title('CD20 initial distribution on tumor cells')


figure;
plot(valuesCD16, pdCD16,'b-', 'LineWidth', 3)
set(gca,'FontSize',16, 'Xscale', 'log')
xlabel('number of CD16 receptors per NK cell')
ylabel('pdf')
xlim([5e2 10e3])
title('CD16 initial distribution on NK cells')



% %% This was only for generating simulated data from a distribution
% CD20dist0 = makedist('Normal','mu',100,'sigma',5); % estimate 100 receptors per T 
% CD16dist0 = makedist('Normal', 'mu',150', 'sigma', 5); % estimate 150 receptors per E 
% % From your data you will have x as values, y as frequencies
% valuesCD20 = 25:1:175;
% valuesCD16 = 100:1:200;
% pdCD20 = pdf(CD20dist0, valuesCD20);
% pdCD16 = pdf(CD16dist0, valuesCD16);
% 
% % Check that data is truly normalized pdf
% if sum(pdCD20)~=1
%     pdCD20 = pdCD20./sum(pdCD20);
% end
% if sum(pdCD16)~=1
%     pdCD16 = pdCD16./sum(pdCD16);
% end

%% Set your initial parameters arbitrarily
% time parameters
tf_et   = 100;
nr_t_et = 100;
tf_mol  = 80;
nr_t_mol = 80;

tvec_et  = linspace(0,tf_et,nr_t_et); % simulate 100 hours
tvec_mol = linspace(0,tf_mol,nr_t_mol); % simulate 100 hours

T0     = 1;     %1e6; 
E0toT0 = 2;     %0.5e6;
E0     = E0toT0*T0;
Estar0 = 1;

%ADXC parameters
% Doubling time from website on Boxnote: 20-30 hours (say, 25)
tdouble = 25;
g     = log(2)/25;
r     = 5;
kexp  = 1;
gamma = 0.5;
adcxpars  = [g,r,kexp,gamma];

% CPX pars
RTX  = 0.1;
% kd20 = 1;
% kd16 = 1;
%CD20
kon20  = 1;
koff20 = 0.2;
%CD16
kon16  = 0.5;
koff16 = 0.1;
gamma_perf = 0.5;
CPXpars = [RTX,kon20,koff20,kon16,koff16,gamma_perf];

%% Visualize your initial CD20 and CD16 distributions from literature
figure(1)
plot(valuesCD20, pdCD20, 'r', 'LineWidth',3)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('% Max MFI of T cell (\alpha CD20 receptors per T cell)')
set(gca,'Xscale','log')
title('Hiraga Fig 2B CD20 Histogram')

figure(2)
plot(valuesCD16, pdCD16, 'b','LineWidth',3)
set(gca,'FontSize',20,'LineWidth',1.5)
xlabel('% Max MFI of E cell (\alpha CD16 receptors per E cell)')
ylabel('Frequency (PDF)')
title('Srpan 2018 Fig 2A CD16 Histogram')


%% Sample from your distribution

% Scale CD20 distribution
nsamps = 10000;
[CD16samps]= randsmpl(pdCD16, nsamps, valuesCD16);
[CD20samps] = randsmpl(pdCD20, nsamps, valuesCD20);
% test the histcounts function using the values form the data

NCD160 = histcounts(CD16samps, valuesCD16);
NCD200 = histcounts(CD20samps, valuesCD20);
% % Tabulate the randomly drawn samples
% CD16samptbl = tabulate(CD16samps);
% CD16sampvals = CD16samptbl(:,1);
% CD16sampcts = CD16samptbl(:,2);
% CD16samppdf = CD16samptbl(:,3)/100;

% CD20samptbl = tabulate(CD20samps);
% CD20sampvals = CD20samptbl(:,1);
% CD20sampcts = CD20samptbl(:,2);
% CD20samppdf = CD20samptbl(:,3)/100;


%% Check that your sampled distributions match your simulated data

figure;
hold off
plot(valuesCD16(1:end-1), (NCD160./sum(NCD160)), 'c', 'LineWidth',3)
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
hold off
plot(valuesCD20(1:end-1),(NCD200./sum(NCD200)), 'm', 'LineWidth', 3)
hold on
plot(valuesCD20, pdCD20, 'r','LineWidth',3)
legend('Normalized sampled data', 'Hiraga Fig 2B data pdf')
legend boxoff
set(gca,'FontSize',20,'LineWidth',1.5, 'Xscale', 'log')
xlabel('% Max MFI of T cell (\alpha CD20 receptors per T cell)')
ylabel('Number of T cells')
title('Sample from literature CD20 histogram')

%% Run the loop that iterates through the CD16 and CD20 samples and runs
% the forward adcx and reaction_ss model

Tmat       = zeros(nr_t_et,nsamps);
Emat       = zeros(nr_t_et,nsamps);
LDHmat     = zeros(nr_t_et,nsamps);
perfmat    = zeros(nr_t_et,nsamps);
delCD16    = zeros(nsamps,1);
gamma = max(CD20samps)*1000;
lambda = max(CD16samps)*1e4;
for i = 1:nsamps
    i
%     [T,E,Estar,LDH,perf,CPX] = adcx(tf_mol,tf_et,nr_t_mol,nr_t_et,T0,E0,g,r,kexp,gamma,...
%     CD20,CD16,RTX,kon20,koff20,kon16,koff16,gamma_perf);
    % generate one column of your sampled data by calling the adcx model
    % which calles the reaction_ss function
    [Ti, Ei, Estari, LDHi, perfi, CPXi] = adcx(tf_mol,tf_et,nr_t_mol,nr_t_et,...
        T0,E0toT0,Estar0,g,r,kexp,gamma,...
        (CD20samps(i)/gamma),(CD16samps(i)/lambda),RTX,kon20,koff20,kon16,koff16,gamma_perf);
    i
    Tmat(:,i)     = Ti;
    Emat(:,i)     = Ei;
    LDHmat(:,i)   = LDHi;
    perfmat(:,i)  = perfi;
    delCD16(i,1)  = CPXi;
    
end

% For a given CD16, CD20 distribution, here are the expected LDH and perf
% trajectories over time
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
xlim([1 100])
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

