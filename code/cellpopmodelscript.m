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



%% Sample from your literature distribution to make sure sampling works well

% Scale CD20 distribution
nsamps = 10000;
[CD16samps]= randsmpl(pdCD16, nsamps, valuesCD16);
[CD20samps] = randsmpl(pdCD20, nsamps, valuesCD20);
% test the histcounts function using the values form the data

NCD160 = histcounts(CD16samps, valuesCD16);
NCD200 = histcounts(CD20samps, valuesCD20);

figure;
plot(valuesCD16(1:end-1), NCD160, 'b.')
xlabel('CD16 expression')
ylabel('Number of cells')
xlim([5e2 10e3])
title('Test Sampling CD16')
set(gca,'FontSize',16, 'Xscale', 'log')


figure;
plot(valuesCD20(1:end-1), NCD200, 'r.')
xlabel('CD20 expression')
ylabel('Number of cells')
title('Test Sampling CD20')
set(gca,'FontSize',16, 'Xscale', 'log')

logmeanCD16 = mean(log(CD16samps));
logmeanCD20 = mean(log(CD20samps));

%% Rescale the distribution so that the new mean is equal to the parameter estimation value

expt = 'Z138';
% Scale them by the log of the mean (geometric mean) which is the
% equivalent of dividing by the log of the mean... Try this and see what
% happens
%MAKE THIS CLEANER 
switch expt
    case 'Z138'
    goallogmeanCD20 = log(pstruct.CD20_Z138);
    goalmeanCD20 = pstruct.CD20_Z138;
    case 'SUDHL4'
    goallogmeanCD20 = log(pstruct.CD20_SUDHL4);
    goalmeanCD20 = pstruct.CD20_SUDHL4;
end  
goalmeanCD16 = pstruct.CD16;
goallogmeanCD16 = log(pstruct.CD16);

%% Rescale distributions linearly

currmeanCD20samps = mean(CD20samps);
currmeanCD16samps = mean(CD16samps);
 
CD20new = (CD20samps./mean(CD20samps)).*goalmeanCD20;
valuesCD20new = (valuesCD20./mean(CD20samps)).*goalmeanCD20;
CD20test = mean(CD20new);
NCD20new = histcounts(CD20new, valuesCD20new);

CD16new = (CD16samps./mean(CD16samps)).*goalmeanCD16;
valuesCD16new = (valuesCD16./mean(CD16samps)).*goalmeanCD16;
CD16test = mean(CD16new);
NCD16new = histcounts(CD16new, valuesCD16new);


figure;
plot(valuesCD20new(1:end-1), NCD20new, 'r.')
hold on
plot(CD20test, 1,'r*','LineWidth', 3)
xlabel('CD20 expression')
ylabel('Number of cells')
legend('samples', 'mean', 'Location', 'NorthWest')
legend boxoff
title('Rescaled CD20 distribution linearly')
set(gca,'FontSize',16, 'Xscale', 'log')


figure;
plot(valuesCD16new(1:end-1), NCD16new, 'b.')
hold on
plot(CD16test, 1,'b*','LineWidth', 3)
xlabel('CD16 expression')
ylabel('Number of cells')
xlim([0.5 20])
legend('samples', 'mean', 'Location', 'NorthWest')
legend boxoff
title('Rescaled CD16 distribution linearly')
set(gca,'FontSize',16, 'Xscale', 'log')
%% Rescale distributions logarithmically-- doesn't work because these 
% become negative

CD16lognew = log(CD16samps) + log(goalmeanCD16./currmeanCD16samps);
valuesCD16lognew = log(valuesCD16) + log(goalmeanCD16./currmeanCD16samps);
CD16logtest = mean(CD16lognew)
NCD16lognew = histcounts(CD16lognew, valuesCD16lognew);

CD20lognew = log(CD20samps) + log(goalmeanCD20./currmeanCD20samps);
valuesCD20lognew = log(valuesCD20) + log(goalmeanCD20./currmeanCD20samps);
CD20logtest = mean(CD20lognew)
NCD20lognew = histcounts(CD20lognew, valuesCD20lognew);

figure;
plot(valuesCD16lognew(1:end-1), NCD16lognew, 'b.')
hold on
plot(CD16logtest, 1,'b*','LineWidth', 3)
xlabel('CD16 expression')
ylabel('Number of cells')
title('Rescaled CD16 distribution logarithmically')
set(gca,'FontSize',16, 'Xscale', 'log')

figure;
plot(valuesCD20lognew(1:end-1), NCD20lognew, 'r.')
hold on
plot(CD20logtest, 1,'r*','LineWidth', 3)
xlabel('CD20 expression')
ylabel('Number of cells')
title('Rescaled CD20 distribution logarithmically')
set(gca,'FontSize',16, 'Xscale', 'log')


%% Run the loop that iterates through the CD16 and CD20 samples and runs
% the forward adcx and reaction_ss model

% I Can't remember what Lambda even does :(
lambda = 1e-2;

% Call cellpopmodel to obtain the individual tumor and effector cell
% trajectories
expt = 'Z138';
%% UNIFORM Simulate using mean only: Just set all the samples to be the parameter value
nsamps = 1000
CD20uniform = exp(goalmeanCD20)*ones(nsamps,1);
CD16uniform = exp(goalmeanCD16)*ones(nsamps,1);
CD20distrib = CD20uniform;
CD16distrib = CD16uniform;
%% Use sampled distributions as input to the model

CD20distrib = CD20new;
CD16distrib = CD16new;
%% Run population level model using CD16 and CD20 value from input each time

% Output of this is rows=time, columns = samples. We sum all of the samples
% to get the total population level behavior
[Tmat, Emat, CD20mat, CD16mat, delCD16, LDHmat, perfmat] = cellpopmodel(CD20distrib,...
    CD16distrib,nsamps, lambda, pstruct, expt);
LDHvec = sum(LDHmat,2);
perfvec = sum(perfmat,2);
tvec =  linspace(0,pstruct.tf_et,pstruct.nr_t_et);

figure;
hold off
plot(tvec,sum(Tmat,2),'LineWidth',2)
hold on
plot(tvec, sum(Emat,2),'LineWidth',2)
plot(tvec, LDHvec,'LineWidth',2)
xlabel('Time')
ylabel('Number of cells/molecules')
title('Tumor and Effector Cell Dynamics with mean CD16 and CD20')
legend('Target cells','Effector cells (CD16)','LDH')
legend boxoff
set(gca,'FontSize',16)




%% Use the output CD20mat and CD16 mat to make videos of distributions 
% Here we want to plot the distribution along all the columns, at each of
% the time points given by the row
% We want to plat the average CD20 and CD16 level per cell over time
% This should show that the average levels decrease over time
Tvec = sum(Tmat,2);
Evec = sum(Emat,2);
%Estarvec = sum(Estarmat,2);
CD20vec = sum(CD20mat,2); % total CD20 levels over time
CD16vec = sum(CD16mat,2);


figure;
subplot(1,2,1)
plot(tvec, CD20vec./Tvec,'g','LineWidth', 2)
hold on
xlabel('Time (hours)')
ylabel('Average number of CD20 receptors per cell')
title('Mean CD20 receptors per tumor cell during RTX treatment')
set(gca,'FontSize',16)
%xlim([0 nr_t_et])
subplot(1,2,2)
plot(tvec, Tvec,'b','LineWidth', 2)
hold on
xlabel('Time ')
ylabel('Number of tumor cells')
title('Tumor cells in time with constant RTX')
set(gca,'FontSize',16)
%xlim([0 nr_t_et])

figure;
subplot(1,2,1)
plot(tvec, CD16vec./Evec,'c','LineWidth', 2)
hold on
xlabel('Time')
ylabel('Average number of CD16 receptors per NK cell')
title('Mean CD16 receptors per NK cell during RTX treatment')
set(gca,'FontSize',16)
%xlim([0 nr_t_et])
subplot(1,2,2)
plot(Evec,'b','LineWidth', 2)
hold on
%plot(Estarvec, 'k', 'LineWidth',2)
xlabel('Time (hours)')
%legend('active NK cells', 'depleted NK cells')
ylabel('Number of NK cells')
title('NK cells under constant RTX')
set(gca,'FontSize',16)
%xlim([0 nr_t_et])

%% Plot the total CD20 receptor distribution over time and the pdf
% CD20 distributions
% Make this into a video
figure;
for i = 1:pstruct.nr_t_et
subplot(1,2,1)
hold off
plot(valuesCD20new(1:end-1), NCD20new,'b.', 'LineWidth', 2)
hold on
% at each row (time point) bin the CD20 levels for the individual combos
NCD20i = histcounts(CD20mat(i,:), valuesCD20new);
plot(valuesCD20new(1:end-1), NCD20i, 'c.','LineWidth', 0.5)
legend('t=0', ['t=', num2str(i), 'hr'], 'Location', 'NorthWest')
legend boxoff
ylim([0 90])
%xlim([])
set(gca,'FontSize',16,'LineWidth',1.5,'Xscale', 'log')
ylabel('Number of CD20 receptors')
xlabel('\lambda*number of CD20 receptors per cell')
title(['Total CD20 distribution on tumor cells, t=', num2str(i), 'hours'])

subplot(1,2,2)
hold off
plot(valuesCD16new(1:end-1), NCD16new, 'r.', 'LineWidth', 2)
%plot(CD20edges(1:end-1), pCD20(1,:),'r', 'LineWidth', 2)
hold on
NCD16i = histcounts(CD16mat(i,:), valuesCD16new);
plot(valuesCD16new(1:end-1), NCD16i, 'm.','LineWidth', 0.5)
legend('t=0', ['t=', num2str(i), 'hr'], 'Location', 'NorthWest')
legend boxoff

set(gca,'FontSize',16,'LineWidth',1.5,'Xscale', 'log')
ylabel('Number of CD16 receptors')
xlabel('\lambda*number of CD16 receptors per cell')
title(['Total CD16 distribution on NK cells, t=', num2str(i), 'hours'])
xlim([ 0.1 10])
ylim([ 0 160])
drawnow
pause(0.01)
end


%% Plot the tumor cell trajectories for high and low CD16 and CD20 values

figure;
for i = 1:100
    if CD20new(i)>median(CD20new)
    plot(Tmat(:,i), 'r', 'LineWidth', 1)
    end
    if CD20new(i)<median(CD20new)
    plot(Tmat(:,i), 'b', 'LineWidth', 1)
    end
    hold on
    %legend('CD20- on tumor cell', 'CD20+ on tumor cell')
    xlabel('Time(hours)')
    ylabel('Relative number of tumor cells')
    title('Observed variability in tumor cell trajectories')
    set(gca,'FontSize',16)
    %ylim([0 14])
end

figure;
for i = 2:101
    if CD16new(i)>median(CD16new) 
    plot(Tmat(:,i), 'g', 'LineWidth', 1)
    end
    if CD16new(i)<median(CD16new)
    plot(Tmat(:,i), 'y', 'LineWidth', 1)
    end
    hold on
    %legend('CD16 + on NK Cell', 'CD16 - on NK cell')
    xlabel('Time(hours)')
    ylabel('Relative number of tumor cells')
    title('Observed variability in tumor cell trajectories')
    set(gca,'FontSize',16)
    %ylim([0 14])
end
