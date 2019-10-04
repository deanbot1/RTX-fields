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
nsamps = 1000;
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

currmeanCD16 = mean(CD16samps)
currmeanCD20 = mean(CD20samps)

%% Rescale the distribution so that the new mean is equal to the parameter estimation value
expt = 'Z138';

% Honestly, not exactly sure how to do this... 
switch expt
    case 'Z138'
    goalmeanCD20 = pstruct.CD20_Z138;
    case 'SUDHL4'
    goalmeanCD20 = pstruct.CD20_SUDHL4;
end  
goalmeanCD16 = pstruct.CD16;
currminCD16 = min(CD16samps);
currmaxCD16 = max(CD16samps);
currminCD20 = min(CD20samps);
currmaxCD20 = max(CD20samps);
% First try just shifting the distributions

CD20samps_shift = CD20samps-(currmeanCD20-goalmeanCD20);
valuesCD20_shift = valuesCD20-(currmeanCD20-goalmeanCD20);
checkCD20mean = mean(CD20samps_shift);
CD16samps_shift = CD16samps-(currmeanCD16-goalmeanCD16);
valuesCD16_shift = valuesCD16-(currmeanCD16-goalmeanCD16);
checkCD16mean = mean(CD16samps_shift);

NCD16_shift = histcounts(CD16samps_shift, valuesCD16_shift);
NCD20_shift = histcounts(CD20samps_shift, valuesCD20_shift);

% Plot the shifted, but not scaled, distributions
figure;
plot(valuesCD16_shift(1:end-1), NCD16_shift, 'b.')
hold on
plot(checkCD16mean, 1, 'b*')
xlabel('CD16 expression')
ylabel('Number of cells')
%xlim([5e2 10e3])
title('Shifted CD16 distribution')
%set(gca,'FontSize',16, 'Xscale', 'log')


figure;
plot(valuesCD20_shift(1:end-1), NCD20_shift, 'r.')
hold on
plot(checkCD20mean, 1, 'r*')
xlabel('CD20 expression')
ylabel('Number of cells')
title('Shifted CD20 distribution')
%set(gca,'FontSize',16, 'Xscale', 'log')

% Next we need to scale these but I need some feedback on how to do this...


%% Run the loop that iterates through the CD16 and CD20 samples and runs
% the forward adcx and reaction_ss model

% I Can't remember what Lambda even does :(
lambda = max(CD16samps_shift)*1e4;

% Call cellpopmodel to obtain the individual tumor and effector cell
% trajectories
expt = 'Z138';
% Simulate using mean only: Just set all the samples to be the parameter value
nsamps = 1000
CD20uniform = goalmeanCD20*ones(nsamps,1);
CD16uniform = goalmeanCD16*ones(nsamps,1);

[Tmatmean, Ematmean, CD20matunif, CD16matunif, delCD16, LDHmat, perfmat] = cellpopmodel(CD20uniform,...
    CD16uniform,nsamps, lambda, pstruct, expt);
LDHvec = sum(LDHmat,2);
perfvec = sum(perfmat,2);

figure;
hold off
plot(sum(Tmatmean,2),'LineWidth',2)
hold on
plot(sum(Ematmean,2),'LineWidth',2)
plot(LDHvec,'LineWidth',2)
xlabel('Time (hours)')
ylabel('Number of cells/molecules')
title('Tumor and Effector Cell Dynamics with mean CD16 and CD20')
legend('Target cells','Effector cells (CD16)','LDH')
legend boxoff
set(gca,'FontSize',16)



%% Simulate using the shifted distribution
% Unfortunately, this doesn't seem to be working well (which isn't
% surprising because the samples of CD20 and CD16 are negative which
% doesn't make any sense
[Tmat, Emat, CD20mat, CD16mat, delCD16, LDHmat, perfmat] = cellpopmodel(CD20samps_shift,...
    CD16samps_shift,nsamps, lambda, pstruct, expt);
LDHvec = sum(LDHmat,2);
perfvec = sum(perfmat,2);

figure;
hold off
plot(sum(Tmat,2),'LineWidth',2)
hold on
plot(sum(Emat,2),'LineWidth',2)
plot(LDHvec,'LineWidth',2)
xlabel('Time (hours)')
ylabel('Number of cells/molecules')
title('Tumor and Effector Cell Dynamics with shifted distribution of CD16 and CD20')
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
Estarvec = sum(Estarmat,2);
CD20vec = sum(CD20mat,2); % total CD20 levels over time
CD16vec = sum(CD16mat,2);


figure;
subplot(1,2,1)
plot(CD20vec./Tvec,'g','LineWidth', 2)
hold on
xlabel('Time (hours)')
ylabel('Average number of CD20 receptors per cell')
title('Mean CD20 receptors per tumor cell during RTX treatment')
set(gca,'FontSize',16)
%xlim([0 nr_t_et])
subplot(1,2,2)
plot(Tvec,'b','LineWidth', 2)
hold on
xlabel('Time (hours)')
ylabel('Number of tumor cells')
title('Tumor cells in time with constant RTX')
set(gca,'FontSize',16)
%xlim([0 nr_t_et])

figure;
subplot(1,2,1)
plot(CD16vec./Evec,'c','LineWidth', 2)
hold on
xlabel('Time (hours)')
ylabel('Average number of CD16 receptors per NK cell')
title('Mean CD16 receptors per NK cell during RTX treatment')
set(gca,'FontSize',16)
xlim([0 nr_t_et])
subplot(1,2,2)
plot(Evec,'b','LineWidth', 2)
hold on
plot(Estarvec, 'k', 'LineWidth',2)
xlabel('Time (hours)')
legend('active NK cells', 'depleted NK cells')
ylabel('Number of NK cells')
title('NK cell phenotypes in time with constant RTX')
set(gca,'FontSize',16)
xlim([0 nr_t_et])

%% Plot the total CD20 receptor distribution over time and the pdf
% CD20 distributions
% Make this into a video
figure;
for i = 1:pstruct.nr_t_et
subplot(1,2,1)
hold off
plot(CD20edges(1:end-1), NCD200,'b', 'LineWidth', 2)
hold on
plot(CD20edges(1:end-1), NCD20(i,:), 'm','LineWidth', 0.5)
legend('t=0', ['t=', num2str(i), 'hr'])
%ylim([0 8000])
%xlim([0 0.1])
set(gca,'FontSize',16,'LineWidth',1.5,'Xscale', 'log')
ylabel('Number of CD20 receptors')
xlabel('\lambda*number of CD20 receptors per cell')
title(['Total CD20 distribution on tumor cells, t=', num2str(i), 'hours'])

subplot(1,2,2)
hold off
plot(CD20edges, pdCD20, 'b', 'LineWidth', 2)
%plot(CD20edges(1:end-1), pCD20(1,:),'r', 'LineWidth', 2)
hold on
plot(CD20edges(1:end-1), pCD20(i,:), 'm','LineWidth', 0.5)
legend('t=0', ['t=', num2str(i), 'hr'])
set(gca,'FontSize',16,'LineWidth',1.5,'Xscale', 'log')
ylabel('Probability')
xlabel('\lambda*number of CD20 receptors per cell')
title(['Normalized CD20 distribution on tumor cells, t=', num2str(i), 'hours'])
%xlim([ 0 0.1])
ylim([ 0 0.01])
drawnow
pause(0.01)
end


%% Plot the tumor cell trajectories for high and low CD16 and CD20 values

figure;
for i = 1:100
    if CD20samps_shift(i)>median(CD20samps_shift)
    plot(Tmat(:,i), 'r', 'LineWidth', 1)
    end
    if CD20samps_shift(i)<median(CD20samps_shift)
    plot(Tmat(:,i), 'b', 'LineWidth', 1)
    end
    hold on
    %legend('CD20- on tumor cell', 'CD20+ on tumor cell')
    xlabel('Time(hours)')
    ylabel('Relative number of tumor cells')
    title('Observed variability in tumor cell trajectories')
    set(gca,'FontSize',16)
    ylim([0 14])
end

figure;
for i = 2:101
    if CD16samps_shift(i)>median(CD16samps_shift) 
    plot(Tmat(:,i), 'g', 'LineWidth', 1)
    end
    if CD16samps_shift(i)<median(CD16samps_shift)
    plot(Tmat(:,i), 'y', 'LineWidth', 1)
    end
    hold on
    %legend('CD16 + on NK Cell', 'CD16 - on NK cell')
    xlabel('Time(hours)')
    ylabel('Relative number of tumor cells')
    title('Observed variability in tumor cell trajectories')
    set(gca,'FontSize',16)
    ylim([0 14])
end
