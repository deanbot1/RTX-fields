% In theory: how does drawing from the observed CD20 and CD16 distribution
% effect the results of the tumor control via ADCC?

% Run the cell popmodel script first with drawing from the distribution,
% then run it sampling repeatedly from the mean alone.

%% Run the loop that iterates through the CD16 and CD20 samples but that 
%uses the mean of CD20 and CD16

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
        (mean(CD20samps)/gamma),(mean(CD16samps)/lambda),RTX,kon20,koff20,kon16,koff16,gamma_perf);
    i
    Tmatc(:,i)     = Ti;
    Ematc(:,i)     = Ei;
    LDHmatc(:,i)   = LDHi;
    perfmatc(:,i)  = perfi;
    delCD16c(i,1)  = CPXi;
    
end

% For a given CD16, CD20 distribution, here are the expected LDH and perf
% trajectories over time
LDHvec = sum(LDHmat,2);
perfvec = sum(perfmat,2);

%% Plotting
figure(5)
hold off
plot(sum(Tmatc,2),'LineWidth',2)
hold on
plot(sum(Ematc,2),'LineWidth',2)
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
plot(sum(Tmatc,2),'LineWidth',2)
hold on
plot(sum(Ematc,2),'LineWidth',2)
plot(LDHvec,'LineWidth',2)
xlabel('Time (hours)')
ylabel('Number of cells/molecules')
title('Long Range Dynamics using mean receptor levels')
legend('Target cells','Effector cells (CD16)','LDH')
legend boxoff
set(gca,'FontSize',16)
%% Plot individual tumor trajectories for CD20-


figure;
for i = 1:100
    if CD20samps(i)>median(CD20samps)
    plot(Tmat(:,i), 'r', 'LineWidth', 1)
    end
    if CD20samps(i)<median(CD20samps)
    plot(Tmat(:,i), 'b', 'LineWidth', 1)
    end
    hold on
    legend('CD20- on tumor cell', 'CD20+ on tumor cell')
    xlabel('Time(hours)')
    ylabel('Relative number of tumor cells')
    title('Observed variability in tumor cell trajectories')
    set(gca,'FontSize',16)
    ylim([0 14])
end
%%
figure;
for i = 2:101
    if CD16samps(i)>median(CD16samps) 
    plot(Tmat(:,i), 'g', 'LineWidth', 1)
    end
    if CD16samps(i)<median(CD16samps)
    plot(Tmat(:,i), 'y', 'LineWidth', 1)
    end
    hold on
    legend('CD16 + on NK Cell', 'CD16 - on NK cell')
    xlabel('Time(hours)')
    ylabel('Relative number of tumor cells')
    title('Observed variability in tumor cell trajectories')
    set(gca,'FontSize',16)
    ylim([0 14])
end
%%
figure;
for i = 95:195
    if CD16samps(i)>median(CD16samps) && CD20samps(i)< median(CD20samps)
        i
    plot(Tmat(:,i), 'r', 'LineWidth', 1.5)
    else
    plot(Tmat(:,i), 'k', 'LineWidth', 1)
    end
    hold on
    legend( 'CD16 + on NK cell, CD20- on tumor', 'others')
    xlabel('Time(hours)')
    ylabel('Relative number of tumor cells')
    title('Observed variability in tumor cell trajectories')
    set(gca,'FontSize',16)
    ylim([0 14])
end


