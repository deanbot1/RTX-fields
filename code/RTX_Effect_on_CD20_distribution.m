% Does RTX modify the CD20 distribution in first line DLBCL patients?
% Before running this, run cellpopmodelscript.m with the correct parameter
% values

% We will answer this by looking at how the initial CD20 distribution of a
% single "patient" is predicted to change in time from our model.

% Visualize the CD20 distribution prior to treatment.

% % Tabulate the randomly drawn samples
% At time 0, T = 1 so the CD20 distribution should be the same as the
% CD20samps


% Make a matrix that is CD20 values * number of tumor cells at each time
CD20mat = zeros(size(Tmat));
for i = 1:nsamps
CD20mat(:,i) = Tmat(:,i).*CD20samps(i); 
end

% Make the same matrix for CD19 using E and Estar mat
Estarmat = max(max(Emat))-Emat;

CD16mat = zeros(size(Emat));
for i = 1:nsamps
    CD16matE(:,i) = Emat(:,i).*CD16samps(i);
    CD16matEstar(:,i) = Estarmat(:,i).*CD16samps(i)-delCD16(i);
end
CD16mat = CD16matE;

% eventually make this the normalized range of CD20 values (1 to 1 million)
% Now just let it

[NCD20] = histcounts(CD20samps, valuesCD20);
[NCD16] = histcounts(CD16samps, valuesCD16);

CD20edges = valuesCD20;
CD16edges = valuesCD16;

%% Look at mean CD20 and mean CD16 over time
%
Tvec = sum(Tmat,2);
Evec = sum(Emat,2);
Estarvec = sum(Estarmat,2);
CD20vec = sum(CD20mat,2); % total CD20 levels over time
CD16vec = sum(CD16mat,2);
% We want to plat the average CD20 and CD16 level per cell over time
% This should show that the average levels decrease over time
figure;
subplot(1,2,1)
plot(CD20vec./Tvec,'g','LineWidth', 2)
hold on
xlabel('Time (hours)')
ylabel('Average number of CD20 receptors per cell')
title('Mean CD20 receptors per tumor cell during RTX treatment')
set(gca,'FontSize',16)
xlim([0 nr_t_et])
subplot(1,2,2)
plot(Tvec,'b','LineWidth', 2)
hold on
xlabel('Time (hours)')
ylabel('Number of tumor cells')
title('Tumor cells in time with constant RTX')
set(gca,'FontSize',16)
xlim([0 nr_t_et])

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


%% Use the bin edges found from the initial distribution

for i = 1:nr_t_et

[NCD20(i,:)] = histcounts(CD20mat(i,:), CD20edges);
pCD20(i,:) = NCD20(i,:)./sum(NCD20(i,:));

[NCD16(i,:)] = histcounts(CD16mat(i,:), CD16edges);
pCD16(i,:) = NCD16(i,:)./sum(NCD16(i,:));

end


%% Plot the total CD20 receptor distribution over time and the pdf
% CD20 distributions
% Make this into a video
figure;
for i = 1:nr_t_et
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
%% Plot the final pdf
figure;
hold off
plot(CD20edges(1:end-1), pCD20(1,:),'r', 'LineWidth', 2)
hold on
plot(CD20edges(1:end-1), pCD20(end,:), 'm','LineWidth', 2)
legend('t=0', 't=500 hr')
set(gca,'FontSize',16,'LineWidth',1.5,'Xscale', 'log')
ylabel('Probability')
xlabel('\lambda*number of CD20 receptors per cell')
title('Normalized CD20 distribution on tumor cells, t=500 hours')
%xlim([ 0 0.1])
ylim([ 0 7e-3])
%% Not sure whats wrong here
figure;
for i = 1:2%nr_t_et

subplot(1,2,1)
hold off
plot(valuesCD16, freqCD16,'b-', 'LineWidth', 3)
hold on
plot(CD16edges(1:end-1), NCD16(1,:),'c', 'LineWidth', 2)
%plot(valuesCD16, NCD160,'b-', 'LineWidth', 3)
%plot(CD16edges(1:end-1), NCD16(i,:),'c-', 'LineWidth', 2)
legend('t=0', ['t=', num2str(i), 'hr'])
%ylim([0 8000])
%xlim([0 0.1])
set(gca,'FontSize',16,'LineWidth',1.5,'Xscale', 'log')
ylabel('Number of CD16 receptors')
xlabel('\lambda*number of CD16 receptors per NK cell')
title(['Total CD16 distribution on NK cells, t=', num2str(i), 'hours'])
xlim([5e2 10e3])

subplot(1,2,2)
hold off
plot(valuesCD16, pdCD16,'b-', 'LineWidth', 3)
hold on
plot(CD16edges(1:end-1), pCD16(i,:), 'c-','LineWidth', 2)
legend('t=0', ['t=', num2str(i), 'hr'])
set(gca,'FontSize',16,'LineWidth',1.5, 'Xscale', 'log')
ylabel('Probability')
xlabel('\lambda*number of CD16 receptors per NK cell')
title(['Normalized CD16 distribution on NK cells, t=', num2str(i), 'hours'])
%xlim([ 0 0.1])
%ylim([ 0 1])
xlim([5e2 10e3])
drawnow
pause(0.1)


end

%%
figure;
plot(valuesCD20, pdCD20, 'r', 'LineWidth',3)
hold on
set(gca,'FontSize',20,'LineWidth',1.5', 'Xscale', 'log')
xlim([0 1])

tbl1 = tabulate(CD20mat(1,:));
pdf1 = tbl1(:,3)./100;
ordered_vals1 = tbl1(:,1);

tblend = tabulate(CD20mat(end,:));
pdfend = tblend(:,3)./100;
ordered_valsend = tblend(:,1);

figure;
plot(ordered_vals1, pdf1, 'b', 'LineWidth',3);
hold on
plot(ordered_valsend, pdfend, 'r', 'LineWidth', 3);
legend('t=0 hr', 't=500 hr')
%xlim([10e-8 1])
set(gca,'FontSize',20,'LineWidth',1.5, 'Xscale', 'log')
xlabel('Number of receptors per cell')
ylabel('PDF')
title('Change in CD20 distribution in T cells')

 
