%% Plotting profie likelihoods

par = 1;
% results_fileroot = 'adcx_prof_likelihoods';
results_fileroot = 'adcx_prof_likelihoods_center';
load([results_fileroot '/prof_likelihoods_bw0_runs20_par' num2str(par) '.mat'])

threshold = chi2inv(0.95,length(list_names))/2 + fval;
for k = par
    
figure(k)
hold on
plot(pj(k,:),obj(k,:),'*-','LineWidth',2)
plot(pbestlog(k),fval,'r*','LineWidth',2)
plot([pj(k,1) pj(k,end)],[threshold threshold],'r--','LineWidth',2)
% ylim([800 1000])
xlabel(list_names{k})
ylabel('fval')
set(gca,'FontSize',16)

end

saveas(gcf,[results_fileroot '/Par' num2str(par) '.png'])




