%% tornado_plot of sensitivities of model at the default values for other parameters
clear all; close all; 

par = readtable('adcx_parameter_table.csv');
par.low = par.value/2;
par.high = par.value*2;
par.Properties.RowNames = par.paramnames;
v2s = @(vec)vect2struct(vec,par.Properties.RowNames); % local function converting parameter vector to parameter structure

paramstochange = par.paramnames(find(par.sens));
Npar = length(paramstochange);

R_conc = par{'RTX','value'};
modelfun = @(p)adcx_wrapper(p,R_conc)
p = v2s(par.value);
defval = modelfun(p);

for j = 1:Npar
	plow = p; phigh = p;
	plow.(paramstochange{j}) = par{paramstochange{j},'low'};
	loval(j,1) = modelfun(plow);
	phigh.(paramstochange{j}) = par{paramstochange{j},'high'};
	hival(j,1) = modelfun(phigh);
end

%% now plot the results
figure;
diffs = abs(loval-hival);
[~,jsort] = sort(diffs);
for j = 1:Npar
	plot([loval(jsort(j)) hival(jsort(j))],j*[1 1],'b-o','linewidth',4,'MarkerFaceColor','b');
	plot(loval(jsort(j)),j,'bo','MarkerFaceColor','g');
	plot(hival(jsort(j)),j,'bo','MarkerFaceColor','r');
	plot(defval,j,'ko','MarkerFaceColor','k','MarkerSize',5);
	hold on;
	if loval(jsort(j)) < hival(jsort(j))
		loalign = 'right'; highalign = 'left';
	else
		loalign = 'left';highalign = 'right';
	end
	text(defval,j,num2str(par{paramstochange(jsort(j)),'value'}),'horizontalalignment','center','verticalalignment','bottom');
	text(loval(jsort(j)),j,['  ' num2str(par{paramstochange(jsort(j)),'low'}) '  '],'HorizontalAlignment',loalign,'color','g');
	text(hival(jsort(j)),j,['  ' num2str(par{paramstochange(jsort(j)),'high'}) '  '],'HorizontalAlignment',highalign,'color','r');
end


plot(defval*[1 1],[0 Npar+0.5],'k-');
set(gca,'Ytick',[1:Npar],'Yticklabel',paramstochange(jsort));
set(gca,'Xlim',[-20 120],'Ylim',[0.5 Npar+0.5]);
xlabel('% ADCC @ 4h');
grid on
