%% tornado_plot of sensitivities of model at the default values for other parameters
clear all; close all; 

par = readtable('adcx_parameter_table.csv');
par.low = par.value/2;
par.high = par.value*2;
par.Properties.RowNames = par.paramnames;
v2s = @(vec)vect2struct(vec,par.Properties.RowNames); % local function converting parameter vector to parameter structure

paramstochange = par.paramnames(find(par.sens));
Npar = length(paramstochange);

R_conc = 10.^[-4:3];
modelfun = @(p)adcx_wrapper(p,R_conc)
p = v2s(par.value);
defval = modelfun(p)';

for j = 1:Npar
	plow = p; phigh = p;
	plow.(paramstochange{j}) = par{paramstochange{j},'low'};
	loval(j,:) = modelfun(plow)';
	phigh.(paramstochange{j}) = par{paramstochange{j},'high'};
	hival(j,:) = modelfun(phigh)';
end

%% now plot the results
figure;
sq_diffs_low = (defval-loval).^2;
sq_diffs_high = (hival - defval).^2;
diffs_lows = sum(sq_diffs_low,2);
diffs_high = sum(sq_diffs_high,2);
diffs = diffs_high + diffs_lows ; 
[~,jsort] = sort(diffs);
def_vals_zeros = zeros(1,Npar);
for j = 1:Npar
	%plot([loval(jsort(j)) hival(jsort(j))],j*[1 1],'b-o','linewidth',4,'MarkerFaceColor','b');
	plot([-diffs_lows(jsort(j)) diffs_high(jsort(j))],j*[1 1],'b-o','linewidth',4,'MarkerFaceColor','b');
    plot(-diffs_lows(jsort(j)),j,'bo','MarkerFaceColor','g');
	plot(diffs_high(jsort(j)),j,'bo','MarkerFaceColor','r');
	plot(def_vals_zeros,j,'ko','MarkerFaceColor','k','MarkerSize',5);
	hold on;
	if diffs_lows(jsort(j)) < diffs_high(jsort(j))
		loalign = 'right'; highalign = 'left';
	else
		loalign = 'left';highalign = 'right';
	end
	%text(defval,j,num2str(par{paramstochange(jsort(j)),'value'}),'horizontalalignment','center','verticalalignment','bottom');
	text(diffs_lows(jsort(j)),j,['  ' num2str(par{paramstochange(jsort(j)),'low'}) '  '],'HorizontalAlignment',loalign,'color','g');
	text(diffs_high(jsort(j)),j,['  ' num2str(par{paramstochange(jsort(j)),'high'}) '  '],'HorizontalAlignment',highalign,'color','r');
end

hold on
plot(0*[1 1],[0 Npar+0.5],'k-');
set(gca,'Ytick',[1:Npar],'Yticklabel',paramstochange(jsort));
%set(gca,'Xlim',[-20 120],'Ylim',[0.5 Npar+0.5]);
xlabel('% ADCC @ 4h');
grid on
set(gca,'FontSize',16)