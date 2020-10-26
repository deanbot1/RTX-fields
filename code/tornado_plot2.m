%% tornado_plot of sensitivities of model at the default values for other parameters
clear all; close all; 

Tpar = readtable('adcx_parameter_results_fat.csv');  % need to update this with veronica's latest numbers!

for j = 1:height(Tpar)
	pbest.(Tpar.name{j}) = Tpar.E1_V158onZ138(j); % pick which experiment to center tornado around
end

%% test the main function
dp = struct('CD20',0.01); % for example, dp = struct('CD20',0.01) sets up the question, how much fold increase in 'CD16' is needed to offset 2 logs drop in CD20?
dx = how_much_dx_to_offset_dp(pbest,dp,'gamma',0.1,4);

ofun = @(pstruct)how_much_dx_to_offset_dp(pstruct,dp,'gamma',0.1,4) ; % set up anonymous function to run analysis on

%% next step, read in quantiles and step through them for each paramter in the quantiles table.
Tfull = readtable('quantiles_bootstraps_skinny.csv'); Tfull.Properties.RowNames = Tfull.Var1;
qlevels = Tfull{'name',2:end};
Tq = readtable('quantiles_bootstraps_skinny.csv','header',1);



%% outdated stuff

par = readtable('adcx_parameter_table.csv');
par.low = par.value/2;
par.high = par.value*2;
par.Properties.RowNames = par.paramnames;
v2s = @(vec)pvec2struct(vec,par.Properties.RowNames); % local function converting parameter vector to parameter structure

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
