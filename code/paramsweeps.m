%% parameter sweeps test behavior when sweeping over specified parameters...

close all; clear all; 
Tpar = readtable('adcx_parameter_results_fat.csv');
varname = 'E4_F158onSUDHL4' ;

for j = 1:length(Tpar.name)
	pdef.(Tpar.name{j}) = Tpar.(varname)(j);
end

%% run a bunch of stuff

Tout = handiwrap(pdef,[],'CD20',pdef.CD20*[0:1:10],'CD16',pdef.CD16*[0:.5:2],'RTX',[.05 .5 5 50]);
figure; handiplot(Tout);

%% EtoT dependency
Tout = handiwrap(pdef,[],'E0toT0',50*2.^[-5:0],'RTX',[5]);
figure; handiplot(Tout);
set(gca,'Xscale','log');
% overlay observations
Tdat = readtable('../data/WangetAl/Wang-Fig8-ET-ADCC.csv');
hold on;
plot(Tdat.EtoT,Tdat.pctADCC,'ro','MarkerFaceColor','r');
set(gca,'Xtick',Tdat.EtoT);
set(gca,'Ylim',[0 100]);
axis tight

%% CD20 dependency
Tout = handiwrap(pdef,[],'CD20',[0:10:300],'RTX',[5],'gamma',pdef.gamma*[2 20 200]);
figure; handiplot(Tout);
% now overlay data
Tdat = readtable('../data/TsaiEtAl/Tsai_CD20_ADCC.csv');
hold on; 
plot(Tdat.CD20,Tdat.PCT_LYSIS,'ro','MarkerFaceColor','r');
set(gca,'Ylim',[0 100]);