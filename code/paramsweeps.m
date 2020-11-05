%% parameter sweeps test behavior over specified parameters...

close all; clear all; 
Tpar = readtable('adcx_parameter_results_fat.csv');
varname = 'E4_F158onSUDHL4' ;

for j = 1:length(Tpar.name)
	pdef.(Tpar.name{j}) = Tpar.(varname)(j);
end

%% run a bunch of stuff

Tout = handiwrap(pdef,[],'CD20',pdef.CD20*[0:1:10],'CD16',pdef.CD16*[0:.5:2],'RTX',[.05 .5 5 50]);
figure; handiplot(Tout);
