function [Tpar,Texp] = read_par_expt(parameter_filename) 
%READ_PAR_EXPT reads in parameter_filename
% and returns two tables
% Tpar, the parameter settings table
% Texp, the experimental settings table

Tbig = readtable(parameter_filename);

ires = find(strcmp(Tbig.name,'RESERVED'));
Tpar = Tbig(1:ires-1,:);
Tpar.default = arrayfun(@(s)str2num(s{1}),Tpar.default);



Texp = Tbig(ires+1:end,:);
Texp.Properties.RowNames = Texp.name;