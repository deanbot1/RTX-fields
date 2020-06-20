function [] = sensitivity_analysis(par_name)
par_name = 'kon16' ;
parameter_filename = 'adcx_parameter_table.csv';
par = readtable(parameter_filename);
par.Properties.RowNames = par.paramnames;
v2s = @(vec)vect2struct(vec,par.Properties.RowNames);

data_path = '../data/';
data_filenames = {'Herter_4B_V158onSU-DHL4.csv'};
data = readtable([data_path data_filenames{1}]);
figure;
semilogx(data.Var1,data.Var2,'-k.','MarkerSize',20);
ylim([0,100]);
hold on

R_conc = 10.^[-4:3];
bounds = {'low','high'};
iter = [par{par_name,bounds}];
p = v2s(par.value);
% for i = [iter(1) (iter(2)-iter(1))/2 iter(2)]
%     p.par_name = i
%     dead = adcx_wrapper(p,R_conc);
%     semilogx(R_conc,dead);
% end
legend('data','low','mid','high')
end