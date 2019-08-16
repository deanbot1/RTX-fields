%% how much CD16 upregulation is needed to overcome a given degree of CD20 loss?
clear all; close all

par = readtable('adcx_parameter_table.csv');
par.Properties.RowNames = par.paramnames;
v2s = @(vec)vect2struct(vec,par.Properties.RowNames); % local function converting parameter vector to parameter structure

CD16def = par{'CD16','value'};
CD20def = par{'CD20','value'};

CD16range = CD16def*(10.^[0:.1:3]);
%CD20range = CD20def*(10.^[-10:.25:0]);
CD20range = CD20def*(10.^[-10:.1:-7]);

[XX,YY] = meshgrid(CD20range,CD16range);


for i = 1:length(CD16range)
	for j = 1:length(CD20range)
		p = v2s(par.value);
		p.CD20 = XX(i,j);
		p.CD16 = YY(i,j);
		ZZ(i,j) = adcx_wrapper(p,p.RTX);
	end
end

%% plot the results
figure;
contourf(log10(XX),log10(YY),ZZ);
colorbar;
colormap copper
xlabel('log_{10}CD20');
ylabel('log_{10}CD16');
title('%ADCC');
caxis([0 100]);
grid on;
