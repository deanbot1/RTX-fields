%% how much CD16 upregulation is needed to overcome a given degree of CD20 loss?
clear all; close all

par = readtable('parameter_estimation_results.csv');
par.Properties.RowNames = par.paramnames;
for j = 1:height(par);
	pdef.(par.target{j}) = par.bestest(j);
end

%v2s = @(vec)vect2struct(vec,par.Properties.RowNames); % local function converting parameter vector to parameter structure

CD16def = pdef.CD16;
CD20def = pdef.CD20;

CD16range = CD16def*(10.^[-1:.1:3]);
%CD20range = CD20def*(10.^[-10:.25:0]);
CD20range = CD20def*(10.^[-3:.1:1]);

[XX,YY] = meshgrid(CD20range,CD16range);


for i = 1:length(CD16range)
	for j = 1:length(CD20range)
		p = pdef;
		p.CD20 = XX(i,j);
		p.CD16 = YY(i,j);
		ZZ(i,j) = adcx_wrapper(p,p.RTX);
	end
end

%% plot the results
figure;
[C,h]=contourf(log10(XX),log10(YY),ZZ);
clabel(C,h);
colorbar;
colormap copper
xlabel('log_{10}CD20');
ylabel('log_{10}CD16');
title('%ADCC');
caxis([0 100]);
grid on;
hold on
plot(log10(CD20def*[1 1]),log10([min(CD16range) max(CD16range)]),'w-'); hold on
plot(log10([min(CD20range) max(CD20range)]),log10(CD16def*[1 1]),'w-');

