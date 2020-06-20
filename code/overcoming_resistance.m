%% how much CD16 upregulation is needed to overcome a given degree of CD20 loss?
clear all; close all

Xname = 'CD20'; % specify X parameter name for the sensitivity analysis
Yname = 'gamma'; % Y parameter name
RTXconc = 100; % Rituximab concentration

par = readtable('adcx_parameter_results_fat.csv');
par.Properties.RowNames = par.name;
for j = 1:height(par);
	pdef_V158.(par.name{j}) = par.E2_V158onSUDHL4(j);
	pdef_F158.(par.name{j}) = par.E4_F158onSUDHL4(j);
end

%v2s = @(vec)vect2struct(vec,par.Properties.RowNames); % local function converting parameter vector to parameter structure


for kk = 1:2
	switch kk
		case 1
			pdef=pdef_V158;
			titl='V158';
		case 2
			pdef=pdef_F158;
			titl='F158';
	end

Ydef = pdef.(Yname);
Xdef = pdef.(Xname);

Yrange = Ydef*(10.^[-1:.1:2]);
%Xrange = Xdef*(10.^[-10:.25:0]);
Xrange = Xdef*(10.^[-2:.1:1]);

[XX,YY] = meshgrid(Xrange,Yrange);


for i = 1:length(Yrange)
	for j = 1:length(Xrange)
		p = pdef;
		p.(Xname) = XX(i,j);
		p.(Yname) = YY(i,j);
		ZZ(i,j) = adcx_wrapper(p,RTXconc);
	end
end

% plot the results
figure;
[C,h]=contourf(log10(XX),log10(YY),ZZ);
clabel(C,h);
colorbar;
colormap copper
xlabel(['log_{10}' Xname]);
ylabel(['log_{10}' Yname]);
title(['%ADCC (' titl ')']);
caxis([0 100]);
grid on;
hold on
plot(log10(Xdef*[1 1]),log10([min(Yrange) max(Yrange)]),'w-'); hold on
plot(log10([min(Xrange) max(Xrange)]),log10(Ydef*[1 1]),'w-');
print(sprintf('../out/%s_vs_%s_%s.png',Xname,Yname,titl),'-dpng');
end