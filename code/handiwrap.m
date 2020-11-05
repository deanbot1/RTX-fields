function Tout = handiwrap(pin,tmax,varargin)
% calls adcx repeatedly for each concentration in R_conc and returns % ADCC
% vector. It's hard-wired to simulate for 4 hours unless tspan is passed

if isempty(tmax)
	tspan = [0:.1:4];
else
	tspan = [0:.1:tmax];
end

% first we generate the design matrix
parname = {}; parvec = {}; parnums = [];

Npars = length(varargin)/2; % it had better be even!

for k = 1:Npars 
	parname{k} = varargin{2*k-1};
	parvec{k} = reshape(varargin{2*k},[],1); % make sure it's a column to not anger matlab
	parnums(k) = length(parvec{k});
end

dFF = fullfact(parnums);


p_adcc = NaN(length(dFF),1);

for i = 1: length(dFF)
	p = pin;
	for k = 1:Npars
		p.(parname{k}) = parvec{k}(dFF(i,k));
	end
	
% 	disp([i,dFF(i,:)])
% 	p
    total_num_cells = p.T0*exp(p.g*max(tspan));
    [T,~,~,LDH,~,~]= adcx(tspan,p.T0,p.E0toT0,p.Estar0,p.g,p.r,p.kexp,p.gamma,...
    p.CD20,p.CD16,p.RTX*53.7496/7750,p.kon20,p.koff20,p.kon16,p.koff16,p.h,p.gammaPerf);
    p_adcc(i,1) = 100*LDH(end)/total_num_cells;
end
Tout = table(p_adcc,'VariableNames',{'PCT_ADCC'});
for k = 1:Npars
	Tout.(parname{k}) = parvec{k}(dFF(:,k));
end
