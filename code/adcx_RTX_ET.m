function p_adcc = adcx_RTX_ET(p,RTX_ET,tspan)
% calls adcx repeatedly for each concentration in RTX_ET(:,1) and paired
% EtoT ratio in RTX_ET(:,2)
% and returns % ADCC
% vector. It's hard-wired to simulate for 3 hours unless tspan is passed
% 3 hours is what Akito Nakamura (Takeda) reported using

if nargin < 3
	tspan = [0:.1:3];
end

R_conc = 53.7496/7750 * RTX_ET(:,1); % convert ug/ml to M? nM/L?
EtoT = RTX_ET(:,2); % E to T ratio
p_adcc = zeros(length(R_conc),1);
for i = 1: length(R_conc)
    total_num_cells = p.T0*exp(p.g*max(tspan));
    [T,~,~,LDH,~,~]= adcx(tspan,p.T0,EtoT(i),p.Estar0,p.g,p.r,p.kexp,p.gamma,...
    p.CD20,p.CD16,R_conc(i),p.kon20,p.koff20,p.kon16,p.koff16,p.h,p.gammaPerf);
    p_adcc(i,1) = 100*LDH(end)/total_num_cells;
end
%if max(R_conc)>1 & EtoT(1) > 5, keyboard,end;