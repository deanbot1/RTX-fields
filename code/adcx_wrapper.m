function p_adcc = adcx_wrapper(p,R_conc,tspan)
% calls adcx repeatedly for each concentration in R_conc and returns % ADCC
% vector. It's hard-wired to simulate for 4 hours unless tspan is passed
if nargin < 3
	tspan = [0:.1:4];
end

p_adcc = zeros(length(R_conc),1);
R_conc = 53.7496/7750 * R_conc;  % convert ug/ml to M? nM/L?
for i = 1: length(R_conc)
    total_num_cells = p.T0*exp(p.g*max(tspan));
    [T,~,~,LDH,~,~]= adcx(tspan,p.T0,p.E0toT0,p.Estar0,p.g,p.r,p.kexp,p.gamma,...
    p.CD20,p.CD16,R_conc(i),p.kon20,p.koff20,p.kon16,p.koff16,p.h,p.gammaPerf);
    p_adcc(i,1) = 100*LDH(end)/total_num_cells;
end