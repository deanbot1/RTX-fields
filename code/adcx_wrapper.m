function p_adcc = adcx_wrapper(p,R_conc)
p_adcc = zeros(length(R_conc),1);
R_conc = 53.7496/7750 * R_conc;
for i = 1: length(R_conc)
    total_num_cells = p.T0*exp(p.g*p.tf_et);
    %p.g = 0;
    %tf_mol,tf_et,nr_t_mol,nr_t_et,T0,E0toT0,Estar0,g,r,kexp,gamma,CD20,CD16,RTX,kon20,koff20,kon16,koff16,gamma_perf)
    [T,~,~,LDH,~,~]= adcx(p.tf_mol,p.tf_et,p.nr_t_mol,p.nr_t_et,p.T0,p.E0toT0,p.Estar0,p.g,p.r,p.kexp,p.gamma,...
    p.CD20,p.CD16,R_conc(i),p.kon20,p.koff20,p.kon16,p.koff16,p.h,p.gamma_perf);
    p_adcc(i,1) = 100*LDH(end)/total_num_cells;
    
end