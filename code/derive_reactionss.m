clear all; close all;

syms koff_RA kon_RA koff_FR kon_FR F R A FR RA FRA Ftot Rtot Atot positive

%dRdt = koff_RA*RA + koff_FR*FR - kon_RA*R*A - kon_FR*F*R;
dAdt = koff_RA*(RA + FRA) - kon_RA*(R*A + FR*A);
dFdt = koff_FR*(FR + FRA) - kon_FR*(F*R + F*RA);
dRAdt = kon_RA*R*A - koff_RA*RA - kon_FR*F*RA + koff_FR*FRA;
dFRdt = kon_FR*F*R - koff_FR*FR - kon_RA*FR*A + koff_RA*FRA;
dFRAdt = kon_FR*F*RA + kon_RA*FR*A - (koff_FR+koff_RA)*FRA;


eqns = [dAdt == 0, dFdt == 0, dRAdt == 0, dFRdt == 0, dFRAdt == 0,...
	Ftot == F + FR + FRA, Atot == A + RA + FRA];
vars = [FRA,F,FR,A,RA];

S = solve(eqns,vars,'real',true,'IgnoreAnalyticConstraints',true,'ReturnConditions',true);
save('reactionss.mat');

%% now work with this...