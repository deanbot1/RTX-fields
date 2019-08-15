clear all; close all;

syms koff_RA kon_RA koff_FR kon_FR F R A FR RA FRA Ftot Rtot Atot positive

dRdt = koff_RA*RA + koff_FR*FR - kon_RA*R*A - kon_FR*F*R;
dAdt = koff_RA*(RA + FRA) - kon_RA*(R*A + FR*A);
dFdt = koff_FR*(FR + FRA) - kon_FR*(F*R + F*RA);
dRAdt = kon_RA*R*A - koff_RA*RA - kon_FR*F*RA + koff_FR*FRA;
dFRdt = kon_FR*F*R - koff_FR*FR - kon_RA*FR*A + koff_RA*FRA;
dFRAdt = kon_FR*F*RA + kon_RA*FR*A - (koff_FR+koff_RA)*FRA;

dRdt = subs(dRdt,RA,Atot-A-FRA);
dAdt = subs(dAdt,RA,Atot-A-FRA);
dFdt = subs(dFdt,RA,Atot-A-FRA);
dRAdt = subs(dRAdt,RA,Atot-A-FRA);
dFRdt = subs(dFRdt,RA,Atot-A-FRA);
dFRAdt = subs(dFRAdt,RA,Atot-A-FRA);

dRdt = subs(dRdt,FR,Ftot-F-FRA);
dAdt = subs(dAdt,FR,Ftot-F-FRA);
dFdt = subs(dFdt,FR,Ftot-F-FRA);
dRAdt = subs(dRAdt,FR,Ftot-F-FRA);
dFRdt = subs(dFRdt,FR,Ftot-F-FRA);
dFRAdt = subs(dFRAdt,FR,Ftot-F-FRA);

eqns = [dRdt == 0, dAdt == 0, dFdt == 0, dRAdt == 0, dFRdt == 0, dFRAdt == 0];%
%Ftot == F + FR + FRA, Atot == A + RA + FRA, Rtot == R + RA + FR + FRA];
vars = [FRA,FR,RA,F,R,A];

S = solve(eqns,vars,'real',true,'IgnoreAnalyticConstraints',true,'ReturnConditions',true);