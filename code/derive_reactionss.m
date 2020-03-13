clear all; close all;

syms koff_RA kon_RA koff_FR kon_FR F R A FR RA FRA Ftot Rtot Atot positive

%dRdt = koff_RA*RA + koff_FR*FR - kon_RA*R*A - kon_FR*F*R;   % free RTX (comment if RTX conc is clamped and not affected by binding)
dAdt = koff_RA*(RA + FRA) - kon_RA*(R*A + FR*A);             % free (membrane bound) CD20 ?(A for antigen, I believe)
dFdt = koff_FR*(FR + FRA) - kon_FR*(F*R + F*RA);             % free (membrane bound) FcgRIIIa receptor 
dRAdt = kon_RA*R*A - koff_RA*RA - kon_FR*F*RA + koff_FR*FRA; % RTX:CD20 complex (membrane bound)
dFRdt = kon_FR*F*R - koff_FR*FR - kon_RA*FR*A + koff_RA*FRA; % FcgRIIIa:RTX complex (membrane bound)
dFRAdt = kon_FR*F*RA + kon_RA*FR*A - (koff_FR+koff_RA)*FRA;  % FcgRIIIa:RTX:CD20 complex (membrane bound << both cells >>)


eqns = [dAdt == 0, dFdt == 0, dRAdt == 0, dFRdt == 0, dFRAdt == 0,...
	Ftot == F + FR + FRA, Atot == A + RA + FRA];
vars = [FRA,F,FR,A,RA];

S = solve(eqns,vars,'real',true,'IgnoreAnalyticConstraints',true,'ReturnConditions',true);
save('reactionss.mat');

%% now test all of the solutions to see which are strictly positive and match up with the ODE-based steady state solution
% copy the winning solution below.  And then copy it to SSadcc.m, I
% think...

FRAF = @(R,Atot,Ftot,kon_FR,koff_FR,kon_RA,koff_RA)...
[
 Atot + ((koff_RA + R*kon_RA)*(koff_RA*koff_FR - (Atot^2*R^2*kon_RA^2*kon_FR^2 - 2*Atot*Ftot*R^2*kon_RA^2*kon_FR^2 + 2*Atot*R^3*kon_RA^2*kon_FR^2 + 2*Atot*R^2*koff_RA*kon_RA*kon_FR^2 + 2*Atot*R^2*koff_FR*kon_RA^2*kon_FR + 2*Atot*R*koff_RA*koff_FR*kon_RA*kon_FR + Ftot^2*R^2*kon_RA^2*kon_FR^2 + 2*Ftot*R^3*kon_RA^2*kon_FR^2 + 2*Ftot*R^2*koff_RA*kon_RA*kon_FR^2 + 2*Ftot*R^2*koff_FR*kon_RA^2*kon_FR + 2*Ftot*R*koff_RA*koff_FR*kon_RA*kon_FR + R^4*kon_RA^2*kon_FR^2 + 2*R^3*koff_RA*kon_RA*kon_FR^2 + 2*R^3*koff_FR*kon_RA^2*kon_FR + R^2*koff_RA^2*kon_FR^2 + 4*R^2*koff_RA*koff_FR*kon_RA*kon_FR + R^2*koff_FR^2*kon_RA^2 + 2*R*koff_RA^2*koff_FR*kon_FR + 2*R*koff_RA*koff_FR^2*kon_RA + koff_RA^2*koff_FR^2)^(1/2) + R*koff_RA*kon_FR + R*koff_FR*kon_RA + R^2*kon_RA*kon_FR - Atot*R*kon_RA*kon_FR + Ftot*R*kon_RA*kon_FR))/(2*R*kon_RA*(koff_RA*kon_FR + R*kon_RA*kon_FR))
 ];
