clear all; close all;

syms ooh koff_RA kon_RA koff_FR kon_FR F R A FR RA FRA Ftot Rtot Atot positive

% below, ooh = 1/h the height of the synapse
%dRdt = ooh*(koff_RA*RA + koff_FR*FR - kon_RA*R*A - kon_FR*F*R);   % free RTX (comment if RTX conc is clamped and not affected by binding)
dAdt = koff_RA*(RA + FRA) - kon_RA*(R*A + ooh*FR*A);               % free (membrane bound) CD20 ?(A for antigen, I believe)
dFdt = koff_FR*(FR + FRA) - kon_FR*(F*R + ooh*F*RA);               % free (membrane bound) FcgRIIIa receptor 
dRAdt = kon_RA*R*A - koff_RA*RA - ooh*kon_FR*F*RA + koff_FR*FRA;   % RTX:CD20 complex (membrane bound)
dFRdt = kon_FR*F*R - koff_FR*FR - ooh*kon_RA*FR*A + koff_RA*FRA;   % FcgRIIIa:RTX complex (membrane bound)
dFRAdt = ooh*(kon_FR*F*RA + kon_RA*FR*A) - (koff_FR+koff_RA)*FRA;  % FcgRIIIa:RTX:CD20 complex (membrane bound << both cells >>)


eqns = [dAdt == 0, dFdt == 0, dRAdt == 0, dFRdt == 0, dFRAdt == 0,...
	Ftot == F + FR + FRA, Atot == A + RA + FRA];
vars = [FRA,F,FR,A,RA];

S = solve(eqns,vars,'real',true,'IgnoreAnalyticConstraints',true,'ReturnConditions',true);
save('reactionss.mat');

%% now test all of the solutions to see which are strictly positive and match up with the ODE-based steady state solution
% copy the winning solution below.  
P = readtable('parameter_estimation_results.csv','ReadRowNames',true);
xspan = 10.^[-4:0.5:6]'; % xspan approx from parameter_estimation.m
figure;
Y = zeros(length(S.FRA),length(xspan));

for j = 1:length(S.FRA)
	disp(sprintf('****** trying solution # %d ***********',j));
	SS{j} = @(RTX)eval(subs(S.FRA(j),[R,Ftot,Atot,ooh,koff_RA,kon_RA,koff_FR,kon_FR],...
		[RTX,P{'CD16','value'},P{'CD20_Z138','value'},1,P{'koff20','value'},P{'kon20','value'},P{'koff16','value'},P{'kon16_V158','value'}]));
	for k = 1:length(xspan)
		Y(j,k) = SS{j}(xspan(k));
	end
	h = plot(xspan,Y(j,:),'-'); hold on
	text(xspan(1),Y(j,1),num2str(j),'HorizontalAlignment','right','Color',get(h,'Color'));
end

set(gca,'Xscale','log','Yscale','log');
xlabel('[RTX]');
ylabel('[FRA trimer]');
%set(gca,'YLim',[0 1.1*max(Y(1,:))]);

% when I run it, the only graph that makes sense is for solution #1, so
% let's output that one!

disp(S.FRA(1));


%% final solution
FRAF = @(R,Ftot,Atot,ooh,koff_RA,kon_RA,koff_FR,kon_FR)...
[
Atot + ((koff_RA + R*kon_RA).*(koff_RA*koff_FR - (Atot^2*R.^2*kon_RA^2*kon_FR^2*ooh^2 - 2*Atot*Ftot*R.^2*kon_RA^2*kon_FR^2*ooh^2 + 2*Atot*R.^3*kon_RA^2*kon_FR^2*ooh + 2*Atot*R.^2*koff_RA*kon_RA*kon_FR^2*ooh + 2*Atot*R.^2*koff_FR*kon_RA^2*kon_FR*ooh + 2*Atot*R*koff_RA*koff_FR*kon_RA*kon_FR*ooh + Ftot^2*R.^2*kon_RA^2*kon_FR^2*ooh^2 + 2*Ftot*R.^3*kon_RA^2*kon_FR^2*ooh + 2*Ftot*R.^2*koff_RA*kon_RA*kon_FR^2*ooh + 2*Ftot*R.^2*koff_FR*kon_RA^2*kon_FR*ooh + 2*Ftot*R*koff_RA*koff_FR*kon_RA*kon_FR*ooh + R.^4*kon_RA^2*kon_FR^2 + 2*R.^3*koff_RA*kon_RA*kon_FR^2 + 2*R.^3*koff_FR*kon_RA^2*kon_FR + R.^2*koff_RA^2*kon_FR^2 + 4*R.^2*koff_RA*koff_FR*kon_RA*kon_FR + R.^2*koff_FR^2*kon_RA^2 + 2*R*koff_RA^2*koff_FR*kon_FR + 2*R*koff_RA*koff_FR^2*kon_RA + koff_RA^2*koff_FR^2).^(1/2) + R*koff_RA*kon_FR + R*koff_FR*kon_RA + R.^2*kon_RA*kon_FR - Atot*R*kon_RA*kon_FR*ooh + Ftot*R*kon_RA*kon_FR*ooh))./(2*R*kon_RA.*(koff_RA*kon_FR*ooh + R*kon_RA*kon_FR*ooh))
 ];

YY = FRAF(xspan,P{'CD16','value'},P{'CD20_Z138','value'},1,P{'koff20','value'},P{'kon20','value'},P{'koff16','value'},P{'kon16_V158','value'});

plot(xspan,YY,'k:','LineWidth',2);  % just to double check I didn't mess anything up when I vectorized the expression
