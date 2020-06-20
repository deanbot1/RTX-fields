	function SStrimer = SSadcc(RTX,CD16,CD20,k16RTXon,k16RTXoff,kRTX20on,kRTX20off,h)
% SSadcc returns steady state ADCC trimer levels given RTX levels, CD16 on
% effector cells, CD20 on target cells, and the on/off rates.  Model
% assumes RTX is not significantly depleted by the reaction.
% NOW VECTORIZED, which isn't really used by anyone...

		
	R = RTX;
	Ftot = CD16;
	Atot = CD20;
	kon_FR = k16RTXon;
	koff_FR = k16RTXoff;
	kon_RA = kRTX20on;
	koff_RA = kRTX20off;
	ooh = 1/h; 
	
	% SStrimer expression is derived in derive_reactionss.m, line 53 or so
	SStrimer = Atot + ((koff_RA + R*kon_RA).*(koff_RA*koff_FR - (Atot^2*R.^2*kon_RA^2*kon_FR^2*ooh^2 - 2*Atot*Ftot*R.^2*kon_RA^2*kon_FR^2*ooh^2 + 2*Atot*R.^3*kon_RA^2*kon_FR^2*ooh + 2*Atot*R.^2*koff_RA*kon_RA*kon_FR^2*ooh + 2*Atot*R.^2*koff_FR*kon_RA^2*kon_FR*ooh + 2*Atot*R*koff_RA*koff_FR*kon_RA*kon_FR*ooh + Ftot^2*R.^2*kon_RA^2*kon_FR^2*ooh^2 + 2*Ftot*R.^3*kon_RA^2*kon_FR^2*ooh + 2*Ftot*R.^2*koff_RA*kon_RA*kon_FR^2*ooh + 2*Ftot*R.^2*koff_FR*kon_RA^2*kon_FR*ooh + 2*Ftot*R*koff_RA*koff_FR*kon_RA*kon_FR*ooh + R.^4*kon_RA^2*kon_FR^2 + 2*R.^3*koff_RA*kon_RA*kon_FR^2 + 2*R.^3*koff_FR*kon_RA^2*kon_FR + R.^2*koff_RA^2*kon_FR^2 + 4*R.^2*koff_RA*koff_FR*kon_RA*kon_FR + R.^2*koff_FR^2*kon_RA^2 + 2*R*koff_RA^2*koff_FR*kon_FR + 2*R*koff_RA*koff_FR^2*kon_RA + koff_RA^2*koff_FR^2).^(1/2) + R*koff_RA*kon_FR + R*koff_FR*kon_RA + R.^2*kon_RA*kon_FR - Atot*R*kon_RA*kon_FR*ooh + Ftot*R*kon_RA*kon_FR*ooh))./(2*R*kon_RA.*(koff_RA*kon_FR*ooh + R*kon_RA*kon_FR*ooh));

	
	
% sbioloadproject('ADCC_reactions.sbproj');
% m1.Species(1).ConstantAmount = RTXconst;
% if length(RTX)==1
% 	varin = sbiovariant('v1',{'species','RTX','InitialAmount',RTX});
% 	if nargin > 2
% 		addcontent(varin,{'species','CD16','InitialAmount',CD16});
% 		addcontent(varin,{'species','CD20','InitialAmount',CD20});
% 		addcontent(varin,{'parameter','k16RTXon','Value',k16RTXon});
% 		addcontent(varin,{'parameter','k16RTXoff','Value',k16RTXoff});
% 		addcontent(varin,{'parameter','kRTX20on','Value',kRTX20on});
% 		addcontent(varin,{'parameter','kRTX20off','Value',kRTX20off});
% 	end
% 	[succ2,var2,m2] = sbiosteadystate(m1,varin);
% 	SStrimer = m2.Species(6).InitialAmount;
% else
% 	SSfun = @(RTX)SSadcc(RTX,RTXconst);
% 	SStrimer = arrayfun(SSfun,RTX);
% end


	