	function SStrimer = SSadcc(RTX,RTXconst,CD16,CD20,k16RTXon,k16RTXoff,kRTX20on,kRTX20off)

sbioloadproject('ADCC_reactions.sbproj');
m1.Species(1).ConstantAmount = RTXconst;
if length(RTX)==1
	varin = sbiovariant('v1',{'species','RTX','InitialAmount',RTX});
	if nargin > 2
		addcontent(varin,{'species','CD16','InitialAmount',CD16});
		addcontent(varin,{'species','CD20','InitialAmount',CD20});
		addcontent(varin,{'parameter','k16RTXon','Value',k16RTXon});
		addcontent(varin,{'parameter','k16RTXoff','Value',k16RTXoff});
		addcontent(varin,{'parameter','kRTX20on','Value',kRTX20on});
		addcontent(varin,{'parameter','kRTX20off','Value',kRTX20off});
	end
	[succ2,var2,m2] = sbiosteadystate(m1,varin);
	SStrimer = m2.Species(6).InitialAmount;
else
	SSfun = @(RTX)SSadcc(RTX,RTXconst);
	SStrimer = arrayfun(SSfun,RTX);
end


	