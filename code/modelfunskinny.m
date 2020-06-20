function Ymod = modelfunskinny(p,expt,pxform)
% This function returns the model predicted values for the [RTX]
% concentrations in each of the 5 experiments as a single vector

Ne = length(expt);
Ymod = [];
ttnames = fieldnames(expt(1).pmap); % names of target paramter names
pbig = pvec2struct(p,pxform);

for j = 1:Ne
	for k = 1:length(ttnames)
		tname = ttnames{k};
		pfunc = expt(j).pmap.(tname);
		pstruct.(tname) = pfunc(pbig);
	end
	Ypred = expt(j).model(pstruct,expt(j).xval);
    Ymod = vertcat(Ymod, Ypred);
end
end