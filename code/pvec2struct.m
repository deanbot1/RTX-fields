function pstruct = pvec2struct(pvec,pxform)

fnames = fieldnames(pxform);
for j = 1:length(pvec)
	fname = fnames{j};
	switch pxform.(fname)
		case 'log'
			pstruct.(fname) = exp(pvec(j));
		case 'logit'
			pstruct.(fname) = 1/(1 + exp(-pvec(j))); % DOUBLE CHECK THIS IS INVERSE LOGIT@
		case 'linear'
			pstruct.(fname) = pvec(j);
	end
end