function pvec = pstruct2vec(pstruct,pxform)
% converts parameter structure to a column vector applying transforms as
% specified by pxform (also used for indexing)
fnames = fieldnames(pxform);
pvec = NaN(length(fnames),1);
for j = 1:length(fnames)
	fname = fnames{j};
	switch pxform.(fname)
		case 'log'
			pvec(j) = log(pstruct.(fname));
		case 'logit'
			pvec(j) = log(pstruct.(fname)/(1-pstruct.(fname)));
		case 'linear'
			pvec(j) = pstruct.(fname);
	end
end