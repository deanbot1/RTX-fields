function cvec = cvstruct2vec(cvstruct,pstruct,pxform)
% converts parameter structure to a column vector applying transforms as
% specified by pxform (also used for indexing)
fnames = fieldnames(pxform);
cvec = NaN(length(fnames),1);
for j = 1:length(fnames)
	fname = fnames{j};
	switch pxform.(fname)
		case 'log'
			cvec(j) = sqrt(log((cvstruct.(fname))^2 + 1)); % since CV of lognormal = sqrt(exp(sigma^2)-1), inverse is this
		case 'logit'
			cvec(j) = cvstruct.(fname); % this is arbitrary as there's no analytic expression for variance of logit normal distribution
		case 'linear'
			cvec(j) = cvstruct.(fname)*abs(pstruct.(fname)); % since CV or normal = std/mu, we need pstruct to provide scale
	end
end