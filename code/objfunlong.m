function ofun = objfunlong(p,expt,pxform,errfunlong, pvec0, cvec, bweight, ii)

Ne = length(expt);
funval = 0;
ttnames = fieldnames(expt(1).pmap); % names of target paramter names
pbig = pvec2struct(p,pxform);
funvec = [];

for j = 1:Ne
	for k = 1:length(ttnames)
		tname = ttnames{k};
		pfunc = expt(j).pmap.(tname);
		pstruct.(tname) = pfunc(pbig);
	end
	Ypred = expt(j).model(pstruct,expt(j).xval);
	if isfield(expt,'err')  % if measurement error vector exists, use it! That means errfunlong had better know what to do with it!
		funvec = vertcat(funvec,errfunlong(Ypred,expt(j).obs,expt(j).err));
	else % backwards compatibility
		funvec = vertcat(funvec, errfunlong(Ypred,expt(j).obs));
	end

end

%B = bweight*sum(((p(ii)-pvec0(ii))./cvec(ii)).^2); bayes penalty, scalar
if size(cvec,2) == 1
	B = bweight*(p(ii)-pvec0(ii))./cvec(ii); % bayes penalty, vector
else  % cvec is a COVARIANCE matrix not a STDEV vector
	B = bweight*chol(inv(cvec(ii,ii)))*(p(ii)-pvec0(ii));
end
	
ofun = vertcat(funvec, B);