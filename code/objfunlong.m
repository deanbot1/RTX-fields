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
	funvec = vertcat(funvec, errfunlong(Ypred,expt(j).obs));

end

B = bweight*sum(((p(ii)-pvec0(ii))./cvec(ii)).^2); % +bweight*profiled parameter stuff
ofun = vertcat(funvec, B);