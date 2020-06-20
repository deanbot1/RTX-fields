function funval = objfun2(p,expt,pxform,errfun)

Ne = length(expt);
funval = 0;
ttnames = fieldnames(expt(1).pmap); % names of target paramter names
pbig = pvec2struct(p,pxform);

for j = 1:Ne
	for k = 1:length(ttnames)
		tname = ttnames{k};
		pfunc = expt(j).pmap.(tname);
		pstruct.(tname) = pfunc(pbig);
	end
	Ypred = expt(j).model(pstruct,expt(j).xval);
	funval = funval + errfun(Ypred,expt(j).obs);
end

