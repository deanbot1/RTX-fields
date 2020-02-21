function plot_expt2(expt,pbest,xspan,varargin);

Ne = length(expt);
ttnames = fieldnames(expt(1).pmap); % names of target paramter names
for j = 1:Ne
	subplot(1,Ne,j);
	plot(expt(j).xval,expt(j).obs,'ro','markerfacecolor','r'); hold on;
		for k = 1:length(ttnames)
		tname = ttnames{k};
		pfunc = expt(j).pmap.(tname);
		pstruct.(tname) = pfunc(pbest);
	end
	Ypred = expt(j).model(pstruct,xspan);
	plot(xspan,Ypred,'b.-'); 
	set(gca,varargin{:});
	title(expt(j).name,'Interpreter','none');
	xlabel(expt(j).Xname);
	grid on
	if j==1
		ylabel(expt(j).Ynames);
	end
end