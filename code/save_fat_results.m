function Tout = save_fat_results(Tpar,pbest,expt,filename)

Tout = Tpar;
Npar = length(Tout.name);
for j = 1:length(expt)
	clear ppp
	for k = 1:Npar
		pname = Tout.name{k};
		ppp(k,1) = feval(expt(j).pmap.(pname),pbest);
	end
	Tout.(expt(j).name) = ppp;
end

disp(Tout)
writetable(Tout,filename); % save to disk
disp(sprintf('%s written to: %s','Tfat',filename));