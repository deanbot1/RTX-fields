function Tskinny = save_skinny_results(pxform,cvs,pinit,pbest,filename)

clear name xform cv cvxf initial final
[name,is] = sort(fieldnames(pbest));
for k = 1:length(name)
	xform{k,1} = pxform.(name{k});
	cv(k,1) = cvs.(name{k});
%	cvxf(k,1) = cvec(is(k));
	initial(k,1) = pinit.(name{k});
	final(k,1) = pbest.(name{k});
end

Tskinny = table(name,xform,cv,initial,final);
disp(Tskinny)
writetable(Tskinny,filename); % save to disk
disp(sprintf('%s written to: %s','Tskinny',filename));
