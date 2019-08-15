function p = vect2struct(pvals,pnames)
% converts a vector of pvals to a parameter structure p with fieldnames
% given by pnames
%
% $URL$
% $Author$
% $Rev$
% $Date$

for j = 1:length(pnames)
	p.(pnames{j}) = pvals(j);
end
