function pfull = pfluff(pbig,pdef,prow,pcol,Ns)
% fluffs up pbig into pfull matrix
% pbig is the estimated parameter vector
% pdef is the full size vector of default parameter vector
% prow points from pbig into the rows of pfull
% pcol points from pbig inot the columns of pfull where 0 means all columns
%
% $URL$
% $Author$
% $Rev$
% $Date$

pfull = pdef*ones(1,Ns); % pre-fluff with default parameter values
for j = 1:length(pbig)
    if pcol(j)==0
        pfull(prow(j),:) = pbig(j);
    else
        pfull(prow(j),pcol(j)) = pbig(j);
    end
end

