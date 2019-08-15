function px = fxform(p,xforms)
% calculates the forward transform on parameter vector or matrix p based on user specified transformation
% provided in xforms Npx1 cell array of strings
% currently 'log' and 'logit' are supported.
% anything else is ignored (identity). 
%
% $URL$
% $Author$
% $Rev$
% $Date$

ilog = find(strcmpi('log',xforms));
ilogit = find(strcmpi('logit',xforms));

px = p;
px(ilog,:) = log(p(ilog,:));
px(ilogit,:) = log(p(ilogit,:)./(1-p(ilogit,:)));

