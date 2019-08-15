function obj = objfun(p,expt,p0,prow,pcol,ixf,errfun)
% OBJFUN calculates the objective function value by summing the errfun value over
% all experiments.
%
% INPUTS
% P is a big ??x1 parameter vector that the optimizer likes
% EXPT is the experiment structure array (see help on plot_expt for more
% detail) with .model and .obs fields
% P0 is the initial parameter guess, Npx1
% PROW, PCOL are bookkeeping functions described in pfluff and psquash.m
% IXF is the inverse transform local function defined in glocal_demo.m
% which converts parameters on the real line into their 'natural' parameter
% space
% ERRFUN is a function of Ypred and Yobs matrices that calculates say sum
% squared error between the two matrices in log space
% 
% OUTPUTS
% obj is the objective function value 
%
% $URL$
% $Author$
% $Rev$
% $Date$

Ne = length(expt);

obj = 0;
pmat = ixf(pfluff(p,p0,prow,pcol,Ne));
for i = 1:Ne
	Ypred = expt(i).model(pmat(:,i),expt(i).time);
	obj = obj + errfun(Ypred,expt(i).obs);
end