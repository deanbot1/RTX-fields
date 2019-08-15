function [pbig,prow,pcol] = psquash(pdef,pfit,Ns)
% squashes a pmat or a single parameter vector pdef into one long vector suitable for parameter optimization
% 
% INPUTS
% PDEF is either a single column vector of Np parameters ('default' values
% of parameters, for the purpose of initial guess or unfitted parameter
% values that wouldn't appear in the pbig output.) if PDEF is a Np x Ns
% "pmat" matrix, then psquash distributes the "local" parameter values across
% columns of PDEF appropriately.
% PFIT is a vector of 0,1,-1 corresponding to
%  0 = don't fit this parameter. 
%  1 = fit this parameter "globally" ie one value across all experiments
% -1 = fit this parameter "locally" ie one distinct value for each experiment
% NS is the number of experiments or studies
%
% OUTPUTS
% PBIG is the optimizer-friendly parameter vector
% PROW points from entries of pbig to rows of pmat, ie 1:Np where Np is the
% number of parameters in the original table (see demo.m). 
% PCOL points to which experiment (column in pmat) each entry in pbig
% corresponds to. These are bookkeeping vectors used by pfluff to perform the inverse task to psquash. 
% Don't overthink it. 
%
% $URL$
% $Author$
% $Rev$
% $Date$

i = 0;
Ne = size(pdef,2);
for k = 1:size(pdef,1)
	switch pfit(k)
		case 0 % it doesn’t even go into pbig because it’s not being fit, so optimizer shoudn’t even know about it
		case 1 % one value to rule them all
			i = i+1;
			pbig(i,1)=pdef(k);
			prow(i,1)=k;
			pcol(i,1)=0;
		case -1 % one value for each experiment
			for j = 1:Ns
				i=i+1;
				if Ne == 1
					pbig(i,1)=pdef(k);
				else
					pbig(i,1)=pdef(k,j);
				end
				prow(i,1)=k;
				pcol(i,1)=j;
			end
	end
end
