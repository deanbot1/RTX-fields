function [pmap] = default_pmap(par)
% makes a structure pmap whose fields are parameter names known to the
% model and whose values are functions that 

pnames = par.paramnames;
tnames = unique(par.target,'stable');
for k = 1:length(tnames)
	ibig = find(strcmp(par.target,tnames{k}));
	if length(ibig) == 1
		switch par.fit(ibig)
			case 0
				pmap.(tnames{k}) = @(pfit)par.value(ibig);
			case 1
				pmap.(tnames{k}) = @(pfit)pfit.(tnames{k});
			case -1
				error('-1 value for fit only allowed if one target has multiple paramnames associated with it');
			otherwise
				error('unsupported fit value');
		end
	else
		if max(par.fit(ibig)) ~= min(par.fit(ibig)) % they have to be all equal
			warning('when multiple parameters are mapped to one parameter, fit entries must equal each other else mapping is ambiguous -- actually it might not be let me think about it. I guess you could have zeros and -1s mixed together, the zeros are pass throughs, the -1s are fit');
		end
		pmap.(tnames{k}) = @(pfit)NaN; % all these mappings are determined at the experiment level, there is no default behavior
		switch par.fit(ibig(1))  
 			case 0 % they're going to get passed through by experiment -- do nothing
 			case 1 % not sure what should happen here
				pmap.(tnames{k}) = @(pfit)pfit.(tnames{k});
				%error('not sure what to do when many-->1 mapping and fit==1');
 			case -1 % they get fit depending on experiment number
			otherwise
				error('unsupported fit value');
  		end
	end
end