function dx = how_much_dx_to_offset_dp(pin,dp,xname,modelfun,varargin) 
% a function that returns how much fold change in a particular parameter
% named xname is required to offset a perturbation dp to the input
% parameter vector pin. 


pdelta = pin;
pnames = fieldnames(dp);
for j = 1:length(pnames)
	pdelta.(pnames{j}) = pdelta.(pnames{j})*dp.(pnames{j});
end


adcx_target = modelfun(pin,varargin{:});
adcx_init = modelfun(pdelta,varargin{:}); % not used, just warming up...


[dx,fval,exitflag,output] = fzero(@(ddx)adcx_target - modelfun(pnewfun(ddx,xname,pdelta),varargin{:}),1);
end

function pnew = pnewfun(dddx,xname,pdelta)
pnew = pdelta;
pnew.(xname) = pnew.(xname)*dddx;
end