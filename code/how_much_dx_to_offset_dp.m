function [dx,do] = how_much_dx_to_offset_dp(pin,dp,xname,modelfun,varargin) 
% a function that returns how much fold change in a particular parameter
% named xname is required to offset a perturbation dp to the input
% parameter vector pin. Optionally returns do, which is change in ADCC% caused by dp relative to pin (fraction or fold?) 


pdelta = pin;
pnames = fieldnames(dp);
for j = 1:length(pnames)
	pdelta.(pnames{j}) = pdelta.(pnames{j})*dp.(pnames{j});
end


adcx_target = modelfun(pin,varargin{:});
adcx_init = modelfun(pdelta,varargin{:});
do = adcx_init/adcx_target;

dlx = log(adcx_target/adcx_init); % how many logs to perturb over to see what 'slope' is
adcx_hi = modelfun(pnewfun( dlx,xname,pdelta),varargin{:});
adcx_lo = modelfun(pnewfun(-dlx,xname,pdelta),varargin{:});

dladcx_dldx = (log(adcx_hi) - log(adcx_lo))/(2*dlx);
ldx0 = (log(adcx_target)-log(adcx_init))/dladcx_dldx;

[ldx,fval,exitflag,output] = fzero(@(lddx)adcx_target - modelfun(pnewfun(lddx,xname,pdelta),varargin{:}),ldx0);
dx = exp(ldx);
end

function pnew = pnewfun(ldx,xname,pdelta)
pnew = pdelta;
pnew.(xname) = pnew.(xname)*exp(ldx);
end