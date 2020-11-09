function [dx,do] = how_much_dx_to_offset_dp(pin,dp,xname,dx0,modelfun,varargin) 
% a function that returns how much fold change in a particular parameter
% named xname is required to offset a perturbation dp to the input
% parameter vector pin. Optionally returns do, which is change in ADCC% caused by dp relative to pin (fraction or fold?) 


pdelta = pin;
pnames = fieldnames(dp);
maxabsldp = 0;
for j = 1:length(pnames)
	pdelta.(pnames{j}) = pdelta.(pnames{j})*dp.(pnames{j});
	maxabsldp = max(maxabsldp,abs(log(dp.(pnames{j}))));
end


adcx_target = modelfun(pin,varargin{:});
adcx_init = modelfun(pdelta,varargin{:});
do = adcx_init/adcx_target;

if isempty(dx0)
dlx = log(adcx_target/adcx_init); % how many logs to perturb over to see what 'slope' is
dlx = maxabsldp; % perturb x on same scale as p
adcx_hi = modelfun(pnewfun( dlx,xname,pdelta),varargin{:});
adcx_lo = modelfun(pnewfun(-dlx,xname,pdelta),varargin{:});

dladcx_dldx = (log(adcx_hi) - log(adcx_lo))/(2*dlx);
ldx0 = (log(adcx_target)-log(adcx_init))/dladcx_dldx;

% dfdp = (modelfun(pnewfun(log(pdelta.(pnames{1})),pnames{1},pin),varargin{:}) - modelfun(pin,varargin{:}))/log(pdelta.(pnames{1})); % partial f wrt p
% dfdx = (modelfun(pnewfun(log(pdelta.(pnames{1})),xname,pin),varargin{:}) - modelfun(pin,varargin{:}))/log(pdelta.(pnames{1})); % partial f wrt x
% 
% gfp = [-dfdp;dfdx];
% ldx0 = -dfdp*log(pdelta.(pnames{1}))

else
	ldx0 = log(dx0);
end

[ldx,fval,exitflag,output] = fzero(@(lddx)adcx_target - modelfun(pnewfun(lddx,xname,pdelta),varargin{:}),ldx0);
dx = exp(ldx);
end

function pnew = pnewfun(ldx,xname,pdelta)
pnew = pdelta;
pnew.(xname) = pnew.(xname)*exp(ldx);
end