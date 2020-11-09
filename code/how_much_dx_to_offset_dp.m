function [dx,do] = how_much_dx_to_offset_dp(pin,dp,xname,dx0,modelfun,varargin) 
% a function that returns how much fold change in a particular parameter
% named xname is required to offset a perturbation dp to the input
% parameter vector pin. Optionally returns do, which is change in ADCC% caused by dp relative to pin (fraction or fold?) 


pdelta = pin;
pnames = fieldnames(dp);
dldp = 0;
for j = 1:length(pnames)
	pdelta.(pnames{j}) = pdelta.(pnames{j})*dp.(pnames{j});
	dldp = max(dldp,abs(log(dp.(pnames{j}))));
end

dldp = -1; % override the derived dldp

adcx_target = modelfun(pin,varargin{:});
adcx_init = modelfun(pdelta,varargin{:});
do = adcx_init/adcx_target;

if isempty(dx0)
dlx = log(adcx_target/adcx_init); % how many logs to perturb over to see what 'slope' is
dldx = dldp; % perturb x on same scale as p

%adcx_hi = modelfun(pnewfun( dlx,xname,pdelta),varargin{:});
%adcx_lo = modelfun(pnewfun(-dlx,xname,pdelta),varargin{:});

%dladcx_dldx = (log(adcx_hi) - log(adcx_lo))/(2*dlx);
%ldx0 = (log(adcx_target)-log(adcx_init))/dladcx_dldx;

df_dlp = (modelfun(pnewfun(dldp,pnames{1},pin),varargin{:})-modelfun(pnewfun(-dldp,pnames{1},pin),varargin{:}))/(2*dldp); %center difference partial of f wrt logp
df_dlx = (modelfun(pnewfun(dldx,xname,pin),varargin{:})-modelfun(pnewfun(-dldx,xname,pin),varargin{:}))/(2*dldx); % center difference partial of f wrt logx

ldx0 = -log(pdelta.(pnames{1}))*df_dlp/df_dlx;

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