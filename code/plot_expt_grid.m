function plot_expt_grid(expt,pbest,xdim,xspan,varargin)
% plots experiments/model in a grid
% the columns are the experiments in expt
% the rows are different unique expt(j).xval(:,~xdim) values
% expt is a structure array indexed by experiment number
% pbest is a 'tall' parameter vector, maybe?
% xspan is requested simulation vector for first X variable (column xdim in
% expt(j).xval
% varargin is comma separated list of name,value pairs passed as graphics
% options, ie, 'linewidth','2', for the model curve.

Ne = length(expt);
ttnames = fieldnames(expt(1).pmap); % names of target paramter names

xvalbig = cat(1,expt.xval);
fdim = 3-xdim; % dimension along which to 'facet' plots
ufac = unique(xvalbig(:,fdim));
nfac = length(ufac); % number of vertical 'facets'
subplot2d = @(nrow,ncol,rowi,colj)subplot(nrow,ncol,ncol*(rowi-1)+mod(colj-1,ncol)+1); % 2d subplot


for j = 1:Ne
	for k = 1:nfac
		subplot2d(nfac,Ne,k,j);
		igood = find(expt(j).xval(:,fdim)==ufac(k));
		if ~isempty(igood)
			
			errorbar(expt(j).xval(igood,xdim),expt(j).obs(igood),expt(j).err(igood),'ro','markerfacecolor','r'); hold on
			
			for kk = 1:length(ttnames)
				tname = ttnames{kk};
				pfunc = expt(j).pmap.(tname);
				pstruct.(tname) = pfunc(pbest);
			end
		end
			xmat = zeros(length(xspan),2);
			xmat(:,xdim) = xspan;
			xmat(:,fdim) = ufac(k);
			Ypred = expt(j).model(pstruct,xmat);
			plot(xspan,Ypred,'b.-');
			set(gca,varargin{:});
			grid on
			%title(num2str(max(Ypred)));
		
		if j==1
			ylabel(sprintf('%s,%s=%3.3g',char(expt(j).Ynames),expt(j).Xname{fdim},ufac(k)));
		end
		if k==1
			title(expt(j).name,'Interpreter','none');
		end
	end
	xlabel(expt(j).Xname{xdim});
end