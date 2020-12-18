function [S,AX,BigAx,H,HAx] = plotparammatrix(pboot,pxform)
% makes an annotated matrix plot
% pboot is a Nparams x Nruns matrix of parameter values
% paramnames is a corresponding list of Nparams parameter names

figure;
[S,AX,BigAx,H,HAx] = plotmatrix(pboot');
nparams = size(pboot,1);
paramnames = fieldnames(pxform);
for j = 1:nparams
	ylabel(AX(j,1),paramnames{j},'Interpreter','none')
	xlabel(AX(nparams,j),paramnames{j},'Interpreter','none')
end
hold on

title(BigAx,'raw values')
% for j = 1:nparams
% 	ylabel(AX(j,j),paramnames{j})
% 	xlabel(AX(nparams,j),paramnames{j})
% 	switch pxform.(paramnames{j})
% 		case 'linear'
% 			ixform = @(x)x;
% 		case 'log'
% 			ixform = @(x)exp(x);
% 		otherwise
% 			error('not supported yet');
% 	end
% 	set(AX(nparams,j),'XTickLabel',num2str(ixform(get(AX(nparams,j),'Xtick'))','%2.2g'));
% 	set(AX(j,1),'Yticklabel',num2str(ixform(get(AX(j,1),'Ytick'))','%2.2g'));
% end

covmat = cov(pboot');
cax = max(max(abs(covmat)))*[-1 1];
colormap cool
caxis(cax);
colorbar;
cmap = colormap('cool');
for j = 1:nparams
	for k = 1:nparams
		col = cmap(ceil(64*(0.5+(1/diff(caxis))*covmat(j,k))),:);
		if j==k
			set(HAx(j),'Color',col);
		else
			set(AX(j,k),'Color',col);
		end
	end
end

set(S,'Color','w');