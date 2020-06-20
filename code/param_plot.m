function param_plot(Tskinny,titlestr)

xlim = 10.^[-1 1];
Nfit = height(Tskinny);
semilogx(Tskinny.final./Tskinny.initial,1:Nfit,'ro','MarkerFaceColor','r'); hold on
plot(ones(1,Nfit),1:Nfit,'ko');
hold on
set(gca,'YDir','reverse','YTick',1:Nfit,'YTickLabel',Tskinny.name);
plot([1 1],[0.5 Nfit+0.5],'k:');
for j = 1:Nfit
	if Tskinny.cv(j) < Inf
		plot([1+Tskinny.cv(j) 1./(1+Tskinny.cv(j))],[j j],'k-');
		xlim(1) = min([xlim(1),get(gca,'Xlim')]);
		xlim(2) = max([xlim(2),get(gca,'Xlim')]);
	else
		xlim = get(gca,'Xlim');
		plot(xlim,[j j],'k-');
		plot(min(xlim),j,'k<','MarkerFaceColor','k');
		plot(max(xlim),j,'k>','MarkerFaceColor','k');
	end
	text(Tskinny.final(j)/Tskinny.initial(j),j,sprintf('%3.3g',Tskinny.final(j)),'Color','r','VerticalAlignment','Bottom','HorizontalAlignment','Center');
	text(1,j,sprintf('%3.3g',Tskinny.initial(j)),'Color','k','VerticalAlignment','Top','HorizontalAlignment','Center');
end
set(gca,'TickLabelInterpreter','none','Xlim',xlim)
xlabel('Final Estimates (red) / Initial Estimates');
title(titlestr);