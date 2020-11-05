function h = handiplot(Tout)
% automatically plots a handiwrap output table Tout

varnames = Tout.Properties.VariableNames;

yname = varnames{1}; % name of y variable is always the first
constnames = {};
xname = ''; % name of x variable
lname = ''; % plot curve label names
sname = ''; % subplotname

for j = 2:length(varnames)
	if std(Tout.(varnames{j})) < eps, constnames{end+1} = varnames{j},continue;end
	if isempty(xname),xname = varnames{j},continue;end
	if isempty(lname),lname = varnames{j},continue;end
	if isempty(sname),sname = varnames{j},continue;end
end

if ~isempty(sname)
	Us = unique(Tout.(sname));
	Ns = length(Us);
else
	Ns = 1;
	Us = [];
end

for js = 1:Ns
	subplot(ceil(sqrt(Ns)),ceil(Ns/ceil(sqrt(Ns))),js);
	if ~isempty(sname)
		Ts = Tout(Tout.(sname)==Us(js),:);
	else 
		Ts = Tout;
	end
	if ~isempty(lname)
		Ul = unique(Tout.(lname));
		Nl = length(Ul);
	else
		Nl = 1;
		Ul = [];
	end
	for jl = 1:Nl
		if ~isempty(lname)
			Tsl = Ts(Ts.(lname)==Ul(jl),:);
		else
			Tsl = Ts;
		end
		hl = plot(Tsl.(xname),Tsl.(yname),'-o'); hold on
	
		if ~isempty(lname)
			text(Tsl{end,xname},Tsl{end,yname},sprintf('  %s=%3.3g',lname,Ul(jl)),'Color',get(hl,'Color'));
		end
	end
	%Ux = unique(Ts.(xname));
	%if length(Ux) < 15, set(gca,'Xtick',Ux);end
	xlabel(xname,'Interpreter','none')
	ylabel(yname,'Interpreter','none')
	set(gca,'Xlim',[0.9*min(Ts.(xname)),1.1*max(Ts.(xname))]);
	set(gca,'Ylim',[0.9*min(Tout.(yname)),1.1*max(Tout.(yname))]);
	if ~isempty(constnames)
		title(sprintf('%s=%3.3g,%s=%3.3g',sname,Ts.(sname)(1),constnames{:},Tout{1,constnames}));
	else
		title(sprintf('%s=%3.3g,%s=%3.3g',sname,Ts.(sname)(1)));
	end
	grid on
end