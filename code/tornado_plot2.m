%% 'sharknado plot' of sensitivities of model at the default values for other parameters
clear all; close all; 

Tpar = readtable('adcx_parameter_results_fat.csv');  % these estimates will be the 'centerline' in the plot  ** nb "sens" column is not used!

p_over = struct('E0toT0',25); % OVERRIDEs to PBEST values...  

varname = 'E4_F158onSUDHL4' ; % which experiment to center the sensitivity analysis around
isnp1 = strfind(varname,'_'); isnp2 = strfind(varname,'on'); % this code won't work if there's something after the cell line name in varname!!!!
SNPvar = varname(isnp1+1:isnp2-1);
CellLine = varname(isnp2+2:end);
Rconc = 100; % RTX conc in uM
tend = 4; % hours of ADCC assay
titlstr = sprintf('%s SNP on %s cells:%d\\muM RTX @ %dh',SNPvar,CellLine,Rconc,tend);


for j = 1:height(Tpar)
	pbest.(Tpar.name{j}) = Tpar.(varname)(j); % pick which experiment to center tornado around
end

% overrride pbest values with local p_over values
if exist('p_over')
	poverstr = '';
	povernames = fieldnames(p_over);
	for k = 1:length(povernames)
		pppname = povernames{k};
		if ismember(pppname,fieldnames(pbest))
			if abs(pbest.(pppname)-p_over.(pppname)) > eps
			disp(sprintf('overriding pbest.%s=%d with %d',pppname,pbest.(pppname),p_over.(pppname)));
			pbest.(pppname) = p_over.(pppname);
			poverstr = [poverstr sprintf(' %s=%d',pppname,pbest.(pppname))];
			end
		else
			error(sprintf('no such parameter: %s',pppname));
		end
	end
	if ~isempty(poverstr)
		titlstr = {titlstr,poverstr};
	end
end

%% test the main function
xname = 'gamma';
dp = struct('CD20',0.1); % for example, dp = struct('CD20',0.01) sets up the question, how much fold increase in 'CD16' is needed to offset 2 logs drop in CD20?
dx = how_much_dx_to_offset_dp(pbest,dp,xname,@adcx_wrapper,Rconc,[0:.1:tend]);

%% specify the objective function on which to do sensitivity analysis
ofun = @(pstruct)how_much_dx_to_offset_dp(pstruct,dp,xname,@adcx_wrapper,Rconc,[0:.1:tend]) ; % set up anonymous Output function to run analysis on
ppname = fieldnames(dp); ppname = ppname{1}; % this breaks, or is wrong, if there is more than one parameter listed in dp
xlab = sprintf('%s f.c. to offset %s \\leftarrow %s \\times %3.3g',xname,ppname,ppname,dp.(ppname));

% uncomment below if you want to sensitivity analysis on %ADCC
%ofun = @(pstruct)adcx_wrapper(pstruct,Rconc,[0:.1:tend]);
%xlab = '%ADCC';


%% next step, read in quantiles and step through them for each paramter in the quantiles table.
Tfull = readtable('quantiles_bootstraps_skinny.csv'); Tfull.Properties.RowNames = Tfull.Var1;
qlevels = Tfull{'name',2:end};
Tq = readtable('quantiles_bootstraps_skinny.csv','header',1);
Snames = Tq.Var1;
% next step is to convert Tq into something the same "height" as Tpar,
% etc...

pnames = fieldnames(pbest);
qmat = NaN(length(pnames),length(qlevels));

jgood = [];
for j = 1:length(pnames)
	pname = pnames{j};
	igood = find(strcmp(Snames,pname));
	if isempty(igood) % maybe it's a SNP variant 
		igood = find(strcmp(Snames,sprintf('%s_%s',pname,SNPvar)));
		if isempty(igood) % maybe it's a cell line variant
			igood = find(strcmp(Snames,sprintf('%s_%s',pname,CellLine)));
		end
	end
	% at this point igood better not be empty! Don't forget to inverse transform!
	if isempty(igood) % then it's not on the list and should be fixed to default
		qmat(j,:) = Tpar.default(j);
	else
		jgood(end+1)=j;
	switch(Tpar.xform{j})
		case 'log'
			qmat(j,:) = exp(Tq{igood,2:end});
		case 'linear'
			qmat(j,:) = Tq{igood,2:end};
		case 'logit'
			1./(1 + exp(-Tq{igood,2:end}));
		otherwise
			error(sprintf('unsupported transformation %s',Tpar.xform{j}));
	end
	end
end

%% now we run the senstivity main loop

%jgood = find(std(qmat,[],2)./mean(qmat,2) > 100*eps);
sensnames = pnames(jgood); % only loop over parameters that have variability over qmat
qmat = qmat(jgood,:);
Nsens = length(sensnames);
Nquan = length(qlevels);

funvals = NaN(Nsens,Nquan); % matrix of 'objective' function values

for j = 1:Nsens
	senspar = sensnames{j}; % which parameter we're perturbing
	for k = 1:Nquan
		psens = pbest;
		psens.(senspar) = qmat(j,k);
		funvals(j,k) = ofun(psens);
	end
end

%% now we plot the results
figure('Position',[560   153   407   424]);
fbest = ofun(pbest);
% pmid = pbest;
% for j = 1:Nsens
% 	pmid.(sensnames{j}) = funvals(j,find(qlevels==0.5));
% end
% fmid = ofun(pmid);
plot(fbest*[1 1],[0 Nsens+1],'w-'); hold on

spans = max(funvals,[],2)-min(funvals,[],2);
[~,isort] = sort(spans);

barcol = 0.5*[1 1 1];
dj = 0.05; % delta j for shark thickness 
for j = 1:Nsens
	jj = isort(j);
	
	xvec = funvals(jj,:);
	yvec = dj*[1:1:floor(Nquan/2) ceil(Nquan/2) floor(Nquan/2):-1:1]; 
% 	for k = 1:floor(Nquan/2)
% 		kint = funvals(jj,[k,Nquan-k]);
% 		plot(kint,j*[1 1],'-','LineWidth',2*k,'Color',barcol);
% 	end
    fill([xvec flip(xvec)],j+[-yvec flip(yvec)],barcol,'edgecolor','none'); % shark 'body'
	dk = diff(xvec(ceil(Nquan/2)+[-1 1]))/8;
	fill(funvals(jj,find(qlevels==0.5))+dk*[-1 1 -1],0.2*[0 0 1]+j+dj*ceil(Nquan/2),barcol,'edgecolor',barcol); % shark dorsal 'fin'
	fill(funvals(jj,find(qlevels==0.5))+dk*[-1 1 -1],-0.1*[0 0 1]+j,barcol,'edgecolor','w'); % shark side 'fin'
	plot(funvals(jj,end-1),j+dj,'w.'); % eye
	fill(funvals(jj,1)+dk*[0 0 2],dj*ceil(Nquan/2)*[1 -1 0]+j,barcol,'edgecolor',barcol); % tailfin
% 	switch sign(dk)
% 		case 1
% 			plot(funvals(jj,1),j,'>','Color',barcol,'MarkerFaceColor',barcol);
% 		case -1
% 			plot(funvals(jj,1),j,'<','Color',barcol,'MarkerFaceColor',barcol);
% 	end
% 	switch(sign(xvec(end)-xvec(1)))
% 		case 1
% 			plot(funvals(jj,ceil(Nquan/2)),j,'w>');
% 		case -1
% 			plot(funvals(jj,ceil(Nquan/2)),j,'w<');
% 	end
% 	
	hfirst = text(xvec(1)-dk/2,j,num2str(qmat(jj,1),'%2.3g'),'FontSize',8,'HorizontalAlignment','right','Rotation',0);
	hlast  = text(xvec(end)+dk/2,j,num2str(qmat(jj,end),'%2.3g'),'FontSize',8,'HorizontalAlignment','left','Rotation',0);
	hmid = text(xvec(ceil(Nquan/2)),j-dj*ceil(Nquan/2),num2str(qmat(jj,ceil(Nquan/2)),'%2.3g'),'FontSize',8,'HorizontalAlignment','center','Rotation',0,'VerticalAlignment','top');
	
% 	if xvec(end) > xvec(1)
% 		set(hlast,'VerticalAlignment','top');
% 		set(hfirst,'VerticalAlignment','bottom');
% 	else
% 		set(hlast,'VerticalAlignment','bottom');
% 		set(hfirst,'VerticalAlignment','top');
% 	end

if xvec(end) < xvec(1)
	set(hlast,'HorizontalAlignment','right');
	set(hfirst,'HorizontalAlignment','left');
end
	%plot(funvals(jj,ceil(Nquan/2)),j,'ro');
end

set(gca,'Ytick',[1:Nsens],'YtickLabel',sensnames(isort));
set(gca,'Color',[.65 .78 .93]);
%set(gca,'Color','c')
xlabel(xlab);
title(titlstr);
grid on
