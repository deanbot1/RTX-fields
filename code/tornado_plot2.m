%% 'sharknado plot' of sensitivities of model at the default values for other parameters
clear all; close all; 

Tpar = readtable('adcx_parameter_results_fat.csv');  % these estimates will be the 'centerline' in the plot  ** nb "sens" column is not used!


SNPvar = 'F158'; % which SNP to be 'centered' around
CellLine = 'SUDHL4'; % which cell line to be 'centered' around
varname = sprintf('%son%s',SNPvar,CellLine);
Rconc = 1000; % RTX conc in uM?
tend = 4; % hours of ADCC assay
titlstr = sprintf('%s SNP on %s cells:%d\\muM RTX @ %dh',SNPvar,CellLine,Rconc,tend);


for j = 1:height(Tpar)
	pbest.(Tpar.name{j}) = Tpar.E1_V158onZ138(j); % pick which experiment to center tornado around
end

%% test the main function
xname = 'gamma';
dp = struct('CD20',0.01); % for example, dp = struct('CD20',0.01) sets up the question, how much fold increase in 'CD16' is needed to offset 2 logs drop in CD20?
dx = how_much_dx_to_offset_dp(pbest,dp,xname,@adcx_wrapper,Rconc,[0:.1:tend]);

%% specify the objective function on which to do sensitivity analysis
ofun = @(pstruct)how_much_dx_to_offset_dp(pstruct,dp,'gamma',@adcx_wrapper,Rconc,[0:.1:tend]) ; % set up anonymous Output function to run analysis on
xlab = 'how much gamma fc required to offset CD20 \leftarrow CD20*0.01'; % this needs to be hand curated!  Make it automatic...
ppname = fieldnames(dp); ppname = ppname{1};
xlab = sprintf('%s f.c. to offset %s \\leftarrow %s \\times %3.3g',xname,ppname,ppname,dp.(ppname));


ofun = @(pstruct)adcx_wrapper(pstruct,Rconc,[0:.1:tend]);
xlab = '%ADCC';


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
    fill([xvec flip(xvec)],j+[-yvec flip(yvec)],barcol,'edgecolor','none');
	%dk = diff(kint)/4;
	%fill(funvals(jj,find(qlevels==0.5))+dk*[-2 0 -2],0.1*[0 0 dk]+j,barcol,'edgecolor','none');
	switch(sign(xvec(end)-xvec(1)))
		case 1
			plot(funvals(jj,ceil(Nquan/2)),j,'w>');
		case -1
			plot(funvals(jj,ceil(Nquan/2)),j,'w<');
	end
	%plot(funvals(jj,ceil(Nquan/2)),j,'ro');
end

set(gca,'Ytick',[1:Nsens],'YtickLabel',sensnames(isort));
set(gca,'Color',[.65 .78 .93]);
%set(gca,'Color','c')
xlabel(xlab);
title(titlstr);
grid on
