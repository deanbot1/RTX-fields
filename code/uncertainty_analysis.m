%% uncertainty analysis measures spread of key output functions over the bootstrapped uncertainty in parameters

clear all; close all; 

foo = load('./bootstrap_500runs_lsqn_8pars_r1.mat');
pbigbootall = foo.pbigbootall;
pbootnames = fieldnames(foo.pinit); % VERONICA IS THIS THE RIGHT WAY TO GET NAMES LIGNED UP WITH ROWS OF PBIGBOOTALL?
pxform = foo.pxform;
expt = foo.expt;
[Npar,Nboot] = size(pbigbootall);

%% run the analysis
Tpar = readtable('adcx_parameter_results_fat.csv');  % VERONICA ARE THESE THE RIGHT DEFAULT PARAMETER VALUES FROM PBEST?

varname = 'E4_F158onSUDHL4' ; % which experiment to center the sensitivity analysis around
isnp1 = strfind(varname,'_'); isnp2 = strfind(varname,'on'); % this code won't work if there's something after the cell line name in varname!!!!
SNPvar = varname(isnp1+1:isnp2-1);
CellLine = varname(isnp2+2:end);
Rconc = 10; % RTX conc in uM
tend = 4; % hours of ADCC assay
titlstr = sprintf('%s SNP on %s cells:%d\\muM RTX @ %dh',SNPvar,CellLine,Rconc,tend);

% first read in teh pbest values as defaults
for j = 1:height(Tpar)
	pbest.(Tpar.name{j}) = Tpar.(varname)(j); 
end

[enames{1:length(expt)}] = deal(expt.name);
enum = find(strcmp(enames,varname));

%% test the main function
xname = 'gamma';
dp = struct('CD20',0.01); % for example, dp = struct('CD20',0.01) sets up the question, how much fold increase in 'CD16' is needed to offset 2 logs drop in CD20?
%dx = how_much_dx_to_offset_dp(pbest,dp,xname,@adcx_wrapper,Rconc,[0:.1:tend]);
%[dx,do] = how_much_dx_to_offset_dp(pbest,dp,xname,@handiwrap,[]);

%% specify the objective function on which to do sensitivity analysis
ofun = @(pstruct)how_much_dx_to_offset_dp(pstruct,dp,xname,[],@handiwrap,[]) ; % set up anonymous Output function to run analysis on
%ppname = fieldnames(dp); ppname = ppname{1}; % this breaks, or is wrong, if there is more than one parameter listed in dp
%xlab = sprintf('%s f.c. to offset %3.3g\\times%s \\rightarrow %3.3g\\timesADCC',xname,dp.(ppname),ppname,NaN);

% uncomment below if you want to sensitivity analysis on %ADCC
%ofun = @(pstruct)handiwrap(pstruct,[]);
%xlab = '%ADCC';


%% main loop

xnames = {'CD16','kon16','gamma'};
dpnames = {'CD20','RTX','kon16'};
dpfcs = [0.1, 0.1, 0.1];

ovec = {};
avec = {};



for jp = 1:length(dpnames)
	for jx = 1:length(xnames)
		
		hb = waitbar(0,'working...');
		
		ovec{jp,jx} = zeros(1,Nboot);
		avec{jp,jx} = zeros(1,Nboot);
		xname = xnames{jx}
		dp = struct(dpnames{jp},dpfcs(jp)) % for example, dp = struct('CD20',0.01) sets up the question, how much fold increase in 'CD16' is needed to offset 2 logs drop in CD20?\

		bigofun = @(pstruct,dx0)how_much_dx_to_offset_dp(pstruct,dp,xname,dx0,@handiwrap,[])
		
		xbvec = median(pbigbootall,2); % column vector, in transformed space
			ptest = pvec2struct(xbvec,pxform);
			ttnames = fieldnames(expt(enum).pmap);
			for k = 1:length(ttnames)
				tname = ttnames{k};
				pfunc = expt(enum).pmap.(tname);
				pstruct.(tname) = pfunc(ptest);
			end		
			
			pstruct.RTX = Rconc;
			
			dxmed = bigofun(pstruct,[])
			
		ofun = @(p)bigofun(p,dxmed); % set up anonymous Output function to run analysis on
		

		
		
		for jb = 1:Nboot
			xbvec = pbigbootall(:,jb); % column vector, in transformed space
			ptest = pvec2struct(xbvec,pxform);
			ttnames = fieldnames(expt(enum).pmap);
			for k = 1:length(ttnames)
				tname = ttnames{k};
				pfunc = expt(enum).pmap.(tname);
				pstruct.(tname) = pfunc(ptest);
			end
			pstruct.RTX = Rconc;
			[ovec{jp,jx}(jb),avec{jp,jx}(jb)] = ofun(pstruct);
			waitbar(jb/Nboot,hb);
			
		end
		
		close(hb)
		
	end
end

%% save the results
save('../out/uncertainty_analyses.mat','ovec','avec','pbigbootall','xnames','dpnames','dpfcs','varname','expt','pxform','pbootnames','Rconc');

%% report the results

