function [expt,pinit,pxform,cvs] = parse_par_expt(Tpar,Texp)

vnames = Tpar.Properties.VariableNames;
EE = arrayfun(@(s)sscanf(s{1},'E%d_*'),vnames,'UniformOutput',false);
EI = find(arrayfun(@(s)~isempty(s{1}),EE));
Nexp = length(EI);

for j = 1:Nexp
	i = EI(j);
	ename = vnames{i};
	expt(j).name = ename;
	deffun = @(rname)ifelseval(isempty(Texp{rname,ename}{1}),Texp{rname,'default'},Texp{rname,ename});
	expt(j).Xname = deffun('XNAME');
	expt(j).Ynames = {deffun('YNAME')};
	datafile = fullfile(deffun('DATAPATH'),deffun('DATAFILE'));
	dat = readtable(datafile);
	expt(j).xval = dat.(deffun('XHEADER'));
	expt(j).obs = dat.(deffun('YHEADER'));
	expt(j).model = str2func(deffun('MODEL'));
	expt(j).weight = str2num(deffun('WEIGHT'));
end

% now build par in the older way and define defpmap? No, may want to bypass
% old way, which involves parsing strings...or maybe yes... ugh
clear pxform
clear pinit cvs % initial parameters and CV structure

for j = 1:Nexp
	clear pmap
	ej = expt(j);
	for i = 1:height(Tpar)
		pname = Tpar.name{i};
		if Tpar.cv(i) < eps
			% it's fixed, no fitting
			pmap.(pname) = @(pfit)ifelseval(isempty(Tpar.(ej.name){i}),Tpar.default(i),str2num(Tpar.(ej.name){i})); % it's just a numeric value, fixed
%			pinit.(pname) = Tpar.default(i);
%			pxform.(pname) = Tpar.xform{i};
%			cvs.(pname) = Tpar.cv(i);
		elseif isempty(Tpar.(ej.name){i}) % then it's fit "globally' ie one value across all expts
			pmap.(pname) = @(pfit)pfit.(pname);
			pinit.(pname) = Tpar.default(i);
			pxform.(pname) = Tpar.xform{i};
			cvs.(pname) = Tpar.cv(i);
		else
			if isempty(str2num(Tpar.(ej.name){i}))
				% then it's a string and needs to be parsed out somehow
				tname = [pname '_' Tpar.(ej.name){i}];
				pmap.(pname) = @(pfit)pfit.(tname);
				pinit.(tname) = Tpar.default(i);
				pxform.(tname) = Tpar.xform{i};
				cvs.(tname) = Tpar.cv(i);
			else
				% then it's a value and it's fixed somehow
				pmap.(pname) = @(pfit)str2num(Tpar.(ej.name){i}); % fixed to # in expt
%				pinit.(pname) = feval(pmap.(pname)); % hopefully this doesn't get used anywhere
%				pxform.(pname) = Tpar.xform{i};
%				cvs.(pname) = Tpar.cv(i);
			end
		end
	end
	expt(j).pmap = pmap;
end
