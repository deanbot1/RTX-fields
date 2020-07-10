function [expt,pinit,pxform,cvs] = parse_par_expt(Tpar,Texp,varargin)
% takes in Tpar table and Texp experiment description table 
% optional arguments are 'PARNAME',PARVAL name,value pairs that override the parameter values with the PARVAL and doesn't fit that parameter later on... 
% it can handle a PARNAME that's expt-specific, like 'CD16_F158'
% BE CAREFUL, if 'CD16_F158' is a thing and you pass it 'CD16' it probably
% won't complain but it won't be doing what you think it is...
% returns:
% expt structure array with pmap functions set
% pinit, initial parameter guess
% pxform, param transformations
% cvs, a structure of cv values (0 means don't fit), nonzero is used as sigma in
% bayes penalty 

clamplist = {};
clampvals = [];

if nargin > 2
	for j = 1:2:length(varargin)
		clamplist{end+1} = varargin{j};
		clampvals(end+1) = varargin{j+1};
	end
end

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
		ipar = find(strcmp(clamplist,pname));
		if Tpar.cv(i) < eps | ~isempty(ipar)
			% it's fixed, no fitting
			if isempty(ipar)
				pmap.(pname) = @(pfit)ifelseval(isempty(Tpar.(ej.name){i}),Tpar.default(i),str2num(Tpar.(ej.name){i})); % it's just a numeric value, fixed
			else
				pmap.(pname) = @(pfit)clampvals(ipar); % override the Tpar value and keep it off the pfit list...
			end
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
				ipar = find(strcmp(clamplist,tname));
				if isempty(ipar)
					pmap.(pname) = @(pfit)pfit.(tname);
					pinit.(tname) = Tpar.default(i);
					pxform.(tname) = Tpar.xform{i};
					cvs.(tname) = Tpar.cv(i);
				else
					pmap.(pname) = @(pfit)clampvals(ipar); % but if it's in the clamplist just set its 'value' to the passed through value in varargin
				end
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
