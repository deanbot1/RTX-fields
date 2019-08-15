function plot_expt(expt,pmat,T)
% automatically plot experimental data in a nice way
% INPUTS
% EXPT is a 1d structure array of experiments with the following required
% fields:
% .model is a function @(pvec,T) where pvec is Np x 1 vector of parameters
% and T is a vector of time points at which to evaluate the model function.
% The model returns an NT x NY matrix, where NT is the number of time
% points, and NY is the number of distinct prediction variables
% .obs is a NTxNY matrix of observations from the ith experiment. It must
% be this size. It can be padded with NaN's if there are no observations
% available for a particular time/Ytype pair. 
% .time is a vector of times at which observations are available 
% .Ynames is a 1xNY cell array of names for the Y variables. It's used to
% match up variables for data fitting.
% .name is a string specifying a human-readable label for the experiment. Used only for plotting and display. 
% 
% PMAT is a Np x Ne matrix, where Ne = length(expt). Each ith column is a
% parameter vector understandable to the corresponding expt(i).model function. 
%
%
% $URL$
% $Author$
% $Rev$
% $Date$

Ne = length(expt);

allYnames = {};
for i = 1:Ne
	allYnames = union(allYnames,expt(i).Ynames);
end
Ny = length(allYnames);

subplot2 = @(i,j)(subplot(Ny,Ne,(i-1)*Ne+j));

Ymax = zeros(Ny,1);
hi = 0;
for i = 1:Ne
	jj = zeros(1,length(expt(i).Ynames));
	for j = 1:length(expt(i).Ynames)
		Yname = expt(i).Ynames{j};
		hi = hi+1;
		jj(j) = find(strcmp(allYnames,Yname));
		subplot2(jj(j),i);
		if isfield(expt,'obslo');
			YYY = expt(i).obs(:,j);
			errorbar(expt(i).time,YYY,YYY-expt(i).obslo(:,j),expt(i).obshi(:,j)-YYY,'r-'); hold on;
		else
			plot(expt(i).time,expt(i).obs(:,j),'ro','MarkerFaceColor','r'); hold on;
		end
		title([num2str(i) ':' expt(i).name]);
		ylabel(Yname);
	end
	if nargin > 1
		if nargin > 2
			tsim = T;
		else
			tsim = expt(i).time;
		end
		ysim = expt(i).model(pmat(:,i),tsim);
		%Ymax = max([Ymax,max(ysim)'],[],2);
		for j = 1:length(jj)
			subplot2(jj(j),i);
			plot(tsim,ysim(:,j),'b.-');
			Ymax(jj(j))=max(Ymax(jj(j)),max(ysim(:,j)));
		end
	end
end

for i = 1:Ny
	for j = 1:Ne
		subplot2(i,j)
		%set(gca,'YLim',[0 Ymax(i)]);
	end
end

%set(h,'YLim',[0 Ymax]);