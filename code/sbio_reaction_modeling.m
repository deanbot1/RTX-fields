%% script that uses sbiosteadystate to determine steady state rates adn tests how it operates

clear all; close all; 


%% loop thru various RTX concentrations...
tic
RTX = 10.^[-4:.5:3]';  % initial RTX concentrations to consider
SStrimerRTXconst = SSadcc(RTX,true);
SStrimerRTXvar = SSadcc(RTX,false);
toc

%% plot the results
figure;
imax = @(vec)find(vec==max(vec));
labl = @(vec,txt,col)text(RTX(imax(vec)),vec(imax(vec)),txt,'color',col,'HorizontalAlignment','center','VerticalAlignment','bottom');
semilogx(RTX,SStrimerRTXvar,'r-','LineWidth',3); hold on
labl(SStrimerRTXvar,'RTX variable','r');
semilogx(RTX,SStrimerRTXconst,'b-','LineWidth',3);
labl(SStrimerRTXconst,'RTX constant','b');


xlabel('initial RTX');
ylabel('Trimer@SS');
