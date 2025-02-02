
% working directory should be: E:\Hb9\All Hb9 cells

% v is a vector of pdir based on voltage recordings
% m is a vector of pdir based on dendritic morphology (angle  from soma to
% centroid)
% rho is circular correlation coefficient
% pval is the p value (should be lower than 0.05 to show correlation)
clear all
load('pdir calcium vs voltage.mat');

%% correlation for ON DS cells
v=vON;
m=cON;
[rhoON pvalON]=circ_corrcc(deg2rad(v), deg2rad(m));
n=size(v,1);

figure;
[X,Y]=pol2cart(deg2rad(v), 1);
compassSS2p(X,Y,1,'b',1);
hold on;
[X,Y]=pol2cart(deg2rad(m), 1);
compassSS2p(X,Y,1,'r',1);
hold off

%% correlation for ONOFF DS cells
v=vONOFF;
m=cONOFF;
[rhoONOFF pvalONOFF]=circ_corrcc(deg2rad(v), deg2rad(m));
n=size(v,1);

figure;
[X,Y]=pol2cart(deg2rad(v), 1);
compassSS2p(X,Y,1,'b',1);
hold on;
[X,Y]=pol2cart(deg2rad(m), 1);
compassSS2p(X,Y,1,'r',1);
hold off

%save('pDir physiology vs morphology summary.mat');