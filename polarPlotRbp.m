
% plot a polar plot for all DS Rbp cells (from 'AI and angle for DS Rbp cells.mat' file) with the length of vectors proportional to AI

load('AI and angle for DS Rbp cells.mat');
figure;
% OFF arbor
%[X,Y]=pol2cart([deg2rad(angleOff);deg2rad(angleOff(1,1))],[AIoff;0.07]);
[X,Y]=pol2cart(deg2rad(angleOff),AIoff);
compassSS2p(X,Y,1,'g',2);
hold on;
% ON arbor
%[X,Y]=pol2cart([deg2rad(angleOn);deg2rad(angleOn(1,1))],[AIon;0.07]);
[X,Y]=pol2cart(deg2rad(angleOn),AIon);

compassSS2p(X,Y,1,'b',2);