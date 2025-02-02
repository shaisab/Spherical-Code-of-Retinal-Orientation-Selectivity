
% working directory should be: E:\Hb9\All Hb9 cells

% v is a vector of pdir based on voltage recordings
% m is a vector of pdir based on dendritic morphology (angle  from soma to
% centroid)
% rho is circular correlation coefficient
% pval is the p value (should be lower than 0.05 to show correlation)
clear all
load('all Hb9 voltage and morphology.mat');

DSIvCut=0.3;
DSImCut=0.03;
confidenceCut=90;

j=1;
for i=1:size(pdirV,1)
    if DSIv(i)>DSIvCut && DSIm(i)>DSImCut && confidenceM(i)>confidenceCut
        v(j,1)=pdirV(i);
        m(j,1)=pdirM(i);
        j=j+1;
    end
end
[rho pval]=circ_corrcc(deg2rad(v), deg2rad(m));
n=size(v,1);

v=v-14;     % adjust to account for rotation of the retinas 
m=m-14;

figure;
[X,Y]=pol2cart(deg2rad(v), 1);
compassSS2p(X,Y,1,'r',1);
hold on;
[X,Y]=pol2cart(deg2rad(m), 1);
compassSS2p(X,Y,1,'b',1);
hold off

save('pDir physiology vs morphology summary.mat');