clear all
close all

load('pooledMapAlphaCorr_OS data_40_OSI0.20_1.00_created_Nov_22_2022_13_37_On_Off_Trans_No_Cutoffs.mat');

d=[0,45,90,135,180,225,270,315,0];
d2=[0,45,90,135,180,225,270,315];

type=(pooledmap.symmetry_ratio<0.99).*(pooledmap.symmetry_ratio>0.9).*(pooledmap.OSI>0.48).*(pooledmap.OSI<0.5);

for i=1:length(pooledmap.symmetry_ratio)
        if type(i)==1
        figure 
        ap=pooledmap.meanrs(i,:);
        apm=[ap,ap(1)];
        P=polarSS2p(deg2rad(d), apm/max(apm),1,'b',2);
        hold on
        
        hold off
        end
  
end 
