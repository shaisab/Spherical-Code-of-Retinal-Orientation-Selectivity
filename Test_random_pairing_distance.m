clc
clear all

load('pooledMapAlphaCorr_OS data_40_OSI0.30_1.00_SI_0.00_created_Dec_11_2022_13_30_On_Off_Trans.mat');

p = randperm(length(pooledmap.OSI)); %randomly permute the index of OSGC type
p_trans=p';
p_1 = p (1:length(p)/2);
p_2 = p(length(p)/2+1:end);

allpairs=allcomb(p,p_trans);

% p_1 = randperm(length(pooledmap.OSI)); %randomly permute the first half
% p_2 = randi([length(pooledmap.OSI)/2+1 length(pooledmap.OSI)],1,length(pooledmap.OSI)/2); %randomly permute the second half

OSI_diff_RP_1_2 = abs(pooledmap.OSI(allpairs(:,1),:)-pooledmap.OSI(allpairs(:,2),:));

distance_X_1_2 = ((pooledmap.XYZvecComp.X(allpairs(:,1),:)-pooledmap.XYZvecComp.X(allpairs(:,2),:)).^2);

distance_Y_1_2 = ((pooledmap.XYZvecComp.Y(allpairs(:,1),:)-pooledmap.XYZvecComp.Y(allpairs(:,2),:)).^2);

distance_Z_1_2 = ((pooledmap.XYZvecComp.Z(allpairs(:,1),:)-pooledmap.XYZvecComp.Z(allpairs(:,2),:)).^2);

distance_XYZ_1_2 =(distance_X_1_2+distance_Y_1_2+distance_Z_1_2).^(1/2);


% OSI_diff_RP_1_2 = abs(pooledmap.OSI(p_1,:)-pooledmap.OSI(p_2,:));
% 
% distance_X_1_2 = ((pooledmap.XYZvecComp.X(p_1,:)-pooledmap.XYZvecComp.X(p_2,:)).^2);
% 
% distance_Y_1_2 = ((pooledmap.XYZvecComp.Y(p_1,:)-pooledmap.XYZvecComp.Y(p_2,:)).^2);
% 
% distance_Z_1_2 = ((pooledmap.XYZvecComp.Z(p_1,:)-pooledmap.XYZvecComp.Z(p_2,:)).^2);
% 
% distance_XYZ_1_2 =(distance_X_1_2+distance_Y_1_2+distance_Z_1_2).^(1/2);

R=0.505530000000000;
Arc= 2*R.*asin(distance_XYZ_1_2/(2*R));

h=histogram(Arc, OSI_diff_RP_1_2);
figure
scatter(Arc, OSI_diff_RP_1_2);
scatterhistogram(Arc, OSI_diff_RP_1_2);
