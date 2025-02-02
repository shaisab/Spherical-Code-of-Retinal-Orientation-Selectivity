clear all

load('pooledMapAlphaCorr_OS data_40_OSI0.30_1.00_SI_0.00_dorsal_ventral_all_cell_types.mat');

dorsal=ismember(pooledmap.cellID,pooledmap.clusterdata{1, 1}.cellID);
ventral=ismember(pooledmap.cellID,pooledmap.clusterdata{1, 2}.cellID);
dorsal_OSI = pooledmap.OSI(dorsal)  ;
ventral_OSI = pooledmap.OSI(ventral)  ;

[h,p,ci,stats] = ttest2(dorsal_OSI,ventral_OSI)

% dorsal_ventral=double(ismember(pooledmap.cellID,pooledmap.clusterdata{1, 1}.cellID));
% dorsal_ventral=categorical(ismember(pooledmap.cellID,pooledmap.clusterdata{1, 1}.cellID),[true,false],["Dorsal","Ventral"]);
% OSI=pooledmap.OSI;
% 
% boxchart(dorsal_ventral,OSI)
% xlabel("Location")
% ylabel("OSI")