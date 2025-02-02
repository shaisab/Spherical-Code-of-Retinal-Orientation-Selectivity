clear all

load('pooledMapAlphaCorr_OS data_40_OSI0.30_1.00_SI_0.00_dorsal_ventral_all_cell_types.mat');

dorsal=ismember(pooledmap.cellID,pooledmap.clusterdata{1, 1}.cellID);
ventral=ismember(pooledmap.cellID,pooledmap.clusterdata{1, 2}.cellID);
dorsal_OSI = pooledmap.OSI(dorsal)  ;
ventral_OSI = pooledmap.OSI(ventral)  ;

dorsal_OSI_on_sus = dorsal_OSI(ismember((pooledmap.clusterdata{1, 1}.cellID),pooledmap.cellID(pooledmap.typeClust==3)))  ;
ventral_OSI_on_sus = ventral_OSI(ismember((pooledmap.clusterdata{1, 2}.cellID),pooledmap.cellID(pooledmap.typeClust==3)))  ;
dorsal_cellid_on_sus = pooledmap.cellID(ismember((pooledmap.clusterdata{1, 1}.cellID),pooledmap.cellID(pooledmap.typeClust==3)))  ;
ventral_cellid_on_sus = pooledmap.cellID(ismember((pooledmap.clusterdata{1, 2}.cellID),pooledmap.cellID(pooledmap.typeClust==3)))  ;

dorsal_OSI_on_trans = dorsal_OSI(ismember((pooledmap.clusterdata{1, 1}.cellID),pooledmap.cellID(pooledmap.typeClust==8)))  ;
ventral_OSI_on_trans = ventral_OSI(ismember((pooledmap.clusterdata{1, 2}.cellID),pooledmap.cellID(pooledmap.typeClust==8)))  ;

dorsal_OSI_on_off_sus = dorsal_OSI(ismember((pooledmap.clusterdata{1, 1}.cellID),pooledmap.cellID(pooledmap.typeClust==4)))  ;
ventral_OSI_on_off_sus = ventral_OSI(ismember((pooledmap.clusterdata{1, 2}.cellID),pooledmap.cellID(pooledmap.typeClust==4)))  ;

dorsal_OSI_on_off_trans = dorsal_OSI(ismember((pooledmap.clusterdata{1, 1}.cellID),pooledmap.cellID(pooledmap.typeClust==6)))  ;
ventral_OSI_on_off_trans = ventral_OSI(ismember((pooledmap.clusterdata{1, 2}.cellID),pooledmap.cellID(pooledmap.typeClust==6)))  ;

[h,p,ci,stats] = ttest2(dorsal_OSI_on_trans,ventral_OSI_on_trans)

% dorsal_ventral=double(ismember(pooledmap.cellID,pooledmap.clusterdata{1, 1}.cellID));
% dorsal_ventral=categorical(ismember(pooledmap.cellID,pooledmap.clusterdata{1, 1}.cellID),[true,false],["Dorsal","Ventral"]);
% OSI=pooledmap.OSI;
% 
% boxchart(dorsal_ventral,OSI)
% xlabel("Location")
% ylabel("OSI")