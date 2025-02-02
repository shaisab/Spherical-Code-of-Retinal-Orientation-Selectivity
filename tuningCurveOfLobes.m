
clear all;
load('pooledMapAlphaCorr_allMapsAlphaCorr_40_DSI0.30_1.00_created_Jan_04_2016_16_56_ONDS_DSI0.3.mat');

bin=45;

N=find(pooledmap.clusterdata{1,bin}.UVStvecComp.VecUSt>0.8);        %nasal
T=find(pooledmap.clusterdata{1,bin}.UVStvecComp.VecUSt<-0.8);       %temporal
D=find(pooledmap.clusterdata{1,bin}.UVStvecComp.VecVSt>0.8);        %dorsal
V=find(pooledmap.clusterdata{1,bin}.UVStvecComp.VecVSt<-0.8);       %ventral

NcellID=pooledmap.clusterdata{1,bin}.cellID(N);
DcellID=pooledmap.clusterdata{1,bin}.cellID(D);
TcellID=pooledmap.clusterdata{1,bin}.cellID(T);
VcellID=pooledmap.clusterdata{1,bin}.cellID(V);


Tmeanr=[];
Tpdir=[];

celln=27;
t1=mergedmeta.meanRepDG.ROIdata{1,celln}.summary.meanr
t2=mergedmeta.meanRepDG.ROIdata{1,celln}.summary.pdir
Tmeanr=[Tmeanr;t1];
Tpdir=[Tpdir;t2];

save('T45.mat','Tmeanr','Tpdir')


















% Nres=pooledmap.pResponse(Nk);
% Dres=pooledmap.pResponse(Dk);
% Tres=pooledmap.pResponse(Tk);
% Vres=pooledmap.pResponse(Vk);
% 
% Nres=reshape(cell2mat(Nres),[],size(Nres,1));
% meanNres=mean(Nres,2);
% 
% Dres=reshape(cell2mat(Dres),[],size(Dres,1));
% meanDres=mean(Dres,2);
% 
% Tres=reshape(cell2mat(Tres),[],size(Tres,1));
% meanTres=mean(Tres,2);
% 
% Vres=reshape(cell2mat(Vres),[],size(Vres,1));
% meanVres=mean(Vres,2);
% 
% figure;
% plot([meanNres,meanDres,meanTres,meanVres]);









