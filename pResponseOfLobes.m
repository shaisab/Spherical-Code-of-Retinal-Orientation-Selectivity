
clear all;
load('pooledMapAlphaCorr_allMapsAlphaCorr_40_DSI0.30_1.00_created_Jan_04_2016_16_56_ONDS_DSI0.3.mat');

bin=45;

N=find(pooledmap.clusterdata{1,bin}.UVStvecComp.VecUSt>0.8);        %nasal
T=find(pooledmap.clusterdata{1,bin}.UVStvecComp.VecUSt<-0.8);       %temporal
D=find(pooledmap.clusterdata{1,bin}.UVStvecComp.VecVSt>0.8);        %dorsal
V=find(pooledmap.clusterdata{1,bin}.UVStvecComp.VecVSt<-0.8);       %ventral

NcellID=pooledmap.clusterdata{1,bin}.cellID(N);
for i=1:size(N,1)
Nk(i)=find(strcmp(pooledmap.cellID,NcellID{i}),1);
end

DcellID=pooledmap.clusterdata{1,bin}.cellID(D);
for i=1:size(D,1)
Dk(i)=find(strcmp(pooledmap.cellID,DcellID{i}),1);
end

TcellID=pooledmap.clusterdata{1,bin}.cellID(T);
for i=1:size(T,1)
Tk(i)=find(strcmp(pooledmap.cellID,TcellID{i}),1);
end

VcellID=pooledmap.clusterdata{1,bin}.cellID(V);
for i=1:size(V,1)
Vk(i)=find(strcmp(pooledmap.cellID,VcellID{i}),1);
end


Nres=pooledmap.pResponse(Nk);
Dres=pooledmap.pResponse(Dk);
Tres=pooledmap.pResponse(Tk);
Vres=pooledmap.pResponse(Vk);

Nres=reshape(cell2mat(Nres),[],size(Nres,1));
meanNres=mean(Nres,2);

Dres=reshape(cell2mat(Dres),[],size(Dres,1));
meanDres=mean(Dres,2);

Tres=reshape(cell2mat(Tres),[],size(Tres,1));
meanTres=mean(Tres,2);

Vres=reshape(cell2mat(Vres),[],size(Vres,1));
meanVres=mean(Vres,2);

figure;
plot([meanNres,meanDres,meanTres,meanVres]);









