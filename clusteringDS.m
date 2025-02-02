



filedir=cd;
namesMap=getFileList(filedir,'map_',0,'anywhere');
for i=1:size(namesMap,2)
    if strfind(namesMap{i},'.mat')
        load(num2str(namesMap{i}));
    end
end

numPresponse=map.Presponse{1};
DSpResponse=reshape(cell2mat(map.Presponse(map.isDS==1)),numel(numPresponse),[]);
%DSpResponse=reshape(cell2mat(map.Presponse(map.isOS==1)),numel(numPresponse),[]);
for p=1:size(DSpResponse,2)
    DSnorpResponse(:,p)=mat2gray(DSpResponse(:,p));
end


% f11=figure;
% set(f11,'position',[5 150 950 800]);
% subplot(1,2,1);F
% plot(DSnorpResponse,'k');
% xlim([0 size(DSnorpResponse,1)]);
% subplot(1,2,2);
% plot(mean(DSnorpResponse,2));
% xlim([0 size(DSnorpResponse,1)]);


for p=1:size(DSpResponse,2)
    DSpResponseS(:,p)=smooth(DSpResponse(:,p),15,'moving');
end

DSpResponseNorS=(DSpResponseS-repmat(mean(DSpResponseS(10:40,:),1),size(DSpResponseS,1),1))./(repmat(max(DSpResponseS,[],1),size(DSpResponseS,1),1)-repmat(mean(DSpResponseS(10:40,:),1),size(DSpResponseS,1),1));


peakMinMag=DSpResponseNorS(130,:);     %data point number 130 is where the minima occur for both ON and ONOFF cells
for p=1:size(DSpResponse,2)
    peakfinder(DSpResponseNorS(:,p),0.1,0.1,1,false);
    [peakLocMax,peakMagMax]=peakfinder(DSpResponseNorS(:,p),0.1,0.1,1,false);
    
    DSpeakLocMax{p}=peakLocMax;
    DSpeakMagMax{p}=peakMagMax;
    
    %%%% excludes from analysis every cell with a weak and noisy response (more than 4 detected peaks)
    if size(DSpeakLocMax{p},1)>4
        DSpResponseNorS(:,p)=0;
    end
    
    %%% find OFF peak
    offind=intersect(find(peakLocMax>145),find(peakLocMax<170));
    if isempty(offind)
        DSpeakLocOff(p)=0;
    else
        DSpeakLocOff(p)=peakLocMax(offind(1),1);
    end
    
    %%% calculate slope of OFF peak
    if DSpeakLocOff(p)
        yOffData=DSpResponseNorS(DSpeakLocOff(p):DSpeakLocOff(p)+15,:);
        % Set up fittype and options.
        ft = fittype('poly1');
        opts = fitoptions(ft);
        % Fit model to data.
        [fitresultOff, gofOff] = fit((1:size(yOffData,1))', yOffData(:,i), ft, opts );
        slopeOff(p)=fitresultOff.p1;
    end
    
    
    %%% find ON peak
    onind=find(peakLocMax>60);
    if isempty(onind)
        DSpeakLocOn(p)=0;
    else
        DSpeakLocOn(p)=peakLocMax(onind(1),1);
    end
    
    %%% calculate slope of ON peak
    if DSpeakLocOn(p)
        yOnData=DSpResponseNorS(DSpeakLocOn(p):DSpeakLocOn(p)+15,:);
        % Fit model to data.
        [fitresultOn, gofOn] = fit((1:size(yOnData,1))', yOnData(:,i), ft, opts );
        slopeOn(p)=fitresultOn.p1;
    end
end
onPeak=logical(DSpeakLocOn);
offPeak=logical(DSpeakLocOff);
onSumResponse=sum(DSpResponseNorS(60:130,:),1);
offSumResponse=sum(DSpResponseNorS(145:210,:),1);

DSpeakLocOn=DSpeakLocOn';
DSpeakLocOff=DSpeakLocOff';
onSumResponse=onSumResponse';
offSumResponse=offSumResponse';
peakMinMag=peakMinMag';
slopeOn=slopeOn';
slopeOff=slopeOff';

data=[DSpeakLocOn,DSpeakLocOff,onSumResponse,offSumResponse,peakMinMag,slopeOn,slopeOff];


evaCH=evalclusters(data,'kmeans','CalinskiHarabasz','KList',[1:6])
evaDB=evalclusters(data,'kmeans','DaviesBouldin','KList',[1:6])
evaG=evalclusters(data,'kmeans','gap','KList',[1:6])
evaS=evalclusters(data,'kmeans','silhouette','KList',[1:6])


evaCHg=evalclusters(data,'linkage','CalinskiHarabasz','KList',[1:6])
tree = linkage(data,'average');
figure()
dendrogram(tree)

idx2=kmeans(data,2,'Display','final','Replicates',5);
figure;
[silh2,h] = silhouette(data,idx2,'cityblock');
h = gca;
h.Children.EdgeColor = [.8 .8 1];
xlabel 'Silhouette Value';
ylabel 'Cluster';

idx3=kmeans(data,3,'Display','final','Replicates',5);
figure;
[silh3,h] = silhouette(data,idx3,'cityblock');
h = gca;
h.Children.EdgeColor = [.8 .8 1];
xlabel 'Silhouette Value';
ylabel 'Cluster';

idx4=kmeans(data,4,'Display','final','Replicates',5);
figure;
[silh4,h] = silhouette(data,idx4,'cityblock');
h = gca;
h.Children.EdgeColor = [.8 .8 1];
xlabel 'Silhouette Value';
ylabel 'Cluster';

idx5=kmeans(data,5,'Display','final','Replicates',5);
figure;
[silh5,h] = silhouette(data,idx5,'cityblock');
h = gca;
h.Children.EdgeColor = [.8 .8 1];
xlabel 'Silhouette Value';
ylabel 'Cluster';

cluster2=mean(silh2);
cluster3=mean(silh3);
cluster4=mean(silh4);
cluster5=mean(silh5);


idx=idx2;

f11=figure;
set(f11,'position',[5 150 950 800]);
subplot(2,5,1);
plot(DSnorpResponse(:,idx==1),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,5,2);
plot(DSnorpResponse(:,idx==2),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,5,3);
plot(DSnorpResponse(:,idx==3),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,5,4);
plot(DSnorpResponse(:,idx==4),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,5,5);
plot(DSnorpResponse(:,idx==5),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,5,6);
plot(mean(DSnorpResponse(:,idx==1),2),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,5,7);
plot(mean(DSnorpResponse(:,idx==2),2),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,5,8);
plot(mean(DSnorpResponse(:,idx==3),2),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,5,9);
plot(mean(DSnorpResponse(:,idx==4),2),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,5,10);
plot(mean(DSnorpResponse(:,idx==5),2),'k');
xlim([0 size(DSnorpResponse,1)]);


sum(idx==1)
sum(idx==2)










p
%%%%% CLUSTERING %%%%%%
% X=isDSnorpResponse';
% [idx,C] = kmeans(tmp,3);
%
% figure;
% plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',12)
% hold on
% plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12)
% plot(C(:,1),C(:,2),'kx',...
%      'MarkerSize',15,'LineWidth',3)
% legend('Cluster 1','Cluster 2','Centroids',...
%        'Location','NW')
% title 'Cluster Assignments and Centroids'
% hold off