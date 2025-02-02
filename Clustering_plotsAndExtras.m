

%%%Figures for lab meeting

% plots the 15 first PCs of DS cells
figure;
for i=1:15
subplot(3,5,i)
plot(score(:,i))
xlim([0 size(score,1)]); ylim([min(min(score)) max(max(score))]);
title(['PC ',num2str(i)]);
end


% plots the 15 first PCs of SNS cells
figure;
for i=1:15
subplot(3,5,i)
plot(SNSscore(:,i))
xlim([0 size(SNSscore,1)]); ylim([min(min(SNSscore)) max(max(SNSscore))]);
title(['PC ',num2str(i)]);
end

% plots the 15 DS cells in the data set for example
figure;
for i=1:15
subplot(3,5,i)
plot(DSpResponseNorS(:,i))
xlim([0 size(DSpResponseNorS,1)]); ylim([-0.1 max(max(DSpResponseNorS))]);
end


% plots the amount of variance explained
figure;
plot(explained, '-or');


% plots individual and mean responses of all DS cells clustered into 2
% groups.
figure;
set(f11,'position',[5 150 950 800]);
subplot(2,2,1);
plot(DSnorpResponse(:,idx==1),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,2,2);
plot(DSnorpResponse(:,idx==2),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,2,3);
shadedErrorBar(1:size(DSnorpResponse,1),mean(DSnorpResponse(:,idx==1),2),std(DSnorpResponse(:,idx==1),[],2),'k');
xlim([0 size(DSnorpResponse,1)]); ylim([0 1]);
subplot(2,2,4);
shadedErrorBar(1:size(DSnorpResponse,1),mean(DSnorpResponse(:,idx==2),2),std(DSnorpResponse(:,idx==2),[],2),'k');
xlim([0 size(DSnorpResponse,1)]); ylim([0 1]);

% plorts the AIC and BIC
figure;
plot(allBIC,'-ok');
hold on
plot(allAIC,'-ob');

figure;
plot(bic,'-ok');
hold on
plot(aic,'-ob');

figure;
plot(deltaBIC,'-or');

% plots the clustered SNS cells responses
figure;
set(f11,'position',[5 150 950 800]);
subplot(2,2,1);
plot(SNSpResponsenorS(:,SNSidx==1),'k');
xlim([0 size(SNSpResponsenorS,1)]);
subplot(2,2,2);
plot(SNSpResponsenorS(:,SNSidx==2),'k');
xlim([0 size(SNSpResponsenorS,1)]);
subplot(2,2,3);
shadedErrorBar(1:size(SNSpResponsenorS,1),mean(SNSpResponsenorS(:,SNSidx==1),2),std(SNSpResponsenorS(:,SNSidx==1),[],2),'k');
xlim([0 size(SNSpResponsenorS,1)]); ylim([0 1]);
subplot(2,2,4);
shadedErrorBar(1:size(SNSpResponsenorS,1),mean(SNSpResponsenorS(:,SNSidx==2),2),std(SNSpResponsenorS(:,SNSidx==2),[],2),'k');
xlim([0 size(SNSpResponsenorS,1)]); ylim([0 1]);


%%
%%%%%%%%%  EXTRA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use Gaussian mixture models for clusterring
% repeat the evaluation of the number of clusters (every time it uses a different initial value) to determine the optimal number of clusters 

clustFcn=@(X,K) (kmeans(X, K,'Replicates',50,'Options',options));
evaKMCH=evalclusters(data,clustFcn,'CalinskiHarabasz','KList',1:10)
evaKMDB=evalclusters(data,clustFcn,'DaviesBouldin','KList',1:10)
evaKMS=evalclusters(data,clustFcn,'silhouette','KList',1:10)
evaKMG=evalclusters(data,clustFcn,'gap','KList',1:10)

figure; plot(evaKMCH); title(['Optimal K = ',num2str(evaKMCH.OptimalK)]);   %Calinski Harabasz
figure; plot(evaKMDB); title(['Optimal K = ',num2str(evaKMDB.OptimalK)]);   %DaviesBouldin
figure; plot(evaKMS); title(['Optimal K = ',num2str(evaKMS.OptimalK)]);     %silhouette
figure; plot(evaKMG); title(['Optimal K = ',num2str(evaKMG.OptimalK)]);     %gap

idx1CH=sum(evaKMCH.OptimalY==1)
idx2CH=sum(evaKMCH.OptimalY==2)
idx3CH=sum(evaKMCH.OptimalY==3)
idx4CH=sum(evaKMCH.OptimalY==4)
idx5CH=sum(evaKMCH.OptimalY==5)

idxS1=sum(evaGMS.OptimalY==1)
idx2S=sum(evaGMS.OptimalY==2)
idx3S=sum(evaGMS.OptimalY==3)
idx4S=sum(evaGMS.OptimalY==4)
idx5S=sum(evaGMS.OptimalY==5)

idx1G=sum(evaGMS.OptimalY==1)
idx2G=sum(evaGMS.OptimalY==2)
idx3G=sum(evaGMS.OptimalY==3)
idx4G=sum(evaGMS.OptimalY==4)
idx5G=sum(evaGMS.OptimalY==5)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hierarchical Clustering
evaCHg=evalclusters(data,'linkage','CalinskiHarabasz','KList',[1:10])
tree = linkage(data,'average');
figure()
dendrogram(tree)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% kmean clustering

idx2=kmeans(data,2,'Display','final','Replicates',50);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% allOptimalK(nevaGM)=evaGM.OptimalK;
% allOptimalY{nevaGM}=evaGM.OptimalY;
% end
% ind=find(allOptimalK==mode(allOptimalK));
% finalOptimalK=allOptimalK(ind(1));
%finalOptimalY=allOptimalY{finalOptimalK(1)};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GM=fitgmdist(data,2,'CovarianceType','diagonal','Replicates',10,'Options',options);
% GM=fitgmdist(data,2,'Start','plus','Replicates',10,'Options',options);