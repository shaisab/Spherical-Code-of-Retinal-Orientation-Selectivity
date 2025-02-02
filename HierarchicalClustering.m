load('clustering summary_DownsampledPCA_Clustering_1.mat');

%Calculate means of dataT by idxT

means=grpstats(dataT,idxT);

%Remove the irrelevant clusters

means(10,:)=[];
means(8,:)=[];
means(4,:)=[];

%Hierarchical cluster tree
%Single linkage, also called nearest neighbor, uses the smallest distance between objects in the two clusters.
%Complete linkage, also called farthest neighbor, uses the largest distance between objects in the two clusters.
%Average linkage uses the average distance between all pairs of objects in any two clusters.

% Z=linkage(means,'single');
% Z=linkage(means,'complete');
Z=linkage(means,'average');
dendrogram(Z)
