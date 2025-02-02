% clear all
% load('pooledMapAlphaCorr_OS data_40_OSI0.20_1.00_SI_0.90_created_Dec_04_2022_12_50.mat');


for i=1:size(pooledmap.clusterdata{1, 56}.poriSz,1)
    if pooledmap.clusterdata{1, 56}.poriSz(i,1)>180
        pooledmap.clusterdata{1, 56}.poriSz(i,1)=360-pooledmap.clusterdata{1, 56}.poriSz(i,1);
    end

end

for i=1:size(pooledmap.clusterdata{1, 56}.poriSz,1)
    if pooledmap.clusterdata{1, 56}.poriSz(i,1)>45 && pooledmap.clusterdata{1, 56}.poriSz(i,1)<90
        pooledmap.clusterdata{1, 56}.poriSz(i,1)=90-pooledmap.clusterdata{1, 56}.poriSz(i,1);
    elseif pooledmap.clusterdata{1, 56}.poriSz(i,1)>90 && pooledmap.clusterdata{1, 56}.poriSz(i,1)<135
        pooledmap.clusterdata{1, 56}.poriSz(i,1)=pooledmap.clusterdata{1, 56}.poriSz(i,1)-90;
    elseif pooledmap.clusterdata{1, 56}.poriSz(i,1)>135 && pooledmap.clusterdata{1, 56}.poriSz(i,1)<180
        pooledmap.clusterdata{1, 56}.poriSz(i,1)=180-pooledmap.clusterdata{1, 56}.poriSz(i,1);
        
    end

end

mean_distance=mean(pooledmap.clusterdata{1, 56}.poriSz);
sd_distance=std(pooledmap.clusterdata{1, 56}.poriSz);

% figure
% 
% OSI=Meandistancefromcardinalorientations(:,3);
% Meandistancefromcardinalorientation=Meandistancefromcardinalorientations(:,1);
% err=Meandistancefromcardinalorientations(:,2);
% 
% scatter(OSI,Meandistancefromcardinalorientation)
% 
% hold on
% 
% % errorbar(OSI,Meandistancefromcardinalorientation,err,'LineStyle','none', 'Color', 'r','linewidth', 2)
% xlim([0 1])
% 
% xlabel('OSI')
% ylabel('Mean Distance from cardinal orientations')
% 
% hold off
% 
% figure
% 
% SI=Meandistancefromcardinalorientations(:,4);
% Meandistancefromcardinalorientation=Meandistancefromcardinalorientations(:,1);
% err=Meandistancefromcardinalorientations(:,2);
% 
% scatter(SI,Meandistancefromcardinalorientation)
% 
% hold on
% 
% % errorbar(SI,Meandistancefromcardinalorientation,err,'LineStyle','none', 'Color', 'r','linewidth', 2)
% xlim([0 1])
% 
% xlabel('SI')
% ylabel('Mean Distance from cardinal orientations')