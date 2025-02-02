load('clustering summary_DownsampledPCA_Clustering_1.mat');

load('allCellsAlphaCorr_OS data_40_PPC_0.5_DSI0.17_1.00_created_Jul_13_2020_14_30_Clustering_1.mat');

%Calculate means of OSI mean and std by idx (cluster)

[meansOSI,stdsOSI]=grpstats(allCells.OSI,idx,{'mean','std'});

%Remove the irrelevant clusters

meansOSI(10,:)=[]; stdsOSI(10,:)=[];
meansOSI(8,:)=[]; stdsOSI(8,:)=[];
meansOSI(4,:)=[]; stdsOSI(4,:)=[];

figure

cluster = 1:7;
meanOSI = [0.494619773305063 0.541770932713327 0.514971231893537 0.554929172356129 0.518994644659273 0.569520561152293 0.496502098496190]';
errhigh = [0.0817880524130258 0.0918433093158456 0.0830732656487595 0.0995684843554752 0.0952782710854866 0.0952262185673227 0.0835871252954975];
errlow  = [-0.0817880524130258 -0.0918433093158456 -0.0830732656487595 -0.0995684843554752 -0.0952782710854866 -0.0952262185673227 -0.0835871252954975];

bar(cluster,meanOSI)                

hold on

er = errorbar(cluster,meanOSI,errlow,errhigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off
