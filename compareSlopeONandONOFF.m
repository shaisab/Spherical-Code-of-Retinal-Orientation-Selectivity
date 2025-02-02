

%% sort file order by ascending ND
folderName=cd;
fileName=getFileList(folderName,'summary_ON_DS',0,'anywhere');
ON=load(char(fileName));
fileName=getFileList(folderName,'summary_ONOFF_DS',0,'anywhere');
ONOFF=load(char(fileName));

%% plot line histograms of slope of ON and ONOFF DS

figure;
hold on
plot(ON.hslopeN.BinEdges,mat2gray(ON.valueHistslopeN),'Color',[0/255 93/255 207/255],'LineStyle',':','LineWidth',2);
plot(ON.hslopeT.BinEdges,mat2gray(ON.valueHistslopeT),'Color',[0/255 129/255 0/255],'LineStyle',':','LineWidth',2);
plot(ON.hslopeD.BinEdges,mat2gray(ON.valueHistslopeD),'Color',[255/255 99/255 27/255],'LineStyle',':','LineWidth',2);
plot(ON.hslopeV.BinEdges,mat2gray(ON.valueHistslopeV),'Color',[255/255 7/255 122/255],'LineStyle',':','LineWidth',2);
xlabel('Slope following ON peak','FontSize',18); ylabel('Proportion of cells','FontSize',18);
set(gcf, 'color', [1 1 1]);
box off

plot(ONOFF.hslopeN.BinEdges,mat2gray(ONOFF.valueHistslopeN),'Color',[0/255 93/255 207/255],'LineStyle','-','LineWidth',2);
plot(ONOFF.hslopeT.BinEdges,mat2gray(ONOFF.valueHistslopeT),'Color',[0/255 129/255 0/255],'LineStyle','-','LineWidth',2);
plot(ONOFF.hslopeD.BinEdges,mat2gray(ONOFF.valueHistslopeD),'Color',[255/255 99/255 27/255],'LineStyle','-','LineWidth',2);
plot(ONOFF.hslopeV.BinEdges,mat2gray(ONOFF.valueHistslopeV),'Color',[255/255 7/255 122/255],'LineStyle','-','LineWidth',2);
xlabel('Slope following ON peak','FontSize',18); ylabel('Proportion of cells','FontSize',18);
set(gcf, 'color', [1 1 1]);
box off
ax = gca;
ax.FontSize=16;
xlim([-0.15 0.05]);

%% plot line histograms of time to ON peak of ON and ONOFF DS
figure;
hold on
plot(ON.hlocN.BinEdges,mat2gray(ON.valueHistlocN),'Color',[0/255 93/255 207/255],'LineStyle',':','LineWidth',2);
plot(ON.hlocT.BinEdges,mat2gray(ON.valueHistlocT),'Color',[0/255 129/255 0/255],'LineStyle',':','LineWidth',2);
plot(ON.hlocD.BinEdges,mat2gray(ON.valueHistlocD),'Color',[255/255 99/255 27/255],'LineStyle',':','LineWidth',2);
plot(ON.hlocV.BinEdges,mat2gray(ON.valueHistlocV),'Color',[255/255 7/255 122/255],'LineStyle',':','LineWidth',2);
xlabel('Time to ON peak (sec)','FontSize',18); ylabel('Number of cells','FontSize',18);
set(gcf, 'color', [1 1 1]);
box off

plot(ONOFF.hlocN.BinEdges,mat2gray(ONOFF.valueHistlocN),'Color',[0/255 93/255 207/255],'LineStyle','-','LineWidth',2);
plot(ONOFF.hlocT.BinEdges,mat2gray(ONOFF.valueHistlocT),'Color',[0/255 129/255 0/255],'LineStyle','-','LineWidth',2);
plot(ONOFF.hlocD.BinEdges,mat2gray(ONOFF.valueHistlocD),'Color',[255/255 99/255 27/255],'LineStyle','-','LineWidth',2);
plot(ONOFF.hlocV.BinEdges,mat2gray(ONOFF.valueHistlocV),'Color',[255/255 7/255 122/255],'LineStyle','-','LineWidth',2);
xlabel('Time to ON peak (sec)','FontSize',18); ylabel('Number of cells','FontSize',18);
set(gcf, 'color', [1 1 1]);
box off
ax = gca;
ax.FontSize=16;
xlim([1 6.5]);
    
%% plot medians and first and third quartiles of slope of ON and ONOFF DS
figure; 
hold on;
errorbar(1,ON.medianSlope(1),ON.medianSlope(1)'-ON.qSlope(1,1),ON.qSlope(1,2)-ON.medianSlope(1)','Color',[0/255 93/255 207/255],'Marker','o','MarkerSize',10,'LineWidth',2);
errorbar(2,ON.medianSlope(3),ON.medianSlope(3)'-ON.qSlope(3,1),ON.qSlope(3,2)-ON.medianSlope(3)','Color',[0/255 129/255 0/255],'Marker','o','MarkerSize',10,'LineWidth',2)
errorbar(3,ON.medianSlope(2),ON.medianSlope(2)'-ON.qSlope(2,1),ON.qSlope(2,2)-ON.medianSlope(2)','Color',[255/255 99/255 27/255],'Marker','o','MarkerSize',10,'LineWidth',2)
errorbar(4,ON.medianSlope(4),ON.medianSlope(4)'-ON.qSlope(4,1),ON.qSlope(4,2)-ON.medianSlope(4)','Color',[255/255 7/255 122/255],'Marker','o','MarkerSize',10,'LineWidth',2)
ylabel('Slope following ON peak','FontSize',18); xlabel('Channel','FontSize',18);
set(gcf, 'color', [1 1 1]);
box off

errorbar(1,ONOFF.medianSlope(1),ONOFF.medianSlope(1)'-ONOFF.qSlope(1,1),ONOFF.qSlope(1,2)-ONOFF.medianSlope(1)','Color',[0/255 93/255 207/255],'Marker','o','MarkerSize',10,'MarkerFaceColor',[0/255 93/255 207/255],'LineWidth',2);
errorbar(2,ONOFF.medianSlope(3),ONOFF.medianSlope(3)'-ONOFF.qSlope(3,1),ONOFF.qSlope(3,2)-ONOFF.medianSlope(3)','Color',[0/255 129/255 0/255],'Marker','o','MarkerSize',10,'MarkerFaceColor',[0/255 129/255 0/255],'LineWidth',2)
errorbar(3,ONOFF.medianSlope(2),ONOFF.medianSlope(2)'-ONOFF.qSlope(2,1),ONOFF.qSlope(2,2)-ONOFF.medianSlope(2)','Color',[255/255 99/255 27/255],'Marker','o','MarkerSize',10,'MarkerFaceColor',[255/255 99/255 27/255],'LineWidth',2)
errorbar(4,ONOFF.medianSlope(4),ONOFF.medianSlope(4)'-ONOFF.qSlope(4,1),ONOFF.qSlope(4,2)-ONOFF.medianSlope(4)','Color',[255/255 7/255 122/255],'Marker','o','MarkerSize',10,'MarkerFaceColor',[255/255 7/255 122/255],'LineWidth',2)
ylabel('Slope following ON peak','FontSize',18); xlabel('Channel','FontSize',18);
set(gcf, 'color', [1 1 1]);
box off
ax = gca;
ax.FontSize=16;
ax.XTick = 1:4;
ax.XTickLabel={'N','T','D','V'};

%% plot medians and first and third quartiles of timr to ON peak of ON and ONOFF DS
figure; 
hold on;
errorbar(1,ON.medianLoc(1),ON.medianLoc(1)'-ON.qLoc(1,1),ON.qLoc(1,2)-ON.medianLoc(1)','Color',[0/255 93/255 207/255],'Marker','o','MarkerSize',10,'LineWidth',2);
errorbar(2,ON.medianLoc(3),ON.medianLoc(3)'-ON.qLoc(3,1),ON.qLoc(3,2)-ON.medianLoc(3)','Color',[0/255 129/255 0/255],'Marker','o','MarkerSize',10,'LineWidth',2)
errorbar(3,ON.medianLoc(2),ON.medianLoc(2)'-ON.qLoc(2,1),ON.qLoc(2,2)-ON.medianLoc(2)','Color',[255/255 99/255 27/255],'Marker','o','MarkerSize',10,'LineWidth',2)
errorbar(4,ON.medianLoc(4),ON.medianLoc(4)'-ON.qLoc(4,1),ON.qLoc(4,2)-ON.medianLoc(4)','Color',[255/255 7/255 122/255],'Marker','o','MarkerSize',10,'LineWidth',2)
ylabel('Time to ON peak (sec)','FontSize',18); xlabel('Channel','FontSize',18);
set(gcf, 'color', [1 1 1]);
box off

errorbar(1,ONOFF.medianLoc(1),ONOFF.medianLoc(1)'-ONOFF.qLoc(1,1),ONOFF.qLoc(1,2)-ONOFF.medianLoc(1)','Color',[0/255 93/255 207/255],'Marker','o','MarkerSize',10,'MarkerFaceColor',[0/255 93/255 207/255],'LineWidth',2);
errorbar(2,ONOFF.medianLoc(3),ONOFF.medianLoc(3)'-ONOFF.qLoc(3,1),ONOFF.qLoc(3,2)-ONOFF.medianLoc(3)','Color',[0/255 129/255 0/255],'Marker','o','MarkerSize',10,'MarkerFaceColor',[0/255 129/255 0/255],'LineWidth',2)
errorbar(3,ONOFF.medianLoc(2),ONOFF.medianLoc(2)'-ONOFF.qLoc(2,1),ONOFF.qLoc(2,2)-ONOFF.medianLoc(2)','Color',[255/255 99/255 27/255],'Marker','o','MarkerSize',10,'MarkerFaceColor',[255/255 99/255 27/255],'LineWidth',2)
errorbar(4,ONOFF.medianLoc(4),ONOFF.medianLoc(4)'-ONOFF.qLoc(4,1),ONOFF.qLoc(4,2)-ONOFF.medianLoc(4)','Color',[255/255 7/255 122/255],'Marker','o','MarkerSize',10,'MarkerFaceColor',[255/255 7/255 122/255],'LineWidth',2)
ylabel('Time to ON peak (sec)','FontSize',18); xlabel('Channel','FontSize',18);
set(gcf, 'color', [1 1 1]);
box off
ax = gca;
ax.FontSize=16;
ax.XTick = 1:4;
ax.XTickLabel={'N','T','D','V'};

%% prepare data for R analysis
ONdata=[ON.dataN;ON.dataD;ON.dataT;ON.dataV];
ONOFFdata=[ONOFF.dataN;ONOFF.dataD;ONOFF.dataT;ONOFF.dataV];










