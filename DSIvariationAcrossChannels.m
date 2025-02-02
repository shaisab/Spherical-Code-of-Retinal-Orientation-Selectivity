
% need to open the appropriate pooledmap file (with Itrans)
%% import data
DSIN=pooledmap.DSI(pooledmap.UVStvecComp.Itrans==1);
DSID=pooledmap.DSI(pooledmap.UVStvecComp.Itrans==2);
DSIT=pooledmap.DSI(pooledmap.UVStvecComp.Itrans==3);
DSIV=pooledmap.DSI(pooledmap.UVStvecComp.Itrans==4);

%% plot line histograms
        figure;
        hDSIN=histogram(DSIN,10);
        valueHistDSIN=[hDSIN.Values,0];
%         hold on
%         plot(hDSIN.BinEdges,valueHistDSIN);
        
        figure;
        hDSID=histogram(DSID,10);
        valueHistDSID=[hDSID.Values,0];
%         hold on
%         plot(hDSID.BinEdges,valueHistDSID);
        
        figure;
        hDSIT=histogram(DSIT,10);
        valueHistDSIT=[hDSIT.Values,0];
%         hold on
%         plot(hDSIT.BinEdges,valueHistDSIT);
        
        figure;
        hDSIV=histogram(DSIV,10);
        valueHistDSIV=[hDSIV.Values,0];
%         hold on
%         plot(hDSIV.BinEdges,valueHistDSIV);
        
figure;
hold on
plot(hDSIN.BinEdges,valueHistDSIN,'b','LineWidth',2);
plot(hDSID.BinEdges,valueHistDSID,'g','LineWidth',2);
plot(hDSIT.BinEdges,valueHistDSIT,'r','LineWidth',2);
plot(hDSIV.BinEdges,valueHistDSIV,'m','LineWidth',2);
xlabel('DSI','FontSize',16); ylabel('Number of cells','FontSize',16);
set(gcf, 'color', [1 1 1]);
box off
ax = gca;
ax.FontSize=14;
        
%% prepare data for R code and plot errorbar plot
data=[DSIN;DSID;DSIT;DSIV];
group=[repmat('v1',size(DSIN,1),1);repmat('v2',size(DSID,1),1);repmat('v3',size(DSIT,1),1);repmat('v4',size(DSIV,1),1)];
Group=cellstr(group);

SD=[std(DSIN),std(DSID),std(DSIT),std(DSIV)];
SEM=[std(DSIN)/sqrt(length(DSIN)),std(DSID)/sqrt(length(DSID)),std(DSIT)/sqrt(length(DSIT)),std(DSIV)/sqrt(length(DSIV))];
meanDSI=[mean(DSIN),mean(DSID),mean(DSIT),mean(DSIV)];
medianDSI=[median(DSIN),median(DSID),median(DSIT),median(DSIV)];
%figure;
%errorbar([1,2,3,4],meanDSI,SD,'ro','MarkerSize',10,'LineWidth',2);

% figure;
% hold on;
% errorbar(1,meanDSI(1),SD(1),'bo','MarkerSize',10,'LineWidth',2);
% errorbar(2,meanDSI(2),SD(2),'go','MarkerSize',10,'LineWidth',2);
% errorbar(3,meanDSI(3),SD(3),'ro','MarkerSize',10,'LineWidth',2);
% errorbar(4,meanDSI(4),SD(4),'mo','MarkerSize',10,'LineWidth',2);
% ylabel('DSI','FontSize',16); xlabel('Channel','FontSize',16);
% set(gcf, 'color', [1 1 1]);
% box off
% ax = gca;
% ax.FontSize=14;
% ax.XTick = 1:4;
% ax.XTickLabel={'nasal','dorsal','temporal','ventral'};

qDSIN=[quantile(DSIN,[0.25 0.75]);quantile(DSID,[0.25 0.75]);quantile(DSIT,[0.25 0.75]);quantile(DSIV,[0.25 0.75])];
figure; 
hold on;
errorbar(1,medianDSI(1),medianDSI(1)'-qDSIN(1,1),qDSIN(1,2)-medianDSI(1)','bo','MarkerSize',10,'LineWidth',2);
errorbar(2,medianDSI(3),medianDSI(3)'-qDSIN(3,1),qDSIN(3,2)-medianDSI(3)','ro','MarkerSize',10,'LineWidth',2)
errorbar(3,medianDSI(2),medianDSI(2)'-qDSIN(2,1),qDSIN(2,2)-medianDSI(2)','go','MarkerSize',10,'LineWidth',2)
errorbar(4,medianDSI(4),medianDSI(4)'-qDSIN(4,1),qDSIN(4,2)-medianDSI(4)','mo','MarkerSize',10,'LineWidth',2)
ylabel('DSI','FontSize',16); xlabel('Channel','FontSize',16);
set(gcf, 'color', [1 1 1]);
box off
ax = gca;
ax.FontSize=14;
ax.XTick = 1:4;
ax.XTickLabel={'nasal','temporal','dorsal','ventral'};

%% plot raw data along with medians or means
% figure
% groupn=[repmat(1,size(DSIN,1),1);repmat(2,size(DSID,1),1);repmat(3,size(DSIT,1),1);repmat(4,size(DSIV,1),1)];
% notBoxPlot(data,groupn)

%%
% % check normality
% normality=[kstest(DSIN),kstest(DSID),kstest(DSIT),kstest(DSIV)];
% % check homogenity of variance
% vartestn(data,group,'TestType','LeveneAbsolute')
% 
% parametric ANOVA
[p,tbl,stats]=anova1(data,group)
figure;
[c,~,~,gnames]=multcompare(stats)
% 
% %%% USE R code instead
% % ANOVA permutation
% [F df]=anova1_cell({data,group})
% signif=1-fcdf(F,df(1),df(2))
% 
% [pval,Factual,Fdist]=randanova1(data,group,1000)
% 
% %ts=[tinv([0.025  0.975],length(DSIN)-1);tinv([0.025  0.975],length(DSID)-1);tinv([0.025  0.975],length(DSIT)-1);tinv([0.025  0.975],length(DSIV)-1)];      % T-Score
% %CI=[mean(DSIN)+ts(1,:).*SEM(1,1);mean(DSID)+ts(2,:).*SEM(1,2);mean(DSIT)+ts(3,:).*SEM(1,3);mean(DSIV)+ts(4,:).*SEM(1,4)];                      % Confidence Intervals
% 
% %figure; errorbar([1,2,3,4],meanDSI,CI(:,1),CI(:,2),'ro')
% figure; errorbar([1,2,3,4],stats.means,SEM,'ro')
% 

save('DSI as a function of translation channels_ON.mat');