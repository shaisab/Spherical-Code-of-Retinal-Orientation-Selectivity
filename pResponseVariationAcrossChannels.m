
% need to open the appropriate pooledmap file (with Itrans)

%% define parameters
slopePoints=5;

%% import data
pResN=pooledmap.pResponse(pooledmap.UVStvecComp.Itrans==1);
pResD=pooledmap.pResponse(pooledmap.UVStvecComp.Itrans==2);
pResT=pooledmap.pResponse(pooledmap.UVStvecComp.Itrans==3);
pResV=pooledmap.pResponse(pooledmap.UVStvecComp.Itrans==4);

onset1N=pooledmap.onset1(pooledmap.UVStvecComp.Itrans==1);
onset1D=pooledmap.onset1(pooledmap.UVStvecComp.Itrans==2);
onset1T=pooledmap.onset1(pooledmap.UVStvecComp.Itrans==3);
onset1V=pooledmap.onset1(pooledmap.UVStvecComp.Itrans==4);


%% prepare data 
channel={pResN,pResD,pResT,pResV};
for cn=1:size(channel,2)
    pResponse=channel{cn};
% determine the minimum number of elements
elements=[];
for nu=1:numel(pResponse)
    elements=[elements,numel(pResponse{nu})];
end

% standardize the length of responses % I manually set the number of elements to 243
for nu=1:numel(pResponse)
    numPresponseT{nu}=pResponse{nu}(1:243);
end
DSpResponse=reshape(cell2mat(numPresponseT'),243,[]);
clear numPresponseT

% normalize the responses
for p=1:size(DSpResponse,2)
    DSnorpResponse(:,p)=mat2gray(DSpResponse(:,p));
end

% smooth the reponses
for p=1:size(DSpResponse,2)
    DSpResponseS(:,p)=smooth(DSpResponse(:,p),15,'moving');
end

% normalize the smoothed responses
DSpResponseNorS=(DSpResponseS-repmat(mean(DSpResponseS(10:40,:),1),size(DSpResponseS,1),1))./(repmat(max(DSpResponseS,[],1),size(DSpResponseS,1),1)-repmat(mean(DSpResponseS(10:40,:),1),size(DSpResponseS,1),1));

DSresponse=DSpResponseNorS;

clear DSnorpResponse
clear DSpResponseNorS
clear DSpResponseS
%% extract this merged data file from Mar. 16 just to get the time vector
load('E:\Calcium imaging data\clusteringTest\mergedmeta_fov#11');
time=mergedmeta.meanRepDB.ROIdata{1,1}.summary.ordfluotime;
time=time(1:243,1);   % adjust the length of time to be the same as of norpResponseRS

%% calculate parameters for clustering all DS cells
figure;
DSpeakMinMag=DSresponse(130,:);     %data point number 130 is where the minima occur for both ON and ONOFF cells
for p=1:size(DSresponse,2)
    peakfinder(DSresponse(:,p),0.1,0.1,1,false); drawnow expose
    [peakLocMax,peakMagMax]=peakfinder(DSresponse(:,p),0.1,0.1,1,false);
    
    DSpeakLocMax{p}=peakLocMax;
    DSpeakMagMax{p}=peakMagMax;
    
%     %%%% excludes from analysis every cell with a weak and noisy response (more than 4 detected peaks)
%     if size(DSpeakLocMax{p},1)>4
%         DSresponse(:,p)=0;
%     end
    
    %%% find OFF peak
    offind=intersect(find(peakLocMax>145),find(peakLocMax<170));
    if isempty(offind)
        DSpeakLocOff(p)=NaN;
    else
        DSpeakLocOff(p)=peakLocMax(offind(1),1);
    end
    
    %%% calculate slope of OFF peak
    if ~isnan(DSpeakLocOff(p))
        yOffData=DSresponse(DSpeakLocOff(p):DSpeakLocOff(p)+15,p);
        % Set up fittype and options.
        ft = fittype('poly1');
        opts = fitoptions(ft);
        % Fit model to data.
        [fitresultOff, gofOff] = fit((1:size(yOffData,1))', yOffData(:,1), ft, opts );
        DSslopeOff(p)=fitresultOff.p1;
    else
        DSslopeOff(p)=NaN;
    end
       
    %%% find ON peak
    [peakLocMaxON,peakMagMaxON]=peakfinder(DSresponse(1:150,p),0.1,0.1,1,false);
    onind=find(peakLocMaxON>60);
    if isempty(onind)
        [peakLocMaxON,peakMagMaxON]=peakfinder(DSresponse(1:155,p),0.01,0.01,1,false);
        onind=find(peakLocMaxON>60);
        if isempty(onind) 
        DSpeakLocOn(p,1)=130;
        else
            DSpeakLocOn(p,1)=peakLocMaxON(onind(end),1);
        end
    else
        DSpeakLocOn(p,1)=peakLocMaxON(onind(1),1);
    end
    
    %%% calculate slope of ON peak
    if ~isnan(DSpeakLocOn(p,1))
        yOnData=DSresponse(DSpeakLocOn(p,1):DSpeakLocOn(p,1)+slopePoints,p);
        % Fit model to data.
        ft = fittype('poly1');
        opts = fitoptions(ft);
        [fitresultOn, gofOn] = fit((1:size(yOnData,1))', yOnData(:,1), ft, opts );
        DSslopeOn(p,1)=fitresultOn.p1;
    else
        DSslopeOn(p,1)=NaN;
    end
end

DSonSumResponse=sum(DSresponse(60:130,:),1);
DSoffSumResponse=sum(DSresponse(145:210,:),1);

DSpeakLocOn=time(DSpeakLocOn);      % convert DSpeakLocOn from number of data points to time

switch cn
    case 1
        onset1=onset1N;
    case 2
        onset1=onset1D;
    case 3
        onset1=onset1T;
    case 4
        onset1=onset1V;
end

save('tmp.mat');
DSpeakLocOn=DSpeakLocOn-onset1;    % adjust DSpeakLocOn relative to the onset of the stimulus

data=[DSpeakLocOn,DSslopeOn];

switch cn
    case 1
        dataN=data;
        DSresponseN=DSresponse;
    case 2
        dataD=data;
        DSresponseD=DSresponse;
    case 3
        dataT=data;
        DSresponseT=DSresponse;
    case 4
        dataV=data;
        DSresponseV=DSresponse;
end

clear data
clear DSpeakLocOn
clear DSslopeOn

end

%% adjust the DS response for the onset of stimulus for each cell

tmpN=repmat(time,1,size(DSresponseN,2))-repmat(onset1N',size(DSresponseN,1),1);
tmpD=repmat(time,1,size(DSresponseD,2))-repmat(onset1D',size(DSresponseD,1),1);
tmpT=repmat(time,1,size(DSresponseT,2))-repmat(onset1T',size(DSresponseT,1),1);
tmpV=repmat(time,1,size(DSresponseV,2))-repmat(onset1V',size(DSresponseV,1),1);


meanDSresponseN=mean(DSresponseN,2);
meanDSresponseD=mean(DSresponseD,2);
meanDSresponseT=mean(DSresponseT,2);
meanDSresponseV=mean(DSresponseV,2);
meanDSresponseNS=(meanDSresponseN-repmat(mean(meanDSresponseN(10:40,:),1),size(meanDSresponseN,1),1))./(repmat(max(meanDSresponseN,[],1),size(meanDSresponseN,1),1)-repmat(mean(meanDSresponseN(10:40,:),1),size(meanDSresponseN,1),1));
meanDSresponseDS=(meanDSresponseD-repmat(mean(meanDSresponseD(10:40,:),1),size(meanDSresponseD,1),1))./(repmat(max(meanDSresponseD,[],1),size(meanDSresponseD,1),1)-repmat(mean(meanDSresponseD(10:40,:),1),size(meanDSresponseD,1),1));
meanDSresponseTS=(meanDSresponseT-repmat(mean(meanDSresponseT(10:40,:),1),size(meanDSresponseT,1),1))./(repmat(max(meanDSresponseT,[],1),size(meanDSresponseT,1),1)-repmat(mean(meanDSresponseT(10:40,:),1),size(meanDSresponseT,1),1));
meanDSresponseVS=(meanDSresponseV-repmat(mean(meanDSresponseV(10:40,:),1),size(meanDSresponseV,1),1))./(repmat(max(meanDSresponseV,[],1),size(meanDSresponseV,1),1)-repmat(mean(meanDSresponseV(10:40,:),1),size(meanDSresponseV,1),1));

figure; 
hold on;
plot(mean(tmpN,2),meanDSresponseNS,'Color',[0/255 93/255 207/255],'LineStyle','-','LineWidth',2);
plot(mean(tmpT,2),meanDSresponseTS,'Color',[0/255 129/255 0/255],'LineStyle','-','LineWidth',2);
plot(mean(tmpD,2),meanDSresponseDS,'Color',[255/255 99/255 27/255],'LineStyle','-','LineWidth',2);
plot(mean(tmpV,2),meanDSresponseVS,'Color',[255/255 7/255 122/255],'LineStyle','-','LineWidth',2);
set(gcf, 'color', [1 1 1]);
box off
ax = gca;
ax.FontSize=16;
ylim([-0.1 1]);
legend(ax,'N','T','D','V');
legend('boxoff');
xlabel('Time (sec)','FontSize',20);ylabel('\DeltaF/F','FontSize',18);



% iTime=mean([tmpN,tmpD,tmpT,tmpV],2);
% for t=1:size(DSresponseN,2)
% iDSresponseN(:,t)=interp1(tmpN(:,t),DSresponseN(:,t),iTime);
% end
% for t=1:size(DSresponseD,2)
% iDSresponseD(:,t)=interp1(tmpD(:,t),DSresponseD(:,t),iTime);
% end
% for t=1:size(DSresponseT,2)
% iDSresponseT(:,t)=interp1(tmpT(:,t),DSresponseT(:,t),iTime);
% end
% for t=1:size(DSresponseV,2)
% iDSresponseV(:,t)=interp1(tmpV(:,t),DSresponseV(:,t),iTime);
% end
% 
% meaniDSresponse=mean(iDSresponse,2);
% figure;
% plot(iTime,iDSresponse);
% figure;
% hold on;
% plot(iTime,meaniDSresponse,'r');
% figure;
% plot(tmpN,DSresponseN)

%% plot pResponse for each channel
% f11=figure;
% set(f11,'position',[5 150 950 400]);
% subplot(1,4,1);
% shadedErrorBar(time,mean(DSresponseN,2),std(iDSresponseN,[],2),'k');
% xlim([0 16]); ylim([0 1]);
% subplot(1,4,2);
% shadedErrorBar(time,mean(DSresponseD,2),std(iDSresponseD,[],2),'k');
% xlim([0 16]); ylim([0 1]);
% subplot(1,4,3);
% shadedErrorBar(time,mean(DSresponseT,2),std(iDSresponseT,[],2),'k');
% xlim([0 16]); ylim([0 1]);
% subplot(1,4,4);
% shadedErrorBar(time,mean(DSresponseV,2),std(iDSresponseV,[],2),'k');
% xlim([0 16]); ylim([0 1]);


%% plot line histograms for time to ON peak
        figure;
        hlocN=histogram(dataN(:,1),20);
        valueHistlocN=[hlocN.Values,0];
%         hold on
%         plot(hlocN.BinEdges,valueHistlocN);

        figure;
        hlocD=histogram(dataD(:,1),20);
        valueHistlocD=[hlocD.Values,0];
%         hold on
%         plot(hlocN.BinEdges,valueHistlocN);
        
        figure;
        hlocT=histogram(dataT(:,1),20);
        valueHistlocT=[hlocT.Values,0];
%         hold on
%         plot(hlocN.BinEdges,valueHistlocN); 

        figure;
        hlocV=histogram(dataV(:,1),20);
        valueHistlocV=[hlocV.Values,0];
%         hold on
%         plot(hlocN.BinEdges,valueHistlocN);
        
figure;
hold on
plot(hlocN.BinEdges,valueHistlocN,'b','LineWidth',2);
plot(hlocD.BinEdges,valueHistlocD,'g','LineWidth',2);
plot(hlocT.BinEdges,valueHistlocT,'r','LineWidth',2);
plot(hlocV.BinEdges,valueHistlocV,'m','LineWidth',2);
xlabel('Time to ON peak (sec)','FontSize',16); ylabel('Number of cells','FontSize',16);
set(gcf, 'color', [1 1 1]);
box off
ax = gca;
ax.FontSize=14;
xlim([1 2.5]);

figure;
hold on
plot(hlocN.BinEdges,mat2gray(valueHistlocN),'b','LineWidth',2);
plot(hlocD.BinEdges,mat2gray(valueHistlocD),'g','LineWidth',2);
plot(hlocT.BinEdges,mat2gray(valueHistlocT),'r','LineWidth',2);
plot(hlocV.BinEdges,mat2gray(valueHistlocV),'m','LineWidth',2);
xlabel('Time to ON peak (sec)','FontSize',16); ylabel('Number of cells','FontSize',16);
set(gcf, 'color', [1 1 1]);
box off
ax = gca;
ax.FontSize=14;
xlim([1 2.5]);
        
%% plot line histograms for slope following ON peak
        figure;
        hslopeN=histogram(dataN(:,2),20);
        valueHistslopeN=[hslopeN.Values,0];
%         hold on
%         plot(hslopeN.BinEdges,valueHistslopeN);

        figure;
        hslopeD=histogram(dataD(:,2),20);
        valueHistslopeD=[hslopeD.Values,0];
%         hold on
%         plot(hslopeN.BinEdges,valueHistslopeN);
        
        figure;
        hslopeT=histogram(dataT(:,2),20);
        valueHistslopeT=[hslopeT.Values,0];
%         hold on
%         plot(hslopeN.BinEdges,valueHistslopeN); 

        figure;
        hslopeV=histogram(dataV(:,2),20);
        valueHistslopeV=[hslopeV.Values,0];
%         hold on
%         plot(hslopeN.BinEdges,valueHistslopeN);
        
figure;
hold on
plot(hslopeN.BinEdges,valueHistslopeN,'b','LineWidth',2);
plot(hslopeD.BinEdges,valueHistslopeD,'g','LineWidth',2);
plot(hslopeT.BinEdges,valueHistslopeT,'r','LineWidth',2);
plot(hslopeV.BinEdges,valueHistslopeV,'m','LineWidth',2);
xlabel('Slope following ON peak','FontSize',16); ylabel('Proportion of cells','FontSize',16);
set(gcf, 'color', [1 1 1]);
box off
ax = gca;
ax.FontSize=14;

figure;
hold on
plot(hslopeN.BinEdges,mat2gray(valueHistslopeN),'b','LineWidth',2);
plot(hslopeD.BinEdges,mat2gray(valueHistslopeD),'g','LineWidth',2);
plot(hslopeT.BinEdges,mat2gray(valueHistslopeT),'r','LineWidth',2);
plot(hslopeV.BinEdges,mat2gray(valueHistslopeV),'m','LineWidth',2);
xlabel('Slope following ON peak','FontSize',16); ylabel('Proportion of cells','FontSize',16);
set(gcf, 'color', [1 1 1]);
box off
ax = gca;
ax.FontSize=14;

%% prepare data for R code and plot errorbar plot
dataLoc=[dataN(:,1);dataD(:,1);dataT(:,1);dataV(:,1)];
dataSlope=[dataN(:,2);dataD(:,2);dataT(:,2);dataV(:,2)];
group=[repmat('v1',size(dataN,1),1);repmat('v2',size(dataD,1),1);repmat('v3',size(dataT,1),1);repmat('v4',size(dataV,1),1)];
Group=cellstr(group);


%% calculate statistics
SDloc=[std(dataN(:,1)),std(dataD(:,1)),std(dataT(:,1)),std(dataV(:,1))];
SEMloc=[std(dataN(:,1))/sqrt(length(dataN(:,1))),std(dataD(:,1))/sqrt(length(dataD(:,1))),std(dataT(:,1))/sqrt(length(dataT(:,1))),std(dataV(:,1))/sqrt(length(dataV(:,1)))];
meanLoc=[mean(dataN(:,1)),mean(dataD(:,1)),mean(dataT(:,1)),mean(dataV(:,1))];
medianLoc=[median(dataN(:,1)),median(dataD(:,1)),median(dataT(:,1)),median(dataV(:,1))];
qLoc=[quantile(dataN(:,1),[0.25 0.75]);quantile(dataD(:,1),[0.25 0.75]);quantile(dataT(:,1),[0.25 0.75]);quantile(dataV(:,1),[0.25 0.75])];

SDslope=[std(dataN(:,2)),std(dataD(:,2)),std(dataT(:,2)),std(dataV(:,2))];
SEMslope=[std(dataN(:,2))/sqrt(length(dataN(:,2))),std(dataD(:,2))/sqrt(length(dataD(:,2))),std(dataT(:,2))/sqrt(length(dataT(:,2))),std(dataV(:,2))/sqrt(length(dataV(:,2)))];
meanSlope=[mean(dataN(:,2)),mean(dataD(:,2)),mean(dataT(:,2)),mean(dataV(:,2))];
medianSlope=[median(dataN(:,2)),median(dataD(:,2)),median(dataT(:,2)),median(dataV(:,2))];
qSlope=[quantile(dataN(:,2),[0.25 0.75]);quantile(dataD(:,2),[0.25 0.75]);quantile(dataT(:,2),[0.25 0.75]);quantile(dataV(:,2),[0.25 0.75])];


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

% figure; 
% hold on;
% errorbar(1,medianLoc(1),SEMloc(1,1),SEMloc(1,1),'bo','MarkerSize',10,'LineWidth',2);
% errorbar(2,medianLoc(2),SEMloc(1,2),SEMloc(1,2),'go','MarkerSize',10,'LineWidth',2)
% errorbar(3,medianLoc(3),SEMloc(1,3),SEMloc(1,3),'ro','MarkerSize',10,'LineWidth',2)
% errorbar(4,medianLoc(4),SEMloc(1,4),SEMloc(1,4),'mo','MarkerSize',10,'LineWidth',2)
% ylabel('Time to ON peak (sec)','FontSize',16); xlabel('Channel','FontSize',16);
% set(gcf, 'color', [1 1 1]);
% box off
% ax = gca;
% ax.FontSize=14;
% ax.XTick = 1:4;
% ax.XTickLabel={'nasal','dorsal','temporal','ventral'};


figure; 
hold on;
errorbar(1,medianLoc(1),medianLoc(1)'-qLoc(1,1),qLoc(1,2)-medianLoc(1)','bo','MarkerSize',10,'LineWidth',2);
errorbar(2,medianLoc(2),medianLoc(2)'-qLoc(2,1),qLoc(2,2)-medianLoc(2)','go','MarkerSize',10,'LineWidth',2)
errorbar(3,medianLoc(3),medianLoc(3)'-qLoc(3,1),qLoc(3,2)-medianLoc(3)','ro','MarkerSize',10,'LineWidth',2)
errorbar(4,medianLoc(4),medianLoc(4)'-qLoc(4,1),qLoc(4,2)-medianLoc(4)','mo','MarkerSize',10,'LineWidth',2)
ylabel('Time to ON peak (sec)','FontSize',16); xlabel('Channel','FontSize',16);
set(gcf, 'color', [1 1 1]);
box off
ax = gca;
ax.FontSize=14;
ax.XTick = 1:4;
ax.XTickLabel={'nasal','dorsal','temporal','ventral'};

figure; 
hold on;
errorbar(1,medianSlope(1),medianSlope(1)'-qSlope(1,1),qSlope(1,2)-medianSlope(1)','bo','MarkerSize',10,'LineWidth',2);
errorbar(2,medianSlope(2),medianSlope(2)'-qSlope(2,1),qSlope(2,2)-medianSlope(2)','go','MarkerSize',10,'LineWidth',2)
errorbar(3,medianSlope(3),medianSlope(3)'-qSlope(3,1),qSlope(3,2)-medianSlope(3)','ro','MarkerSize',10,'LineWidth',2)
errorbar(4,medianSlope(4),medianSlope(4)'-qSlope(4,1),qSlope(4,2)-medianSlope(4)','mo','MarkerSize',10,'LineWidth',2)
ylabel('Slope following ON peak','FontSize',16); xlabel('Channel','FontSize',16);
set(gcf, 'color', [1 1 1]);
box off
ax = gca;
ax.FontSize=14;
ax.XTick = 1:4;
ax.XTickLabel={'nasal','dorsal','temporal','ventral'};

%% plot raw data along with medians or means
% figure
% groupn=[repmat(1,size(DSIN,1),1);repmat(2,size(DSID,1),1);repmat(3,size(DSIT,1),1);repmat(4,size(DSIV,1),1)];
% notBoxPlot(data,groupn)

%%
% check normality for loc
normality=[kstest(dataN(:,1)),kstest(dataD(:,1)),kstest(dataT(:,1)),kstest(dataV(:,1))];
% check homogenity of variance for loc
vartestn(dataLoc,group,'TestType','LeveneAbsolute')

% check normality for slope
normality=[kstest(dataN(:,2)),kstest(dataD(:,2)),kstest(dataT(:,2)),kstest(dataV(:,2))];
% check homogenity of variance for slope
vartestn(dataSlope,group,'TestType','LeveneAbsolute')


% 
% % parametric ANOVA
[p,tbl,stats]=anova1(dataLoc,group)
figure;
[c,~,~,gnames]=multcompare(stats)

[p,tbl,stats]=anova1(dataSlope,group)
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

save('tmp.mat');

