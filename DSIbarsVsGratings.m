

load('allCells bars and gratings.mat');
j=1;
for i=1:size(allCellsBars.cellID,1)
    if find(strcmp(allCellsGrating.cellID,allCellsBars.cellID{i}),1)
        tGrating(j)=find(strcmp(allCellsGrating.cellID,allCellsBars.cellID{i}),1);
        tBars(j)=i;
        j=j+1;
    end
end

DSIgrating=allCellsGrating.DSI(tGrating);
DSIbars=allCellsBars.DSI(tBars);
TYPEgrating=allCellsGrating.ONorONOFForOFF(tGrating);

% extract DSI for ON and ONOFF
j=1;
k=1;
for i=1:size(DSIgrating,1)
if TYPEgrating(i)==1
    ONgrating(j)=DSIgrating(i);
    ONbars(j)=DSIbars(i);
    j=j+1;
elseif TYPEgrating(i)==2
    ONOFFgrating(k)=DSIgrating(i);
    ONOFFbars(k)=DSIbars(i);
    k=k+1;
end
end

%% plot line histograms
        figure;
        hONgrating=histogram(ONgrating,10);
        valueHistONgrating=[hONgrating.Values,0];
        
        figure;
        hONbars=histogram(ONbars,10);
        valueHistONbars=[hONbars.Values,0];
 
        figure;
        hONOFFgrating=histogram(ONOFFgrating,10);
        valueHistONOFFgrating=[hONOFFgrating.Values,0];

        figure;
        hONOFFbars=histogram(ONOFFbars,10);
        valueHistONOFFbars=[hONOFFbars.Values,0];
      
               
%% prepare data for R code and plot errorbar plot
% only ON
dataON=[ONgrating';ONbars'];
groupON=[repmat('g',size(ONgrating',1),1);repmat('b',size(ONbars',1),1)];
GroupON=cellstr(groupON);

% only ONOFF
dataONOFF=[ONOFFgrating';ONOFFbars'];
groupONOFF=[repmat('g',size(ONOFFgrating',1),1);repmat('b',size(ONOFFbars',1),1)];
GroupONOFF=cellstr(groupONOFF);

% ON and ONOFF pooled
data=[DSIgrating;DSIbars];
group=[repmat('g',size(DSIgrating',1),1);repmat('b',size(DSIbars',1),1)];
GroupDSI=cellstr(group);


%% plot
medianON=[median(ONgrating),median(ONbars)];
qON=[quantile(ONgrating,[0.25 0.75]);quantile(ONbars,[0.25 0.75])];
medianONOFF=[median(ONOFFgrating),median(ONOFFbars)];
qONOFF=[quantile(ONOFFgrating,[0.25 0.75]);quantile(ONOFFbars,[0.25 0.75])];

figure; 
hold on;
errorbar([1,2],medianON,medianON'-qON(:,1),qON(:,2)-medianON','bo','MarkerSize',10,'LineWidth',2);
errorbar([1,2],medianONOFF,medianONOFF'-qONOFF(:,1),qONOFF(:,2)-medianONOFF','ro','MarkerSize',10,'LineWidth',2);
xlim([0.5 2.5]);
ylabel('DSI','FontSize',16); xlabel('Stimulus','FontSize',16);
set(gcf, 'color', [1 1 1]);
box off
ax = gca;
ax.FontSize=14;
ax.XTick = 1:2;
ax.XTickLabel={'gratings','bars'};

