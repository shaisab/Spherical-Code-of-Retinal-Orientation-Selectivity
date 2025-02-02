
function []= plotAverageResponsePerType(ONorONOFForOFF,pResponse)

% extract this merged data file from Mar. 16 just to get the time vector
load('clusteringTest\mergedmeta_fov#11');
time=mergedmeta.meanRepDB.ROIdata{1,1}.summary.ordfluotime;

% plot subplots of pResponse of DS cells by type

elements=[];
for nu=1:length(pResponse)   
    elements=[elements,numel(pResponse{nu})];
end
for nu=1:length(elements)
    numPresponseT{nu}=pResponse{nu}(1:min(elements));
end
pResponseRS=reshape(cell2mat(numPresponseT'),min(elements),[]);
for p=1:size(pResponseRS,2)
    norpResponseRS(:,p)=mat2gray(pResponseRS(:,p));
end
clear numPresponseT


time=time(1:min(elements),1);   % adjust the length of time to th same as of norpResponseRS


f1=figure;
set(f1,'position',[5 150 1400 800]);
subplot(2,5,1);
plot(norpResponseRS(:,ONorONOFForOFF==1),'color',[0.5 0.5 0.5]);
xlim([0 size(norpResponseRS,1)]);
xlabel('time (sec)','FontSize',20);ylabel('normalized response','FontSize',20);
set(gca,'box','off','color','white');
set(gca,'FontSize',16)
subplot(2,5,2);
plot(norpResponseRS(:,ONorONOFForOFF==2),'color',[0.5 0.5 0.5]);
xlim([0 size(norpResponseRS,1)]);
xlabel('time (sec)','FontSize',20);ylabel('normalized response','FontSize',20);
set(gca,'box','off','color','white');
set(gca,'FontSize',16)

try
subplot(2,5,3);
plot(norpResponseRS(:,ONorONOFForOFF==3),'k');
xlim([0 size(norpResponseRS,1)]);
xlabel('time (sec)','FontSize',20);ylabel('normalized response','FontSize',20);
set(gca,'box','off','color','white');
set(gca,'FontSize',16)
subplot(2,5,4);
plot(norpResponseRS(:,ONorONOFForOFF==4),'k');
xlim([0 size(norpResponseRS,1)]);
xlabel('time (sec)','FontSize',20);ylabel('normalized response','FontSize',20);
set(gca,'box','off','color','white');
set(gca,'FontSize',16)
subplot(2,5,5);
plot(norpResponseRS(:,ONorONOFForOFF==5),'k');
xlim([0 size(norpResponseRS,1)]);
xlabel('time (sec)','FontSize',20);ylabel('normalized response','FontSize',20);
set(gca,'box','off','color','white');
set(gca,'FontSize',16)
end
subplot(2,5,6);
shadedErrorBar(time,mean(norpResponseRS(:,ONorONOFForOFF==1),2),std(norpResponseRS(:,ONorONOFForOFF==1),[],2),'k');
xlim([0 time(end)]); ylim([0 1]);
xlabel('time (sec)','FontSize',20);ylabel('normalized response','FontSize',20);
set(gca,'box','off','color','white');
set(gca,'FontSize',16)
subplot(2,5,7);
shadedErrorBar(time,mean(norpResponseRS(:,ONorONOFForOFF==2),2),std(norpResponseRS(:,ONorONOFForOFF==2),[],2),'k');
xlim([0 time(end)]); ylim([0 1]);
xlabel('time (sec)','FontSize',20);ylabel('normalized response','FontSize',20);
set(gca,'box','off','color','white');
set(gca,'FontSize',16)
subplot(2,5,8);
try
shadedErrorBar(time,mean(norpResponseRS(:,ONorONOFForOFF==3),2),std(norpResponseRS(:,ONorONOFForOFF==3),[],2),'k');
xlim([0 time(end)]); ylim([0 1]);
xlabel('time (sec)','FontSize',20);ylabel('normalized response','FontSize',20);
set(gca,'box','off','color','white');
set(gca,'FontSize',16)
subplot(2,5,9);
shadedErrorBar(time,mean(norpResponseRS(:,ONorONOFForOFF==4),2),std(norpResponseRS(:,ONorONOFForOFF==4),[],2),'k');
xlim([0 time(end)]); ylim([0 1]);
xlabel('time (sec)','FontSize',20);ylabel('normalized response','FontSize',20);
set(gca,'box','off','color','white');
set(gca,'FontSize',16)
subplot(2,5,10);
shadedErrorBar(time,mean(norpResponseRS(:,ONorONOFForOFF==5),2),std(norpResponseRS(:,ONorONOFForOFF==5),[],2),'k');
xlim([0 time(end)]); ylim([0 1]);
xlabel('time (sec)','FontSize',20);ylabel('normalized response','FontSize',20);
set(gca,'box','off','color','white');
set(gca,'FontSize',16)
end
set(gcf,'color','w');





end