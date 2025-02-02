clear all
clc

load('allCellsAlphaCorr_OS data_40_PPC_0.9_OSI0.20_1.00_created_Aug_09_2023_12_28.mat');
pResponse=allCells.pResponse(allCells.OSI>0.4);
% determine the minimum number of elements
elements=[];
for nu=1:numel(pResponse) %nu=1:numel(pResponse/2)
    elements=[elements,numel(pResponse{nu})];
end

% standardize the length of responses
for nu=1:numel(pResponse) %nu=1:numel(pResponse/2)
    numPresponseT{nu}=pResponse{nu}(1:min(elements));
end
DSpResponse=reshape(cell2mat(numPresponseT'),min(elements),[]);
clear numPresponseT

% normalize the responses
for p=1:size(DSpResponse,2)
    DSnorpResponse(:,p)=mat2gray(DSpResponse(:,p));
end

%% extract this merged data file from Mar. 16 just to get the time vector
load('clusteringTest\mergedmeta_fov#11');
time=mergedmeta.meanRepDB.ROIdata{1,1}.summary.ordfluotime;
time=time(1:min(elements),1);   % adjust the length of time to be the same as of norpResponseRS


z=750;      % at non-right angles, if the bar start to move within the screnn (rather than outside of it, as is the case for right angles).                     
screenWidthum=3500;  % 3 mm at the plane of the retina
pixelprojwidth=screenWidthum/600;   % 3 mm / 600 pixels
radDendTreeum=300;
corrSpeedFactor=1.2;
delay=0;
                              
speedterm=pixelprojwidth*corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed;
c=radDendTreeum;

%% on-off cells 
f13=figure('Color',[1 1 1],'Renderer','painters');
set(f13,'position',[5 150 1300 400]);

subplot(2,2,1,'Parent',f13,'YColor',[1 1 1],'XColor',[1 1 1]);
title('On-Off Sustained','Interpreter','none');

onset1=mean(allCells.onset1(allCells.typeClust(allCells.OSI>0.4)==6));
onset2=onset1+(c/speedterm);
offset1=onset1+(mergedmeta.metadataDB{1,1}.log.width/(corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed));
offset2=onset2+(mergedmeta.metadataDB{1,1}.log.width/(corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed));

subplot(2,2,2,'Parent',f13,'YColor',[1 1 1],'XColor',[1 1 1]);
shadedErrorBar(time,mean(DSnorpResponse(:,allCells.typeClust(allCells.OSI>0.4)==6),2),std(DSnorpResponse(:,allCells.typeClust(allCells.OSI>0.4)==6),[],2),'k');
hold on
plot(time,mean(DSnorpResponse(:,allCells.typeClust(allCells.OSI>0.4)==6),2),'b','LineWidth',1);
xline(onset1);
xline(offset2);
xlim([0 time(end)]); ylim([0 1]);
xlabel('Time (sec)','FontSize',16);
ylabel('\DeltaF/F','FontSize',16);
title(['\itn\rm = ',num2str(sum(allCells.typeClust(allCells.OSI>0.4)==6)/2)],'FontSize',16);
box off

Non_normalized_On_Off_Sus=DSpResponse(:,allCells.typeClust(allCells.OSI>0.4)==6);
normalized_On_Off_Sus=DSnorpResponse(:,allCells.typeClust(allCells.OSI>0.4)==6);
mean_Non_normalized_On_Off_Sus=mean(DSpResponse(:,allCells.typeClust(allCells.OSI>0.4)==6));
mean_normalized_On_Off_Sus=mean(DSnorpResponse(:,allCells.typeClust(allCells.OSI>0.4)==6));
stderr_Non_normalized_On_Off_Sus=std(DSpResponse(:,allCells.typeClust(allCells.OSI>0.4)==6))/sqrt(length(DSpResponse(:,allCells.typeClust(allCells.OSI>0.4)==6)));
stderr_normalized_On_Off_Sus=std(DSnorpResponse(:,allCells.typeClust(allCells.OSI>0.4)==6))/sqrt(length(DSpResponse(:,allCells.typeClust(allCells.OSI>0.4)==6)));

subplot(2,2,3,'Parent',f13,'YColor',[1 1 1],'XColor',[1 1 1]);
title('On-Off Transient','Interpreter','none');

onset1=mean(allCells.onset1(allCells.typeClust(allCells.OSI>0.4)==3));
onset2=onset1+(c/speedterm);
offset1=onset1+(mergedmeta.metadataDB{1,1}.log.width/(corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed));
offset2=onset2+(mergedmeta.metadataDB{1,1}.log.width/(corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed));

subplot(2,2,4,'Parent',f13,'YColor',[1 1 1],'XColor',[1 1 1]);
shadedErrorBar(time,mean(DSnorpResponse(:,allCells.typeClust(allCells.OSI>0.4)==3),2),std(DSnorpResponse(:,allCells.typeClust(allCells.OSI>0.4)==3),[],2),'k');
hold on
plot(time,mean(DSnorpResponse(:,allCells.typeClust(allCells.OSI>0.4)==3),2),'b','LineWidth',1);
xline(onset1);
xline(offset2);
xlim([0 time(end)]); ylim([0 1]);
xlabel('Time (sec)','FontSize',16);
ylabel('\DeltaF/F','FontSize',16);
%axis off
title(['\itn\rm = ',num2str(sum(allCells.typeClust(allCells.OSI>0.4)==3)/2)],'FontSize',16);
box off

Non_normalized_On_Off_trans=DSpResponse(:,allCells.typeClust(allCells.OSI>0.4)==3);
normalized_On_Off_trans=DSnorpResponse(:,allCells.typeClust(allCells.OSI>0.4)==3);
mean_Non_normalized_On_Off_trans=mean(DSpResponse(:,allCells.typeClust(allCells.OSI>0.4)==3));
mean_normalized_On_Off_trans=mean(DSnorpResponse(:,allCells.typeClust(allCells.OSI>0.4)==3));
stderr_Non_normalized_On_Off_trans=std(DSpResponse(:,allCells.typeClust(allCells.OSI>0.4)==3))/sqrt(length(DSpResponse(:,allCells.typeClust(allCells.OSI>0.4)==3)));
stderr_normalized_On_Off_trans=std(DSnorpResponse(:,allCells.typeClust(allCells.OSI>0.4)==3))/sqrt(length(DSpResponse(:,allCells.typeClust(allCells.OSI>0.4)==3)));

sgtitle('On-Off Cells');

%on cells
f14=figure('Color',[1 1 1],'Renderer','painters');
set(f14,'position',[5 150 1300 400]);

subplot(2,2,1,'Parent',f14,'YColor',[1 1 1],'XColor',[1 1 1]);
title('On Sustained','Interpreter','none');

onset1=mean(allCells.onset1(allCells.typeClust(allCells.OSI>0.4)==2));
onset2=onset1+(c/speedterm);
offset1=onset1+(mergedmeta.metadataDB{1,1}.log.width/(corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed));
offset2=onset2+(mergedmeta.metadataDB{1,1}.log.width/(corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed));

subplot(2,2,2,'Parent',f14,'YColor',[1 1 1],'XColor',[1 1 1]);
shadedErrorBar(time,mean(DSnorpResponse(:,allCells.typeClust(allCells.OSI>0.4)==2),2),std(DSnorpResponse(:,allCells.typeClust(allCells.OSI>0.4)==2),[],2),'k');
hold on
plot(time,mean(DSnorpResponse(:,allCells.typeClust(allCells.OSI>0.4)==2),2),'b','LineWidth',1);
xline(onset1);
xline(offset2);
xlim([0 time(end)]); ylim([0 1]);
xlabel('Time (sec)','FontSize',16);
ylabel('\DeltaF/F','FontSize',16);
%axis off
title(['\itn\rm = ',num2str(sum(allCells.typeClust(allCells.OSI>0.4)==2)/2)],'FontSize',16);
box off

Non_normalized_On_sus=DSpResponse(:,allCells.typeClust(allCells.OSI>0.4)==2);
normalized_On_sus=DSnorpResponse(:,allCells.typeClust(allCells.OSI>0.4)==2);
mean_Non_normalized_On_sus=mean(DSpResponse(:,allCells.typeClust(allCells.OSI>0.4)==2));
mean_normalized_On_sus=mean(DSnorpResponse(:,allCells.typeClust(allCells.OSI>0.4)==2));
stderr_Non_normalized_On_sus=std(DSpResponse(:,allCells.typeClust(allCells.OSI>0.4)==2))/sqrt(length(DSpResponse(:,allCells.typeClust(allCells.OSI>0.4)==2)));
stderr_normalized_On_sus=std(DSnorpResponse(:,allCells.typeClust(allCells.OSI>0.4)==2))/sqrt(length(DSpResponse(:,allCells.typeClust(allCells.OSI>0.4)==2)));

subplot(2,2,3,'Parent',f14,'YColor',[1 1 1],'XColor',[1 1 1]);
title('On Transient','Interpreter','none');

onset1=mean(allCells.onset1(allCells.typeClust(allCells.OSI>0.4)==5));
onset2=onset1+(c/speedterm);
offset1=onset1+(mergedmeta.metadataDB{1,1}.log.width/(corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed));
offset2=onset2+(mergedmeta.metadataDB{1,1}.log.width/(corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed));

subplot(2,2,4,'Parent',f14,'YColor',[1 1 1],'XColor',[1 1 1]);
shadedErrorBar(time,mean(DSnorpResponse(:,allCells.typeClust(allCells.OSI>0.4)==5),2),std(DSnorpResponse(:,allCells.typeClust(allCells.OSI>0.4)==5),[],2),'k');
hold on
plot(time,mean(DSnorpResponse(:,allCells.typeClust(allCells.OSI>0.4)==5),2),'b','LineWidth',1);
xline(onset1);
xline(offset2);
xlim([0 time(end)]); ylim([0 1]);
xlabel('Time (sec)','FontSize',16);
ylabel('\DeltaF/F','FontSize',16);
%axis off
title(['\itn\rm = ',num2str(sum(allCells.typeClust(allCells.OSI>0.4)==5)/2)],'FontSize',16);
box off

Non_normalized_On_trans=DSpResponse(:,allCells.typeClust(allCells.OSI>0.4)==5);
normalized_On_trans=DSnorpResponse(:,allCells.typeClust(allCells.OSI>0.4)==5);
mean_Non_normalized_On_trans=mean(DSpResponse(:,allCells.typeClust(allCells.OSI>0.4)==5));
mean_normalized_On_trans=mean(DSnorpResponse(:,allCells.typeClust(allCells.OSI>0.4)==5));
stderr_Non_normalized_On_trans=std(DSpResponse(:,allCells.typeClust(allCells.OSI>0.4)==5))/sqrt(length(DSpResponse(:,allCells.typeClust(allCells.OSI>0.4)==5)));
stderr_normalized_On_trans=std(DSnorpResponse(:,allCells.typeClust(allCells.OSI>0.4)==5))/sqrt(length(DSpResponse(:,allCells.typeClust(allCells.OSI>0.4)==5)));

sgtitle('On Cells');
