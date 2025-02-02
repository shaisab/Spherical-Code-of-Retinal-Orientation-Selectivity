% clear all
% clc

load('pooledMapAlphaCorr_OS data_40_OSI0.30_1.00_SI_0.00_created_Dec_11_2022_11_41_On_Sus.mat');
pResponse=pooledmap.pResponse;
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

Non_normalized_On_Sus=DSpResponse;
Normalized_On_Sus=DSnorpResponse;
Shifted_meanr_On_Sus=pooledmap.meanrs;

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

onset1=mean(pooledmap.onset1);
onset2=onset1+(c/speedterm);
offset1=onset1+(mergedmeta.metadataDB{1,1}.log.width/(corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed));
offset2=onset2+(mergedmeta.metadataDB{1,1}.log.width/(corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed));
%% on cells 
f13=figure('Color',[1 1 1],'Renderer','painters');
set(f13,'position',[5 150 1300 400]);

subplot(2,2,1,'Parent',f13,'YColor',[1 1 1],'XColor',[1 1 1]);
title('On Sustained','Interpreter','none');

subplot(2,2,2,'Parent',f13,'YColor',[1 1 1],'XColor',[1 1 1]);
shadedErrorBar(time,mean(DSnorpResponse,2),std(DSnorpResponse,[],2),'k');
hold on
plot(time,mean(DSnorpResponse,2),'b','LineWidth',1);
xline(onset1);
xline(offset2);
xlim([0 time(end)]); ylim([0 1]);
xlabel('Time (sec)','FontSize',16);
ylabel('\DeltaF/F','FontSize',16);
title(['\itn\rm = ',num2str(size(DSnorpResponse,2)/2)],'FontSize',16);
box off

load('pooledMapAlphaCorr_OS data_40_OSI0.30_1.00_SI_0.00_created_Dec_11_2022_12_06_On_Trans.mat');
pResponse=pooledmap.pResponse;
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
clear DSnorpResponse
for p=1:size(DSpResponse,2)
    DSnorpResponse(:,p)=mat2gray(DSpResponse(:,p));
end

Non_normalized_On_Trans=DSpResponse;
Normalized_On_Trans=DSnorpResponse;
Shifted_meanr_On_Trans=pooledmap.meanrs;

time=mergedmeta.meanRepDB.ROIdata{1,1}.summary.ordfluotime;
time=time(1:min(elements),1);   % adjust the length of time to be the same as of norpResponseRS

onset1=mean(pooledmap.onset1);
onset2=onset1+(c/speedterm);
offset1=onset1+(mergedmeta.metadataDB{1,1}.log.width/(corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed));
offset2=onset2+(mergedmeta.metadataDB{1,1}.log.width/(corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed));

subplot(2,2,3,'Parent',f13,'YColor',[1 1 1],'XColor',[1 1 1]);
title('On Transient','Interpreter','none');

subplot(2,2,4,'Parent',f13,'YColor',[1 1 1],'XColor',[1 1 1]);
shadedErrorBar(time,mean(DSnorpResponse,2),std(DSnorpResponse,[],2),'k');
hold on
plot(time,mean(DSnorpResponse,2),'b','LineWidth',1);
xline(onset1);
xline(offset2);
xlim([0 time(end)]); ylim([0 1]);
xlabel('Time (sec)','FontSize',16);
ylabel('\DeltaF/F','FontSize',16);
%axis off
title(['\itn\rm = ',num2str(size(DSnorpResponse,2)/2)],'FontSize',16);
box off


sgtitle('On Cells');