
%% load data
load('allTreatments.mat');

%% extract this merged data file from Mar. 16 just to get the time vector
load('E:\Calcium imaging data\clusteringTest\mergedmeta_fov#11');
time=mergedmeta.meanRepDB.ROIdata{1,1}.summary.ordfluotime;
time=time(1:243,1);   % adjust the length of time to be the same as of norpResponseRS

%% extract data for ON and ONOFF cells (based on clustering)
j=1;
k=1;
for i=1:size(predrug.allCells.featureTypeClust,1)
if predrug.allCells.featureTypeClust(i,1)==1
    ON.cellID(j,1)=predrug.allCells.cellID(i,1);
    ON.onset1(j,1)=predrug.allCells.onset1(i,1);
    ON.pResponse{j,1}=predrug.allCells.pResponse{i,1};
    ON.DSI(j,1)=predrug.allCells.DSI(i,1);
    j=j+1;
elseif predrug.allCells.featureTypeClust(i,1)==2
    ONOFF.cellID(k,1)=predrug.allCells.cellID(i,1);
    ONOFF.onset1(k,1)=predrug.allCells.onset1(i,1);
    ONOFF.pResponse{k,1}=predrug.allCells.pResponse{i,1};
    ONOFF.DSI(k,1)=predrug.allCells.DSI(i,1);
    k=k+1;
end
end


%% remove cells 1,2,3,4,7 from the ON cells
% ON.cellID=ON.cellID([1,8,9,12],1);
% ON.onset1=ON.onset1([1,8,9,12],1);
% ON.DSI=ON.DSI([1,8,9,12],1);
% ON.pResponse={ON.pResponse{1};ON.pResponse{8};ON.pResponse{9};ON.pResponse{12}};
% 
% %% remove cells 1,2,3,4,7 from the ONOFF cells
% ONOFF.cellID=ONOFF.cellID([1,14,23,26,27,41,44,46,48,54,55],1);
% ONOFF.onset1=ONOFF.onset1([1,14,23,26,27,41,44,46,48,54,55],1);
% ONOFF.DSI=ONOFF.DSI([1,14,23,26,27,41,44,46,48,54,55],1);
% ONOFF.pResponse={ONOFF.pResponse{1};ONOFF.pResponse{14};ONOFF.pResponse{23};ONOFF.pResponse{26};ONOFF.pResponse{27};ONOFF.pResponse{41};...
%     ONOFF.pResponse{44};ONOFF.pResponse{46};ONOFF.pResponse{48};ONOFF.pResponse{54};ONOFF.pResponse{55}};

%% prepare data 
channel={ON,ONOFF};
for cn=1:size(channel,2)
    pResponse=channel{cn}.pResponse;
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
%DSpResponseNorS=(DSpResponseS-repmat(mean(DSpResponseS(10:40,:),1),size(DSpResponseS,1),1))./(repmat(max(DSpResponseS,[],1),size(DSpResponseS,1),1)-repmat(mean(DSpResponseS(10:40,:),1),size(DSpResponseS,1),1));

% normalize the smoothed responses to start from zero
DSpResponseNorS=(DSpResponseS-repmat(mean(DSpResponseS(10:40,:),1),size(DSpResponseS,1),1));


DSresponse=DSpResponseNorS;

clear DSnorpResponse
clear DSpResponseNorS
clear DSpResponseS

switch cn
    case 1 
        DSresponseON=DSresponse;
    case 2
        DSresponseONOFF=DSresponse;
end


end

%% adjust the DS response for the onset of stimulus for each cell
tmpON=repmat(time,1,size(DSresponseON,2))-repmat(ON.onset1',size(DSresponseON,1),1);
tmpONOFF=repmat(time,1,size(DSresponseONOFF,2))-repmat(ONOFF.onset1',size(DSresponseONOFF,1),1);
meanDSresponseON=mean(DSresponseON,2);
meanDSresponseONOFF=mean(DSresponseONOFF,2);
stdDSresponseON=std(DSresponseON,[],2);
stdDSresponseONOFF=std(DSresponseONOFF,[],2);

% meanDSresponseONS=(meanDSresponseON-repmat(mean(meanDSresponseON(10:40,:),1),size(meanDSresponseON,1),1))./(repmat(max(meanDSresponseON,[],1),size(meanDSresponseON,1),1)-repmat(mean(meanDSresponseON(10:40,:),1),size(meanDSresponseON,1),1));
% meanDSresponseONOFFS=(meanDSresponseONOFF-repmat(mean(meanDSresponseONOFF(10:40,:),1),size(meanDSresponseONOFF,1),1))./(repmat(max(meanDSresponseONOFF,[],1),size(meanDSresponseONOFF,1),1)-repmat(mean(meanDSresponseONOFF(10:40,:),1),size(meanDSresponseONOFF,1),1));

%% plot responses for predrug treatment
figure; 
set(gcf,'position',[450 150 450 800]);
subplot(3,2,1)
shadedErrorBar(mean(tmpON,2),meanDSresponseON,stdDSresponseON,'k');
xlim([min(mean(tmpON,2)) max(mean(tmpON,2))]); ylim([-0.1 1.2]);
subplot(3,2,2)
shadedErrorBar(mean(tmpONOFF,2),meanDSresponseONOFF,stdDSresponseONOFF,'k');
xlim([min(mean(tmpONOFF,2)) max(mean(tmpONOFF,2))]); ylim([-0.1 1.2]);

% plot(mean(tmpON,2),meanDSresponseON);
% plot(mean(tmpONOFF,2),meanDSresponseONOFF);
% set(gcf, 'color', [1 1 1]);
% box off
% ax = gca;
% ax.FontSize=14;
% ylim([-0.1 1]);

%% extract the paired curves for the L-AP-4 treatment for ON cells
j=1;
for i=1:size(ON.cellID,1)
    if strfind(ON.cellID{i},'fov#03')
        target=find(strncmp(LAP4.allCells.cellID,[ON.cellID{i}(1:12),'fov#09'],numel([ON.cellID{i}(1:12),'fov#09']))==1);
    elseif strfind(ON.cellID{i},'fov#04')
        target=find(strncmp(LAP4.allCells.cellID,[ON.cellID{i}(1:12),'fov#10'],numel([ON.cellID{i}(1:12),'fov#10']))==1);
    elseif strfind(ON.cellID{i},'fov#05')
        target=find(strncmp(LAP4.allCells.cellID,[ON.cellID{i}(1:12),'fov#11'],numel([ON.cellID{i}(1:12),'fov#11']))==1);
    elseif strfind(ON.cellID{i},'fov#06')
        target=find(strncmp(LAP4.allCells.cellID,[ON.cellID{i}(1:12),'fov#12'],numel([ON.cellID{i}(1:12),'fov#12']))==1);
    elseif strfind(ON.cellID{i},'fov#07')
        target=find(strncmp(LAP4.allCells.cellID,[ON.cellID{i}(1:12),'fov#13'],numel([ON.cellID{i}(1:12),'fov#13']))==1);
    elseif strfind(ON.cellID{i},'fov#08')
        target=find(strncmp(LAP4.allCells.cellID,[ON.cellID{i}(1:12),'fov#14'],numel([ON.cellID{i}(1:12),'fov#14']))==1);
    end
    for b=1:size(target,1)
        if strfind(ON.cellID{i}(28+1:end),LAP4.allCells.cellID{target(b)}(28+1:end))
            ON.LAP4.cellID(j,1)=LAP4.allCells.cellID(target(b),1);
            ON.LAP4.onset1(j,1)=LAP4.allCells.onset1(target(b),1);
            ON.LAP4.DSI(j,1)=LAP4.allCells.DSI(target(b),1);
            ON.LAP4.pResponse{j,1}=LAP4.allCells.pResponse{target(b),1};
        end
    end
    j=j+1;
end

%% extract the paired curves for the L-AP-4 treatment for ONOFF cells
j=1;
for i=1:size(ONOFF.cellID,1)
    if strfind(ONOFF.cellID{i},'fov#03')
        target=find(strncmp(LAP4.allCells.cellID,[ONOFF.cellID{i}(1:12),'fov#09'],numel([ONOFF.cellID{i}(1:12),'fov#09']))==1);
    elseif strfind(ONOFF.cellID{i},'fov#04')
        target=find(strncmp(LAP4.allCells.cellID,[ONOFF.cellID{i}(1:12),'fov#10'],numel([ONOFF.cellID{i}(1:12),'fov#10']))==1);
    elseif strfind(ONOFF.cellID{i},'fov#05')
        target=find(strncmp(LAP4.allCells.cellID,[ONOFF.cellID{i}(1:12),'fov#11'],numel([ONOFF.cellID{i}(1:12),'fov#11']))==1);
    elseif strfind(ONOFF.cellID{i},'fov#06')
        target=find(strncmp(LAP4.allCells.cellID,[ONOFF.cellID{i}(1:12),'fov#12'],numel([ONOFF.cellID{i}(1:12),'fov#12']))==1);
    elseif strfind(ONOFF.cellID{i},'fov#07')
        target=find(strncmp(LAP4.allCells.cellID,[ONOFF.cellID{i}(1:12),'fov#13'],numel([ONOFF.cellID{i}(1:12),'fov#13']))==1);
    elseif strfind(ONOFF.cellID{i},'fov#08')
        target=find(strncmp(LAP4.allCells.cellID,[ONOFF.cellID{i}(1:12),'fov#14'],numel([ONOFF.cellID{i}(1:12),'fov#14']))==1);
    end
    for b=1:size(target,1)
        if strfind(ONOFF.cellID{i}(28+1:end),LAP4.allCells.cellID{target(b)}(28+1:end))
            ONOFF.LAP4.cellID(j,1)=LAP4.allCells.cellID(target(b),1);
            ONOFF.LAP4.onset1(j,1)=LAP4.allCells.onset1(target(b),1);
            ONOFF.LAP4.DSI(j,1)=LAP4.allCells.DSI(target(b),1);
            ONOFF.LAP4.pResponse{j,1}=LAP4.allCells.pResponse{target(b),1};
        end
    end
    j=j+1;
end

%% prepare data of LAP4
channel={ON.LAP4,ONOFF.LAP4};
for cn=1:size(channel,2)
    pResponse=channel{cn}.pResponse;
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

% normalize the smoothed responses to start from zero
DSpResponseNorS=(DSpResponseS-repmat(mean(DSpResponseS(10:40,:),1),size(DSpResponseS,1),1));

DSresponse=DSpResponseNorS;

clear DSnorpResponse
clear DSpResponseNorS
clear DSpResponseS

switch cn
    case 1 
        DSresponseONLAP4=DSresponse;
    case 2
        DSresponseONOFFLAP4=DSresponse;
end
end

%% adjust the DS response for the onset of stimulus for each cell for LAP4
tmpONLAP4=repmat(time,1,size(DSresponseONLAP4,2))-repmat(ON.LAP4.onset1',size(DSresponseONLAP4,1),1);
tmpONOFFLAP4=repmat(time,1,size(DSresponseONOFFLAP4,2))-repmat(ONOFF.LAP4.onset1',size(DSresponseONOFFLAP4,1),1);
meanDSresponseONLAP4=mean(DSresponseONLAP4,2);
meanDSresponseONOFFLAP4=mean(DSresponseONOFFLAP4,2);
stdDSresponseONLAP4=std(DSresponseONLAP4,[],2);
stdDSresponseONOFFLAP4=std(DSresponseONOFFLAP4,[],2);

% meanDSresponseONS=(meanDSresponseON-repmat(mean(meanDSresponseON(10:40,:),1),size(meanDSresponseON,1),1))./(repmat(max(meanDSresponseON,[],1),size(meanDSresponseON,1),1)-repmat(mean(meanDSresponseON(10:40,:),1),size(meanDSresponseON,1),1));
% meanDSresponseONOFFS=(meanDSresponseONOFF-repmat(mean(meanDSresponseONOFF(10:40,:),1),size(meanDSresponseONOFF,1),1))./(repmat(max(meanDSresponseONOFF,[],1),size(meanDSresponseONOFF,1),1)-repmat(mean(meanDSresponseONOFF(10:40,:),1),size(meanDSresponseONOFF,1),1));

%% plot responses for predrug treatment 
subplot(3,2,3)
shadedErrorBar(mean(tmpONLAP4,2),meanDSresponseONLAP4,stdDSresponseONLAP4,'k');
xlim([min(mean(tmpONLAP4,2)) max(mean(tmpONLAP4,2))]); ylim([-0.1 1]);
subplot(3,2,4)
shadedErrorBar(mean(tmpONOFFLAP4,2),meanDSresponseONOFFLAP4,stdDSresponseONOFFLAP4,'k');
xlim([min(mean(tmpONOFFLAP4,2)) max(mean(tmpONOFFLAP4,2))]); ylim([-0.1 1]);

%% extract the paired curves for the ACET treatment for ON cells
j=1;
for i=1:size(ON.cellID,1)
    if strfind(ON.cellID{i},'fov#03')
        target=find(strncmp(ACET.allCells.cellID,[ON.cellID{i}(1:12),'fov#18'],numel([ON.cellID{i}(1:12),'fov#18']))==1);
    elseif strfind(ON.cellID{i},'fov#05')
        target=find(strncmp(ACET.allCells.cellID,[ON.cellID{i}(1:12),'fov#17'],numel([ON.cellID{i}(1:12),'fov#17']))==1);
    elseif strfind(ON.cellID{i},'fov#06')
        target=find(strncmp(ACET.allCells.cellID,[ON.cellID{i}(1:12),'fov#16'],numel([ON.cellID{i}(1:12),'fov#16']))==1);
    elseif strfind(ON.cellID{i},'fov#08')
        target=find(strncmp(ACET.allCells.cellID,[ON.cellID{i}(1:12),'fov#15'],numel([ON.cellID{i}(1:12),'fov#15']))==1);
    end
    for b=1:size(target,1)
        if strfind(ON.cellID{i}(28+1:end),ACET.allCells.cellID{target(b)}(28+1:end))
            ON.ACET.cellID(j,1)=ACET.allCells.cellID(target(b),1);
            ON.ACET.onset1(j,1)=ACET.allCells.onset1(target(b),1);
            ON.ACET.DSI(j,1)=ACET.allCells.DSI(target(b),1);
            ON.ACET.pResponse{j,1}=ACET.allCells.pResponse{target(b),1};
        end
    end
    j=j+1;
end

%% extract the paired curves for the L-AP-4 treatment for ONOFF cells
j=1;
for i=1:size(ONOFF.cellID,1)
    if strfind(ONOFF.cellID{i},'fov#03')
        target=find(strncmp(ACET.allCells.cellID,[ONOFF.cellID{i}(1:12),'fov#18'],numel([ONOFF.cellID{i}(1:12),'fov#18']))==1);
    elseif strfind(ONOFF.cellID{i},'fov#05')
        target=find(strncmp(ACET.allCells.cellID,[ONOFF.cellID{i}(1:12),'fov#17'],numel([ONOFF.cellID{i}(1:12),'fov#17']))==1);
    elseif strfind(ONOFF.cellID{i},'fov#06')
        target=find(strncmp(ACET.allCells.cellID,[ONOFF.cellID{i}(1:12),'fov#16'],numel([ONOFF.cellID{i}(1:12),'fov#16']))==1);
    elseif strfind(ONOFF.cellID{i},'fov#08')
        target=find(strncmp(ACET.allCells.cellID,[ONOFF.cellID{i}(1:12),'fov#15'],numel([ONOFF.cellID{i}(1:12),'fov#15']))==1);
    end
    for b=1:size(target,1)
        if strfind(ONOFF.cellID{i}(28+1:end),ACET.allCells.cellID{target(b)}(28+1:end))
            ONOFF.ACET.cellID(j,1)=ACET.allCells.cellID(target(b),1);
            ONOFF.ACET.onset1(j,1)=ACET.allCells.onset1(target(b),1);
            ONOFF.ACET.DSI(j,1)=ACET.allCells.DSI(target(b),1);
            ONOFF.ACET.pResponse{j,1}=ACET.allCells.pResponse{target(b),1};
        end
    end
    j=j+1;
end

%% prepare data of ACET
channel={ON.ACET,ONOFF.ACET};
for cn=1:size(channel,2)
    pResponse=channel{cn}.pResponse;
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

% normalize the smoothed responses to start from zero
DSpResponseNorS=(DSpResponseS-repmat(mean(DSpResponseS(10:40,:),1),size(DSpResponseS,1),1));

DSresponse=DSpResponseNorS;

clear DSnorpResponse
clear DSpResponseNorS
clear DSpResponseS

switch cn
    case 1 
        DSresponseONACET=DSresponse;
    case 2
        DSresponseONOFFACET=DSresponse;
end
end

%% adjust the DS response for the onset of stimulus for each cell for ACET
tmpONACET=repmat(time,1,size(DSresponseONACET,2))-repmat(ON.ACET.onset1',size(DSresponseONACET,1),1);
tmpONOFFACET=repmat(time,1,size(DSresponseONOFFACET,2))-repmat(ONOFF.ACET.onset1',size(DSresponseONOFFACET,1),1);
meanDSresponseONACET=mean(DSresponseONACET,2);
meanDSresponseONOFFACET=mean(DSresponseONOFFACET,2);
stdDSresponseONACET=std(DSresponseONACET,[],2);
stdDSresponseONOFFACET=std(DSresponseONOFFACET,[],2);

% meanDSresponseONS=(meanDSresponseON-repmat(mean(meanDSresponseON(10:40,:),1),size(meanDSresponseON,1),1))./(repmat(max(meanDSresponseON,[],1),size(meanDSresponseON,1),1)-repmat(mean(meanDSresponseON(10:40,:),1),size(meanDSresponseON,1),1));
% meanDSresponseONOFFS=(meanDSresponseONOFF-repmat(mean(meanDSresponseONOFF(10:40,:),1),size(meanDSresponseONOFF,1),1))./(repmat(max(meanDSresponseONOFF,[],1),size(meanDSresponseONOFF,1),1)-repmat(mean(meanDSresponseONOFF(10:40,:),1),size(meanDSresponseONOFF,1),1));

%% plot responses for predrug treatment 
subplot(3,2,5)
shadedErrorBar(mean(tmpONACET,2),meanDSresponseONACET,stdDSresponseONACET,'k');
xlim([min(mean(tmpONACET,2)) max(mean(tmpONACET,2))]); ylim([-0.1 1]);
subplot(3,2,6)
shadedErrorBar(mean(tmpONOFFACET,2),meanDSresponseONOFFACET,stdDSresponseONOFFACET,'k');
xlim([min(mean(tmpONOFFACET,2)) max(mean(tmpONOFFACET,2))]); ylim([-0.1 1]);

%% plot example ON cell(cell# 4) and ONOFF cell (cell# 23 )
n=12;
m=23;   %45,27, 

figure; 
set(gcf,'position',[450 150 450 800]);
subplot(3,2,1)
plot(tmpON(:,n),DSresponseON(:,n));
xlim([min(mean(tmpON,2)) max(mean(tmpON,2))]); ylim([-0.1 1]);
subplot(3,2,2)
plot(tmpONOFF(:,m),DSresponseONOFF(:,m));
xlim([min(mean(tmpONOFF,2)) max(mean(tmpONOFF,2))]); ylim([-0.1 1]);

subplot(3,2,3)
plot(tmpONLAP4(:,n),DSresponseONLAP4(:,n));
xlim([min(mean(tmpONLAP4,2)) max(mean(tmpONLAP4,2))]); ylim([-0.1 1]);
subplot(3,2,4)
plot(tmpONOFFLAP4(:,m),DSresponseONOFFLAP4(:,m));
xlim([min(mean(tmpONOFFLAP4,2)) max(mean(tmpONOFFLAP4,2))]); ylim([-0.1 1]);

subplot(3,2,5)
plot(tmpONACET(:,n),DSresponseONACET(:,n));
xlim([min(mean(tmpONACET,2)) max(mean(tmpONACET,2))]); ylim([-0.1 1]);
subplot(3,2,6)
plot(tmpONOFFACET(:,m),DSresponseONOFFACET(:,m));
xlim([min(mean(tmpONOFFACET,2)) max(mean(tmpONOFFACET,2))]); ylim([-0.1 1]);









