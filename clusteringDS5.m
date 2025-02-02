function [idxT,inddataT]=clusteringDS5(pResponse,onset1,NORMALRESstate,postPcut)

%% define parameters

slopePoints=5;

%% prepare data 
% determine the minimum number of elements
elements=[];
for nu=1:numel(pResponse)
    elements=[elements,numel(pResponse{nu})];
end

% standardize the length of responses
for nu=1:numel(pResponse)
    numPresponseT{nu}=pResponse{nu}(1:min(elements));
end
DSpResponse=reshape(cell2mat(numPresponseT'),min(elements),[]);
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

% switch between using the normalized responses or not
switch NORMALRESstate
    case 1
        DSresponse=DSpResponseNorS;
    case 0
        DSresponse=DSpResponseS;
end

%% extract this merged data file from Mar. 16 just to get the time vector
load('clusteringTest\mergedmeta_fov#11');
time=mergedmeta.meanRepDB.ROIdata{1,1}.summary.ordfluotime;
time=time(1:min(elements),1);   % adjust the length of time to be the same as of norpResponseRS

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
%     onind=find(peakLocMax>60);
%     if isempty(onind)
%         DSpeakLocOn(p)=NaN;
%     else
%         DSpeakLocOn(p)=peakLocMax(onind(1),1);
%     end
    
    
    %%% find ON peak
    [peakLocMaxON,peakMagMaxON]=peakfinder(DSresponse(1:150,p),0.1,0.1,1,false);
    onind=find(peakLocMaxON>60);
    if isempty(onind)
        [peakLocMaxON,peakMagMaxON]=peakfinder(DSresponse(1:155,p),0.01,0.01,1,false);
        onind=find(peakLocMaxON>60);
        if isempty(onind) 
        DSpeakLocOn(p)=130;
        else
            DSpeakLocOn(p)=peakLocMaxON(onind(end),1);
        end
    else
        DSpeakLocOn(p)=peakLocMaxON(onind(1),1);
    end
    
    %%% calculate slope of ON peak
    if ~isnan(DSpeakLocOn(p))
        yOnData=DSresponse(DSpeakLocOn(p):DSpeakLocOn(p)+slopePoints,p);
        % Fit model to data.
        ft = fittype('poly1');
        opts = fitoptions(ft);
        [fitresultOn, gofOn] = fit((1:size(yOnData,1))', yOnData(:,1), ft, opts );
        DSslopeOn(p)=fitresultOn.p1;
    else
        DSslopeOn(p)=NaN;
    end
end

DSonMeanResponse=mean(DSresponse(60:130,:),1);
DSoffMeanResponse=mean(DSresponse(158:240,:),1);

DSpeakLocOn=time(DSpeakLocOn);      % convert DSpeakLocOn from number of data points to time
DSpeakLocOn=DSpeakLocOn-onset1;    % adjust DSpeakLocOn relative to the onset of the stimulus



% DSdiff=diff(DSresponse);
% [~,diffmax]=max(DSresponse(1:80,:));
% DSabsPeakLockOn=DSpeakLocOn-diffmax;

DSpeakLocOn=DSpeakLocOn;
DSpeakLocOff=DSpeakLocOff';
DSonMeanResponse=DSonMeanResponse';
DSoffMeanResponse=DSoffMeanResponse';
DSpeakMinMag=DSpeakMinMag';
DSslopeOn=DSslopeOn';
DSslopeOff=DSslopeOff';

data=[DSpeakLocOn,DSslopeOn,DSonMeanResponse,DSoffMeanResponse,DSpeakMinMag];
%data=[DSpeakLocOn,DSslopeOn,DSonMeanResponse,DSoffMeanResponse];


figure;histogram(data(:,1),100);
xlabel('time to ON peak (data points)');

figure;histogram(data(:,2),100);
xlabel('slope following ON peak)');

figure;histogram(data(:,3),100);
xlabel('integral of ON response');

figure;histogram(data(:,4),100);
xlabel('integral of OFF response');

figure;histogram(data(:,5),100);
xlabel('Minimum of ON response');


%% calculate the Gaussian mixture models with 2 to 10 clusters, 50 repetitions (every time it uses a different initial value), 
% and determine the AIC, BIC, and negative log likelihood
options=statset('MaxIter',100000,'TolFun',1e-15);
for nK=2:1:10
GM=fitgmdist(data,nK,'RegularizationValue',0.001,'Replicates',50,'Options',options)
idx(:,nK)=cluster(GM,data);
P{nK}=posterior(GM,data);
allLL(nK)=GM.NegativeLogLikelihood;
allBIC(nK)=GM.BIC;
allAIC(nK)=GM.AIC;
end

% calculates the min BIC, delta BIC, and the weight of delta BIC(the probability that the
% candidate model is the best among the set of candidate models). The model with the highest weight is the the best model
% based on a paper named 'AIC and BIC weights' found in the folder 'Analysis cone'.
deltaBIC=allBIC-min(allBIC);
ExpDeltaBIC=exp(-0.5.*deltaBIC);
weightDeltaBIC=ExpDeltaBIC./sum(ExpDeltaBIC);
[~,OptimalK]=min(allBIC(1:end))                                        
[~,bestEvidenceBIC]=max(weightDeltaBIC);

% calculates the min AIC, delta AIC, and the weight of delta AIC(the probability that the
% candidate model is the best among the set of candidate models). The model with the highest weight is the the best model
% based on a paper named 'AIC and BIC weights' found in the folder 'Analysis cone'.
% deltaAIC=allAIC-min(allAIC);
% ExpDeltaAIC=exp(-0.5.*deltaAIC);
% weightDeltaAIC=ExpDeltaAIC./sum(ExpDeltaAIC);   
% [~,OptimalK]=min(allAIC(1:end))                                        
% [~,bestEvidenceAIC]=max(weightDeltaAIC);

figure;
plot(2:10,allBIC(2:10),'-ob');

%OptimalK=5;

% recalculate the Gaussian mixture models with the optimal number of
% clusters (optimalK)
%GM=fitgmdist(data,OptimalK,'RegularizationValue',0.001,'Replicates',50,'Options',options)
GM=fitgmdist(data,OptimalK,'RegularizationValue',0.001,'Replicates',50,'Options',options)
idx=cluster(GM,data);
post=posterior(GM,data);        %calculates the posterior probability of cells


% calculate the number of cells in each cluster
idx1=sum(idx==1);
idx2=sum(idx==2);
idx3=sum(idx==3);
idx4=sum(idx==4);
idx5=sum(idx==5);
idx6=sum(idx==6);
idx7=sum(idx==7);
idx8=sum(idx==8);
idx9=sum(idx==9);
idx10=sum(idx==10);

f11=figure;
set(f11,'position',[5 150 950 800]);
subplot(2,10,1);
plot(DSnorpResponse(:,idx==1),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,10,2);
plot(DSnorpResponse(:,idx==2),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,10,3);
plot(DSnorpResponse(:,idx==3),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,10,4);
plot(DSnorpResponse(:,idx==4),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,10,5);
plot(DSnorpResponse(:,idx==5),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,10,6);
plot(DSnorpResponse(:,idx==6),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,10,7);
plot(DSnorpResponse(:,idx==7),'k');
xlim([0 size(DSnorpResponse,1)]);

subplot(2,10,8);
plot(DSnorpResponse(:,idx==8),'k');
xlim([0 size(DSnorpResponse,1)]);

subplot(2,10,9);
plot(DSnorpResponse(:,idx==9),'k');
xlim([0 size(DSnorpResponse,1)]);

subplot(2,10,10);
plot(DSnorpResponse(:,idx==10),'k');
xlim([0 size(DSnorpResponse,1)]);


subplot(2,10,11);
shadedErrorBar(1:size(DSnorpResponse,1),mean(DSnorpResponse(:,idx==1),2),std(DSnorpResponse(:,idx==1),[],2),'k');
xlim([0 size(DSnorpResponse,1)]); ylim([0 1]);
subplot(2,10,12);
shadedErrorBar(1:size(DSnorpResponse,1),mean(DSnorpResponse(:,idx==2),2),std(DSnorpResponse(:,idx==2),[],2),'k');
xlim([0 size(DSnorpResponse,1)]); ylim([0 1]);
subplot(2,10,13);
shadedErrorBar(1:size(DSnorpResponse,1),mean(DSnorpResponse(:,idx==3),2),std(DSnorpResponse(:,idx==3),[],2),'k');
xlim([0 size(DSnorpResponse,1)]); ylim([0 1]);
subplot(2,10,14);
shadedErrorBar(1:size(DSnorpResponse,1),mean(DSnorpResponse(:,idx==4),2),std(DSnorpResponse(:,idx==4),[],2),'k');
xlim([0 size(DSnorpResponse,1)]); ylim([0 1]);
subplot(2,10,15);
shadedErrorBar(1:size(DSnorpResponse,1),mean(DSnorpResponse(:,idx==5),2),std(DSnorpResponse(:,idx==5),[],2),'k');
xlim([0 size(DSnorpResponse,1)]); ylim([0 1]);
subplot(2,10,16);
shadedErrorBar(1:size(DSnorpResponse,1),mean(DSnorpResponse(:,idx==6),2),std(DSnorpResponse(:,idx==6),[],2),'k');
xlim([0 size(DSnorpResponse,1)]); ylim([0 1]);
subplot(2,10,17);
shadedErrorBar(1:size(DSnorpResponse,1),mean(DSnorpResponse(:,idx==7),2),std(DSnorpResponse(:,idx==7),[],2),'k');
xlim([0 size(DSnorpResponse,1)]); ylim([0 1]);
subplot(2,10,18);
shadedErrorBar(1:size(DSnorpResponse,1),mean(DSnorpResponse(:,idx==8),2),std(DSnorpResponse(:,idx==8),[],2),'k');
xlim([0 size(DSnorpResponse,1)]); ylim([0 1]);
subplot(2,10,19);
shadedErrorBar(1:size(DSnorpResponse,1),mean(DSnorpResponse(:,idx==9),2),std(DSnorpResponse(:,idx==9),[],2),'k');
xlim([0 size(DSnorpResponse,1)]); ylim([0 1]);
subplot(2,10,20);
shadedErrorBar(1:size(DSnorpResponse,1),mean(DSnorpResponse(:,idx==10),2),std(DSnorpResponse(:,idx==10),[],2),'k');
xlim([0 size(DSnorpResponse,1)]); ylim([0 1]);

%% plot gaussian mixture graphs
% figure
% y=idx==1;
% h = gscatter(data(:,1),data(:,2),y);
% hold on
% xaxislim=xlim;
% yaxislim=ylim;
% 
% scatter(GM.mu(1,1),GM.mu(1,2),'ko')
% scatter(GM.mu(2,1),GM.mu(2,2),'ko')
% 
% % normalize sigma
% % maxsigma=max(GM.Sigma(:));
% % absGM.Sigma(:,:,1)=GM.Sigma(:,:,1)./maxsigma;
% % absGM.Sigma(:,:,2)=GM.Sigma(:,:,2)./maxsigma;
% 
% gauss1=gmdistribution(GM.mu(1,:),GM.Sigma(:,:,1));
% gauss2=gmdistribution(GM.mu(2,:),GM.Sigma(:,:,2));
% ezcontour(@(x1,x2)pdf(gauss1,[x1 x2]),xlim,ylim);
% ezcontour(@(x1,x2)pdf(gauss2,[x1 x2]),xlim,ylim);
% legend(h,'Model 0','Model1')
% hold off
% 
% %ezcontour(@(x1,x2)pdf(GM,[x1 x2]),xlim,ylim,7);
% 
% 
% figure;
% scatter(data(idx==1,1),data(idx==1,2),10,post(idx==1,1),'+')
% hold on
% scatter(data(idx==2,1),data(idx==2,2),10,post(idx==2,1),'o')
% hold off
% legend('Cluster 1','Cluster 2','Location','NW')
% clrmap = jet(80); colormap(clrmap(9:72,:))
% ylabel(colorbar,'Cluster 1 posterior probability')
% 
% figure;
% [~,order] = sort(post(:,1));
% plot(1:size(data,1),post(order,1),'r-',1:size(data,1),post(order,2),'b-');
% legend({'Cluster 1' 'Cluster 2'},'location','NW');
% ylabel('Cluster membership score');
% xlabel('Point ranking');
% hold on
% plot([0 2446],[0.95 0.95],'--k');
% plot([0 2446],[0.05 0.05],'--k');



%% apply the posterior probability cutoff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots the posterior probability of cells belonging to each of the
% two first clusters
% uses postPcut as the cutoff for posterior probablity; cutoff of 0.5 means no restriction. 
figure;scatter(1:length(post(idx==1,1)),post(idx==1,1),'b')
figure;scatter(1:length(post(idx==2,2)),post(idx==2,2),'b')
figure;scatter(1:length(post(idx==3,3)),post(idx==3,3),'b')
figure;scatter(1:length(post(idx==4,4)),post(idx==4,4),'b')
try figure;scatter(1:length(post(idx==5,5)),post(idx==5,5),'b'); end
try figure;scatter(1:length(post(idx==6,6)),post(idx==6,6),'b'); end
try figure;scatter(1:length(post(idx==7,7)),post(idx==7,7),'b'); end
try figure;scatter(1:length(post(idx==8,8)),post(idx==8,8),'b'); end
try figure;scatter(1:length(post(idx==9,9)),post(idx==9,9),'b'); end
try figure;scatter(1:length(post(idx==10,10)),post(idx==10,10),'b'); end

indPost1=intersect(find(idx==1),find(post(:,1)<postPcut));
indPost2=intersect(find(idx==2),find(post(:,2)<postPcut));
indPost3=intersect(find(idx==3),find(post(:,3)<postPcut));
indPost4=intersect(find(idx==4),find(post(:,4)<postPcut));
try indPost5=intersect(find(idx==5),find(post(:,5)<postPcut)); end
try indPost6=intersect(find(idx==6),find(post(:,6)<postPcut)); end
try indPost7=intersect(find(idx==7),find(post(:,7)<postPcut)); end
try indPost8=intersect(find(idx==8),find(post(:,8)<postPcut)); end
try indPost9=intersect(find(idx==9),find(post(:,9)<postPcut)); end
try indPost10=intersect(find(idx==10),find(post(:,10)<postPcut)); end

combIndPost=[indPost1;indPost2;indPost3;indPost4;indPost5];
try combIndPost=[indPost1;indPost2;indPost3;indPost4;indPost5;indPost6]; end
try combIndPost=[indPost1;indPost2;indPost3;indPost4;indPost5;indPost6;indPost7]; end
try combIndPost=[indPost1;indPost2;indPost3;indPost4;indPost5;indPost6;indPost7;indPost8]; end
try combIndPost=[indPost1;indPost2;indPost3;indPost4;indPost5;indPost6;indPost7;indPost8;indPost9]; end
try combIndPost=[indPost1;indPost2;indPost3;indPost4;indPost5;indPost6;indPost7;indPost8;indPost9;indPost10]; end

% create dataT, inddataT, idxT, and DSnorpResponseT that include only cells with posterior probabilities higher than postPcut 
% inddarta is the 

dataT=data;
dataT(combIndPost,:)=[];    % remove all cells that do not exceed the posterior prob cutoff
inddata=(1:size(data,1))';
inddataT=inddata;
inddataT(combIndPost,:)=[];
idxT=idx;
idxT(combIndPost,:)=[];
DSnorpResponseT=DSnorpResponse;
DSnorpResponseT(:,combIndPost)=[];
postT=post;
postT(combIndPost,:)=[];

% figure;
% scatter(dataT(idxT==1,1),dataT(idxT==1,2),10,postT(idxT==1,1),'+')
% hold on
% scatter(dataT(idxT==2,1),dataT(idxT==2,2),10,postT(idxT==2,1),'o')
% hold off
% legend('Cluster 1','Cluster 2','Location','NW')
% clrmap = jet(80); colormap(clrmap(9:72,:))
% ylabel(colorbar,'Cluster 1 posterior probability')


% plots individual and mean responses of all DS cells clustered into 2
% groups.
% f12=figure;
% set(f12,'position',[5 150 950 800]);
% subplot(2,2,1);
% plot(DSnorpResponseT(:,idxT==1),'k');
% xlim([0 size(DSnorpResponseT,1)]);
% subplot(2,2,2);
% plot(DSnorpResponseT(:,idxT==2),'k');
% xlim([0 size(DSnorpResponseT,1)]);
% subplot(2,2,3);
% shadedErrorBar(1:size(DSnorpResponseT,1),mean(DSnorpResponseT(:,idxT==1),2),std(DSnorpResponseT(:,idxT==1),[],2),'k');
% xlim([0 size(DSnorpResponseT,1)]); ylim([0 1]);
% title(['n = ',num2str(sum(idxT==1))]);
% subplot(2,2,4);
% shadedErrorBar(1:size(DSnorpResponseT,1),mean(DSnorpResponseT(:,idxT==2),2),std(DSnorpResponseT(:,idxT==2),[],2),'k');
% xlim([0 size(DSnorpResponseT,1)]); ylim([0 1]);
% title(['n = ',num2str(sum(idxT==2))]);
% subtitle('all DS cells');


f12=figure('Color',[1 1 1],'Renderer','painters');
set(f12,'position',[5 150 1250 400]);

subplot(1,10,1,'Parent',f12,'YColor',[1 1 1],'XColor',[1 1 1]);
shadedErrorBar(time,mean(DSnorpResponseT(:,idxT==1),2),std(DSnorpResponseT(:,idxT==1),[],2),'k');
hold on
plot(time,mean(DSnorpResponseT(:,idxT==1),2),'b','LineWidth',1);
xlim([0 time(end)]); ylim([0 1]);
ylabel('\Delta\itF/F','FontSize',16);
xlabel('Time (sec)','FontSize',16);
title(['\itn\rm = ',num2str(sum(idxT==1))],'FontSize',16);
box off

subplot(1,10,2,'Parent',f12,'YColor',[1 1 1],'XColor',[1 1 1]);
shadedErrorBar(time,mean(DSnorpResponseT(:,idxT==2),2),std(DSnorpResponseT(:,idxT==2),[],2),'k');
hold on
plot(time,mean(DSnorpResponseT(:,idxT==2),2),'b','LineWidth',1);
xlim([0 time(end)]); ylim([0 1]);
xlabel('Time (sec)','FontSize',16);
title(['\itn\rm = ',num2str(sum(idxT==2))],'FontSize',16);
box off

subplot(1,10,3,'Parent',f12,'YColor',[1 1 1],'XColor',[1 1 1]);
shadedErrorBar(time,mean(DSnorpResponseT(:,idxT==3),2),std(DSnorpResponseT(:,idxT==3),[],2),'k');
hold on
plot(time,mean(DSnorpResponseT(:,idxT==3),2),'b','LineWidth',1);
xlim([0 time(end)]); ylim([0 1]);
xlabel('Time (sec)','FontSize',16);
title(['\itn\rm = ',num2str(sum(idxT==3))],'FontSize',16);
box off

subplot(1,10,4,'Parent',f12,'YColor',[1 1 1],'XColor',[1 1 1]);
shadedErrorBar(time,mean(DSnorpResponseT(:,idxT==4),2),std(DSnorpResponseT(:,idxT==4),[],2),'k');
hold on
plot(time,mean(DSnorpResponseT(:,idxT==4),2),'b','LineWidth',1);
xlim([0 time(end)]); ylim([0 1]);
xlabel('Time (sec)','FontSize',16);
title(['\itn\rm = ',num2str(sum(idxT==4))],'FontSize',16);
box off

subplot(1,10,5,'Parent',f12,'YColor',[1 1 1],'XColor',[1 1 1]);
shadedErrorBar(time,mean(DSnorpResponseT(:,idxT==5),2),std(DSnorpResponseT(:,idxT==5),[],2),'k');
hold on
plot(time,mean(DSnorpResponseT(:,idxT==5),2),'b','LineWidth',1);
xlim([0 time(end)]); ylim([0 1]);
xlabel('Time (sec)','FontSize',16);
title(['\itn\rm = ',num2str(sum(idxT==5))],'FontSize',16);
box off

subplot(1,10,6,'Parent',f12,'YColor',[1 1 1],'XColor',[1 1 1]);
shadedErrorBar(time,mean(DSnorpResponseT(:,idxT==6),2),std(DSnorpResponseT(:,idxT==6),[],2),'k');
hold on
plot(time,mean(DSnorpResponseT(:,idxT==6),2),'b','LineWidth',1);
xlim([0 time(end)]); ylim([0 1]);
xlabel('Time (sec)','FontSize',16);
title(['\itn\rm = ',num2str(sum(idxT==6))],'FontSize',16);
box off

subplot(1,10,7,'Parent',f12,'YColor',[1 1 1],'XColor',[1 1 1]);
shadedErrorBar(time,mean(DSnorpResponseT(:,idxT==7),2),std(DSnorpResponseT(:,idxT==7),[],2),'k');
hold on
plot(time,mean(DSnorpResponseT(:,idxT==7),2),'b','LineWidth',1);
xlim([0 time(end)]); ylim([0 1]);
xlabel('Time (sec)','FontSize',16);
title(['\itn\rm = ',num2str(sum(idxT==7))],'FontSize',16);
box off

subplot(1,10,8,'Parent',f12,'YColor',[1 1 1],'XColor',[1 1 1]);
shadedErrorBar(time,mean(DSnorpResponseT(:,idxT==8),2),std(DSnorpResponseT(:,idxT==8),[],2),'k');
hold on
plot(time,mean(DSnorpResponseT(:,idxT==8),2),'b','LineWidth',1);
xlim([0 time(end)]); ylim([0 1]);
xlabel('Time (sec)','FontSize',16);
title(['\itn\rm = ',num2str(sum(idxT==8))],'FontSize',16);
box off

subplot(1,10,9,'Parent',f12,'YColor',[1 1 1],'XColor',[1 1 1]);
shadedErrorBar(time,mean(DSnorpResponseT(:,idxT==9),2),std(DSnorpResponseT(:,idxT==9),[],2),'k');
hold on
plot(time,mean(DSnorpResponseT(:,idxT==9),2),'b','LineWidth',1);
xlim([0 time(end)]); ylim([0 1]);
xlabel('Time (sec)','FontSize',16);
title(['\itn\rm = ',num2str(sum(idxT==9))],'FontSize',16);
box off

subplot(1,10,10,'Parent',f12,'YColor',[1 1 1],'XColor',[1 1 1]);
shadedErrorBar(time,mean(DSnorpResponseT(:,idxT==10),2),std(DSnorpResponseT(:,idxT==10),[],2),'k');
hold on
plot(time,mean(DSnorpResponseT(:,idxT==10),2),'b','LineWidth',1);
xlim([0 time(end)]); ylim([0 1]);
xlabel('Time (sec)','FontSize',16);
title(['\itn\rm = ',num2str(sum(idxT==10))],'FontSize',16);
box off

%% cluster new data sets - SNS cells (10 ONDS and 10 ONOFFDS)
% % extract responses from the file: summary_All SNS cells
% filedir=cd;
% filedir2=[filedir,'\All SNS cells'];
% namesSNS=getFileList(filedir2,'summary_All SNS cells.mat',0,'anywhere');
% load(num2str(namesSNS{1}));
% 
% SNSpResponseONDS=summary.pResponseONDS(1:min(elements),:);
% SNSpResponseONOFFDS=summary.pResponseONOFFDS(1:min(elements),:);
% % standardize the length of responses
% SNSpResponseONDS(min(elements):size(SNSpResponseONDS,1),:)=[];
% SNSpResponseONOFFDS(min(elements):size(SNSpResponseONOFFDS,1),:)=[];
% % normalize the responses
% for p=1:size(SNSpResponseONDS,2)
%     SNSpResponseONDSnor(:,p)=mat2gray(SNSpResponseONDS(:,p));
% end
% for p=1:size(SNSpResponseONOFFDS,2)
%     SNSpResponseONOFFDSnor(:,p)=mat2gray(SNSpResponseONOFFDS(:,p));
% end
% % smooth the responses
% for p=1:size(SNSpResponseONDS,2)
%     SNSpResponseONDSs(:,p)=smooth(SNSpResponseONDS(:,p),15,'moving');
% end
% for p=1:size(SNSpResponseONOFFDS,2)
%     SNSpResponseONOFFDSs(:,p)=smooth(SNSpResponseONOFFDS(:,p),15,'moving');
% end
% % normalize the smoothed responses
% SNSpResponseONDSnorS=(SNSpResponseONDSs-repmat(mean(SNSpResponseONDSs(10:40,:),1),size(SNSpResponseONDSs,1),1))./(repmat(max(SNSpResponseONDSs,[],1),size(SNSpResponseONDSs,1),1)-repmat(mean(SNSpResponseONDSs(10:40,:),1),size(SNSpResponseONDSs,1),1));
% SNSpResponseONOFFDSnorS=(SNSpResponseONOFFDSs-repmat(mean(SNSpResponseONOFFDSs(10:40,:),1),size(SNSpResponseONOFFDSs,1),1))./(repmat(max(SNSpResponseONOFFDSs,[],1),size(SNSpResponseONOFFDSs,1),1)-repmat(mean(SNSpResponseONOFFDSs(10:40,:),1),size(SNSpResponseONOFFDSs,1),1));
% SNSpResponseS=[SNSpResponseONDSs,SNSpResponseONOFFDSs];
% SNSpResponsenorS=[SNSpResponseONDSnorS,SNSpResponseONOFFDSnorS];
% SNSonset1=[summary.onset1ONDS,summary.onset1ONOFFDS];
% 
% % switch between using the normalized responses or not
% switch NORMALRESstate
%     case 1
%         SNSresponse=[SNSpResponseONDSnorS,SNSpResponseONOFFDSnorS];
%     case 0
%         SNSresponse=[SNSpResponseONDSs,SNSpResponseONOFFDSs];
% end
% 
% 
% %% calculate parameters for clustering SNS cells
% figure;
% SNSpeakMinMag=SNSresponse(130,:);     %data point number 130 is where the minima occur for both ON and ONOFF cells
% for p=1:size(SNSresponse,2)
%     peakfinder(SNSresponse(:,p),0.1,0.1,1,false); drawnow expose
%     [peakLocMax,peakMagMax]=peakfinder(SNSresponse(:,p),0.1,0.1,1,false);
%     
%     SNSpeakLocMax{p}=peakLocMax;
%     SNSpeakMagMax{p}=peakMagMax;
%     
%     %%%% excludes from analysis every cell with a weak and noisy response (more than 4 detected peaks)
%     if size(SNSpeakLocMax{p},1)>4
%         SNSresponse(:,p)=0;
%     end
%     
%     %%% find OFF peak
%     offind=intersect(find(peakLocMax>145),find(peakLocMax<170));
%     if isempty(offind)
%         SNSpeakLocOff(p)=NaN;
%     else
%         SNSpeakLocOff(p)=peakLocMax(offind(1),1);
%     end
%     
%     %%% calculate slope of OFF peak
%     if ~isnan(SNSpeakLocOff(p))
%         yOffData=SNSresponse(SNSpeakLocOff(p):SNSpeakLocOff(p)+15,p);
%         % Set up fittype and options.
%         ft = fittype('poly1');
%         opts = fitoptions(ft);
%         % Fit model to data.
%         [fitresultOff, gofOff] = fit((1:size(yOffData,1))', yOffData(:,1), ft, opts );
%         SNSslopeOff(p)=fitresultOff.p1;
%     else
%         SNSslopeOff(p)=NaN;
%     end
%     
%     
%     %%% find ON peak
%     onind=find(peakLocMax>60);
%     if isempty(onind)
%         SNSpeakLocOn(p)=NaN;
%     else
%         SNSpeakLocOn(p)=peakLocMax(onind(1),1);
%     end
%     
%     %%% calculate slope of ON peak
%     if ~isnan(SNSpeakLocOn(p))
%         yOnData=SNSresponse(SNSpeakLocOn(p):SNSpeakLocOn(p)+slopePoints,p);
%         % Fit model to data.
%         [fitresultOn, gofOn] = fit((1:size(yOnData,1))', yOnData(:,1), ft, opts );
%         SNSslopeOn(p)=fitresultOn.p1;
%     else
%         SNSslopeOn(p)=NaN;
%     end
% end
% 
% SNSonSumResponse=sum(SNSresponse(60:130,:),1);
% SNSoffSumResponse=sum(SNSresponse(145:210,:),1);
% 
% SNSpeakLocOn=time(SNSpeakLocOn);      % convert DSpeakLocOn from number of data points to time
% SNSpeakLocOn=SNSpeakLocOn-SNSonset1';    % adjust DSpeakLocOn relative to the onset of the stimulus
% 
% SNSpeakLocOn=SNSpeakLocOn;
% SNSpeakLocOff=SNSpeakLocOff';
% SNSonSumResponse=SNSonSumResponse';
% SNSoffSumResponse=SNSoffSumResponse';
% SNSpeakMinMag=SNSpeakMinMag';
% SNSslopeOn=SNSslopeOn';
% SNSslopeOff=SNSslopeOff';
% 
% %data=[SNSpeakLocOn,SNSpeakLocOff,SNSonSumResponse,SNSoffSumResponse,SNSpeakMinMag,SNSslopeOn,SNSslopeOff];
% SNSdata=[SNSpeakLocOn,SNSslopeOn];
% 
% 
% 
% SNSidx = cluster(GM,SNSdata);
% SNSpost=posterior(GM,SNSdata);        %calculates the posterior probability of cells
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;scatter(1:length(SNSpost(SNSidx==1,1)),SNSpost(SNSidx==1,1),'b')
% figure;scatter(1:length(SNSpost(SNSidx==2,2)),SNSpost(SNSidx==2,2),'r')
% SNSindPost1=intersect(find(SNSidx==1),find(SNSpost(:,1)<postPcut));
% SNSindPost2=intersect(find(SNSidx==2),find(SNSpost(:,2)<postPcut));
% SNScombIndPost=[SNSindPost1;SNSindPost2];
% % create dataT, inddataT, SNSidxT, and DSnorpResponseT that include only cells with posterior probabilities higher than postPcut 
% % inddarta is the index to be used to recover the original indexing 
% SNSdataT=SNSdata;
% SNSdataT(SNScombIndPost,:)=[];
% SNSinddata=(1:size(SNSdata,1))';
% SNSinddataT=SNSinddata;
% SNSinddataT(SNScombIndPost,:)=[];
% SNSidxT=SNSidx;
% SNSidxT(SNScombIndPost,:)=[];
% SNSpResponsenorST=SNSpResponsenorS;
% SNSpResponsenorST(:,SNScombIndPost)=[];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % cluster SNS cells
% ONDSn=10;
% ONOFFDSn=10;
% SNSONclust1=0;
% SNSONclust2=0;
% SNSONOFFclust1=0;
% SNSONOFFclust2=0;
% for s=1:size(SNSinddataT,1)
%     if SNSinddataT(s)<=ONDSn
%         if SNSidxT(s)==1
%             SNSONclust1=SNSONclust1+1;
%         elseif SNSidxT(s)==2
%             SNSONclust2=SNSONclust2+1;
%         end
%     end
%     if SNSinddataT(s)>ONDSn
%         if SNSidxT(s)==1
%             SNSONOFFclust1=SNSONOFFclust1+1;
%         elseif SNSidxT(s)==2
%             SNSONOFFclust2=SNSONOFFclust2+1;
%         end
%     end
% end
% 
% % if SNS2clust2<SNS2clust1
% %     SNSmisclass=100*((SNS1clust1+SNS2clust2)/(SNS1clust1+SNS2clust2+SNS1clust2+SNS2clust1));
% % else
% %     SNSmisclass=100*((SNS1clust2+SNS2clust1)/(SNS1clust1+SNS2clust2+SNS1clust2+SNS2clust1));
% % end
% 
% % plots the clustered SNS cells responses
% % f13=figure;
% % set(f13,'position',[5 150 950 800]);
% % subplot(2,2,1);
% % plot(SNSpResponsenorST(:,SNSidxT==1),'k');
% % xlim([0 size(SNSpResponsenorST,1)]); ylim([0 1]);
% % subplot(2,2,2);
% % plot(SNSpResponsenorST(:,SNSidxT==2),'k');
% % xlim([0 size(SNSpResponsenorST,1)]); ylim([0 1]);
% subplot(2,2,3);
% shadedErrorBar(1:size(SNSpResponsenorST,1),mean(SNSpResponsenorST(:,SNSidxT==1),2),std(SNSpResponsenorST(:,SNSidxT==1),[],2),'k');
% xlim([0 size(SNSpResponsenorST,1)]); ylim([0 1]);
% title(['n = ',num2str(sum(SNSidxT==1))]);
% subplot(2,2,4);
% shadedErrorBar(1:size(SNSpResponsenorST,1),mean(SNSpResponsenorST(:,SNSidxT==2),2),std(SNSpResponsenorST(:,SNSidxT==2),[],2),'k');
% xlim([0 size(SNSpResponsenorST,1)]); ylim([0 1]);
% title(['n = ',num2str(sum(SNSidxT==2))]);
% % subtitle('SNS cells');
% 
% f13=figure('Color',[1 1 1],'Renderer','painters');
% set(f13,'position',[5 150 950 400]);
% subplot(1,2,1,'Parent',f13,'YColor',[1 1 1],'XColor',[1 1 1]);
% shadedErrorBar(time(1:end-1,1),mean(SNSpResponsenorST(:,SNSidxT==1),2),std(SNSpResponsenorST(:,SNSidxT==1),[],2),'k');
% hold on
% plot(time(1:end-1,1),mean(SNSpResponsenorST(:,SNSidxT==1),2),'r','LineWidth',1);
% xlim([0 time(end-1)]); ylim([0 1]);
% ylabel('\Delta\itF/F','FontSize',16);
% xlabel('Time (sec)','FontSize',16);
% title(['\itn\rm = ',num2str(sum(SNSidxT==1))],'FontSize',16);
% box off
% subplot(1,2,2,'Parent',f13,'YColor',[1 1 1],'XColor',[1 1 1]);
% shadedErrorBar(time(1:end-1,1),mean(SNSpResponsenorST(:,SNSidxT==2),2),std(SNSpResponsenorST(:,SNSidxT==2),[],2),'k');
% hold on
% plot(time(1:end-1,1),mean(SNSpResponsenorST(:,SNSidxT==2),2),'b','LineWidth',1);
% xlim([0 time(end-1)]); ylim([0 1]);
% xlabel('Time (sec)','FontSize',16);
% title(['\itn\rm = ',num2str(sum(SNSidxT==2))],'FontSize',16);
% box off


%% cluster new data sets - CART ONandONOFFDS cells
% % extract responses from the following file:
% % pooledMapAlphaCorr_Calcium imaging data_40_DSI0.30_1.00_created_May_15_2015_14_47_CARTONOFFDS
% filedir=cd;
% namesCART=getFileList(filedir,'pooledMapAlphaCorr_Calcium imaging data_40_DSI0.30_1.00_created_Oct_08_2015_16_42_ONandONOFF_CART',0,'anywhere');
% CART=load(num2str(namesCART{1}));
% % standardize the length of CART responses
% % determine the minimum number of elements
% 
% % standardize the length of responses
% for nu=1:numel(CART.pooledmap.pResponse)
%     CARTnumPresponseT{nu}=CART.pooledmap.pResponse{nu}(1:min(elements));
% end
% CARTpResponse=reshape(cell2mat(CARTnumPresponseT'),min(elements),[]);
% % normalize the responses
% for p=1:size(CARTpResponse,2)
%     CARTnorpResponse(:,p)=mat2gray(CARTpResponse(:,p));
% end
% % smooth the reponses
% for p=1:size(CARTpResponse,2)
%     CARTpResponseS(:,p)=smooth(CARTpResponse(:,p),15,'moving');
% end
% % normalize the smoothed responses
% CARTpResponseNorS=(CARTpResponseS-repmat(mean(CARTpResponseS(10:40,:),1),size(CARTpResponseS,1),1))./(repmat(max(CARTpResponseS,[],1),size(CARTpResponseS,1),1)-repmat(mean(CARTpResponseS(10:40,:),1),size(CARTpResponseS,1),1));
% 
% % switch between using the normalized responses or not
% switch NORMALRESstate
%     case 1
%         CARTresponse=CARTpResponseNorS;
%     case 0
%         CARTresponse=CARTpResponseS;
% end
% 
% 
% CARTonset1=CART.pooledmap.onset1;
% 
% %% calculate parameters for clustering CART ONOFF cells
% figure;
% CARTpeakMinMag=CARTresponse(130,:);     %data point number 130 is where the minima occur for both ON and ONOFF cells
% for p=1:size(CARTresponse,2)
%     peakfinder(CARTresponse(:,p),0.1,0.1,1,false); drawnow expose
%     [peakLocMax,peakMagMax]=peakfinder(CARTresponse(:,p),0.1,0.1,1,false);
%     
%     CARTpeakLocMax{p}=peakLocMax;
%     CARTpeakMagMax{p}=peakMagMax;
%     
%     %%%% excludes from analysis every cell with a weak and noisy response (more than 4 detected peaks)
%     if size(CARTpeakLocMax{p},1)>4
%         CARTresponse(:,p)=0;
%     end
%     
%     %%% find OFF peak
%     offind=intersect(find(peakLocMax>145),find(peakLocMax<170));
%     if isempty(offind)
%         CARTpeakLocOff(p)=NaN;
%     else
%         CARTpeakLocOff(p)=peakLocMax(offind(1),1);
%     end
%     
%     %%% calculate slope of OFF peak
%     if ~isnan(CARTpeakLocOff(p))
%         yOffData=CARTresponse(CARTpeakLocOff(p):CARTpeakLocOff(p)+15,p);
%         % Set up fittype and options.
%         ft = fittype('poly1');
%         opts = fitoptions(ft);
%         % Fit model to data.
%         [fitresultOff, gofOff] = fit((1:size(yOffData,1))', yOffData(:,1), ft, opts );
%         CARTslopeOff(p)=fitresultOff.p1;
%     else
%         CARTslopeOff(p)=NaN;
%     end
%     
%     
%     %%% find ON peak
%     onind=find(peakLocMax>60);
%     if isempty(onind)
%         CARTpeakLocOn(p)=NaN;
%     else
%         CARTpeakLocOn(p)=peakLocMax(onind(1),1);
%     end
%     
%     %%% calculate slope of ON peak
%     if ~isnan(CARTpeakLocOn(p))
%         yOnData=CARTresponse(CARTpeakLocOn(p):CARTpeakLocOn(p)+slopePoints,p);
%         % Fit model to data.
%         [fitresultOn, gofOn] = fit((1:size(yOnData,1))', yOnData(:,1), ft, opts );
%         CARTslopeOn(p)=fitresultOn.p1;
%     else
%         CARTslopeOn(p)=NaN;
%     end
% end
% 
% CARTonSumResponse=sum(CARTresponse(60:130,:),1);
% CARToffSumResponse=sum(CARTresponse(145:210,:),1);
% 
% CARTpeakLocOn=time(CARTpeakLocOn);      % convert DSpeakLocOn from number of data points to time
% CARTpeakLocOn=CARTpeakLocOn-CARTonset1;    % adjust DSpeakLocOn relative to the onset of the stimulus
% 
% 
% CARTpeakLocOn=CARTpeakLocOn;
% CARTpeakLocOff=CARTpeakLocOff';
% CARTonSumResponse=CARTonSumResponse';
% CARToffSumResponse=CARToffSumResponse';
% CARTpeakMinMag=CARTpeakMinMag';
% CARTslopeOn=CARTslopeOn';
% CARTslopeOff=CARTslopeOff';
% 
% %data=[CARTpeakLocOn,CARTpeakLocOff,CARTonSumResponse,CARToffSumResponse,CARTpeakMinMag,CARTslopeOn,CARTslopeOff];
% CARTdata=[CARTpeakLocOn,CARTslopeOn];
% 
% 
% 
% CARTidx = cluster(GM,CARTdata);
% CARTpost=posterior(GM,CARTdata);        %calculates the posterior probability of cells
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;scatter(1:length(CARTpost(CARTidx==1,1)),CARTpost(CARTidx==1,1),'b')
% figure;scatter(1:length(CARTpost(CARTidx==2,2)),CARTpost(CARTidx==2,2),'r')
% CARTindPost1=intersect(find(CARTidx==1),find(CARTpost(:,1)<postPcut));
% CARTindPost2=intersect(find(CARTidx==2),find(CARTpost(:,2)<postPcut));
% CARTcombIndPost=[CARTindPost1;CARTindPost2];
% % create CARTdataT, CARTinddataT, CARTidxT, and CARTnorpResponseT that include only cells with posterior probabilities higher than postPcut 
% % inddarta is the index to be used to recover the original indexing 
% CARTdataT=CARTdata;
% CARTdataT(CARTcombIndPost,:)=[];
% CARTinddata=(1:size(CARTdata,1))';
% CARTinddataT=CARTinddata;
% CARTinddataT(CARTcombIndPost,:)=[];
% CARTidxT=CARTidx;
% CARTidxT(CARTcombIndPost,:)=[];
% CARTnorpResponseT=CARTnorpResponse;
% CARTnorpResponseT(:,CARTcombIndPost)=[];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % cluster CART cells
% CARTclust1=numel(find(CARTidxT==1));
% CARTclust2=numel(find(CARTidxT==2));
% if CARTclust2<CARTclust1
%     CARTmisclass=100*(CARTclust2/(CARTclust1+CARTclust2));
% else
%     CARTmisclass=100*(CARTclust1/(CARTclust1+CARTclust2));
% end
%     
% % plots individual and mean responses of all CART positive DS cells clustered into 2
% % groups.
% % f14=figure;
% % set(f14,'position',[5 150 950 800]);
% % subplot(2,2,1);
% % plot(CARTnorpResponseT(:,CARTidxT==1),'k');
% % xlim([0 size(CARTnorpResponseT,1)]);
% % subplot(2,2,2);
% % plot(CARTnorpResponseT(:,CARTidxT==2),'k');
% % xlim([0 size(CARTnorpResponseT,1)]);
% % subplot(2,2,3);
% % shadedErrorBar(1:size(CARTnorpResponseT,1),mean(CARTnorpResponseT(:,CARTidxT==1),2),std(CARTnorpResponseT(:,CARTidxT==1),[],2),'k');
% % xlim([0 size(CARTnorpResponseT,1)]); ylim([0 1]);
% % title(['n = ',num2str(sum(CARTidxT==1))]);
% % subplot(2,2,4);
% % shadedErrorBar(1:size(CARTnorpResponseT,1),mean(CARTnorpResponseT(:,CARTidxT==2),2),std(CARTnorpResponseT(:,CARTidxT==2),[],2),'k');
% % xlim([0 size(CARTnorpResponseT,1)]); ylim([0 1]);
% % title(['n = ',num2str(sum(CARTidxT==2))]);
% % subtitle('CART positive');
% 
% 
% f14=figure('Color',[1 1 1],'Renderer','painters');
% set(f14,'position',[5 150 950 400]);
% subplot(1,2,1,'Parent',f14,'YColor',[1 1 1],'XColor',[1 1 1]);
% shadedErrorBar(time,mean(CARTnorpResponseT(:,CARTidxT==1),2),std(CARTnorpResponseT(:,CARTidxT==1),[],2),'k');
% hold on
% plot(time,mean(CARTnorpResponseT(:,CARTidxT==1),2),'r','LineWidth',1);
% xlim([0 time(end)]); ylim([0 1]);
% ylabel('\Delta\itF/F','FontSize',16);
% xlabel('Time (sec)','FontSize',16);
% title(['\itn\rm = ',num2str(sum(CARTidxT==1))],'FontSize',16);
% box off
% subplot(1,2,2,'Parent',f14,'YColor',[1 1 1],'XColor',[1 1 1]);
% shadedErrorBar(time,mean(CARTnorpResponseT(:,CARTidxT==2),2),std(CARTnorpResponseT(:,CARTidxT==2),[],2),'k');
% hold on
% plot(time,mean(CARTnorpResponseT(:,CARTidxT==2),2),'b','LineWidth',1);
% xlim([0 time(end)]); ylim([0 1]);
% xlabel('Time (sec)','FontSize',16);
% title(['\itn\rm = ',num2str(sum(CARTidxT==2))],'FontSize',16);
% box off
% 
% 
%% cluster new data sets - IF retro DS cells (both ON and ONOFF were included)
% % extract responses from the following file:
% % pooledMapAlphaCorr_Calcium imaging data_40_DSI0.30_1.00_created_May_15_2015_15_30_IFDS
% filedir=cd;
% %namesIF=getFileList(filedir,'pooledMapAlphaCorr_Calcium imaging data_40_DSI0.30_1.00_created_Jul_16_2015_14_07_feature-based clustering_IF',0,'anywhere');
% namesIF=getFileList(filedir,'pooledMapAlphaCorr_Calcium imaging data_40_DSI0.30_1.00_created_Oct_08_2015_16_40_ONandONOFF_IF',0,'anywhere');
% IF=load(num2str(namesIF{1}));
% % standardize the length of IF responses
% % determine the minimum number of elements
% 
% % standardize the length of responses
% for nu=1:numel(IF.pooledmap.pResponse)
%     IFnumPresponseT{nu}=IF.pooledmap.pResponse{nu}(1:min(elements));
% end
% IFpResponse=reshape(cell2mat(IFnumPresponseT'),min(elements),[]);
% % normalize the responses
% for p=1:size(IFpResponse,2)
%     IFnorpResponse(:,p)=mat2gray(IFpResponse(:,p));
% end
% % smooth the reponses
% for p=1:size(IFpResponse,2)
%     IFpResponseS(:,p)=smooth(IFpResponse(:,p),15,'moving');
% end
% % normalize the smoothed responses
% IFpResponseNorS=(IFpResponseS-repmat(mean(IFpResponseS(10:40,:),1),size(IFpResponseS,1),1))./(repmat(max(IFpResponseS,[],1),size(IFpResponseS,1),1)-repmat(mean(IFpResponseS(10:40,:),1),size(IFpResponseS,1),1));
% 
% % switch between using the normalized responses or not
% switch NORMALRESstate
%     case 1
%         IFresponse=IFpResponseNorS;
%     case 0
%         IFresponse=IFpResponseS;
% end
% 
% 
% IFonset1=IF.pooledmap.onset1;
% 
% %% calculate parameters for clustering MTN and SF cells
% figure;
% IFpeakMinMag=IFresponse(130,:);     %data point number 130 is where the minima occur for both ON and ONOFF cells
% for p=1:size(IFresponse,2)
%     peakfinder(IFresponse(:,p),0.1,0.1,1,false); drawnow expose
%     [peakLocMax,peakMagMax]=peakfinder(IFresponse(:,p),0.1,0.1,1,false);
%     
%     IFpeakLocMax{p}=peakLocMax;
%     IFpeakMagMax{p}=peakMagMax;
%     
%     %%%% excludes from analysis every cell with a weak and noisy response (more than 4 detected peaks)
%     if size(IFpeakLocMax{p},1)>4
%         IFresponse(:,p)=0;
%     end
%     
%     %%% find OFF peak
%     offind=intersect(find(peakLocMax>145),find(peakLocMax<170));
%     if isempty(offind)
%         IFpeakLocOff(p)=NaN;
%     else
%         IFpeakLocOff(p)=peakLocMax(offind(1),1);
%     end
%     
%     %%% calculate slope of OFF peak
%     if ~isnan(IFpeakLocOff(p))
%         yOffData=IFresponse(IFpeakLocOff(p):IFpeakLocOff(p)+15,p);
%         % Set up fittype and options.
%         ft = fittype('poly1');
%         opts = fitoptions(ft);
%         % Fit model to data.
%         [fitresultOff, gofOff] = fit((1:size(yOffData,1))', yOffData(:,1), ft, opts );
%         IFslopeOff(p)=fitresultOff.p1;
%     else
%         IFslopeOff(p)=NaN;
%     end
%     
%     
%     %%% find ON peak
%     onind=find(peakLocMax>60);
%     if isempty(onind)
%         IFpeakLocOn(p)=NaN;
%     else
%         IFpeakLocOn(p)=peakLocMax(onind(1),1);
%     end
%     
%     %%% calculate slope of ON peak
%     if ~isnan(IFpeakLocOn(p))
%         yOnData=IFresponse(IFpeakLocOn(p):IFpeakLocOn(p)+slopePoints,p);
%         % Fit model to data.
%         [fitresultOn, gofOn] = fit((1:size(yOnData,1))', yOnData(:,1), ft, opts );
%         IFslopeOn(p)=fitresultOn.p1;
%     else
%         IFslopeOn(p)=NaN;
%     end
% end
% 
% IFonSumResponse=sum(IFresponse(60:130,:),1);
% IFoffSumResponse=sum(IFresponse(145:210,:),1);
% 
% IFpeakLocOn=time(IFpeakLocOn);      % convert DSpeakLocOn from number of data points to time
% IFpeakLocOn=IFpeakLocOn-IFonset1;    % adjust DSpeakLocOn relative to the onset of the stimulus
% 
% IFpeakLocOn=IFpeakLocOn;
% IFpeakLocOff=IFpeakLocOff';
% IFonSumResponse=IFonSumResponse';
% IFoffSumResponse=IFoffSumResponse';
% IFpeakMinMag=IFpeakMinMag';
% IFslopeOn=IFslopeOn';
% IFslopeOff=IFslopeOff';
% 
% %data=[IFpeakLocOn,IFpeakLocOff,IFonSumResponse,IFoffSumResponse,IFpeakMinMag,IFslopeOn,IFslopeOff];
% IFdata=[IFpeakLocOn,IFslopeOn];
% 
% 
% IFidx = cluster(GM,IFdata);
% IFpost=posterior(GM,IFdata);        %calculates the posterior probability of cells
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;scatter(1:length(IFpost(IFidx==1,1)),IFpost(IFidx==1,1),'b')
% figure;scatter(1:length(IFpost(IFidx==2,2)),IFpost(IFidx==2,2),'r')
% IFindPost1=intersect(find(IFidx==1),find(IFpost(:,1)<postPcut));
% IFindPost2=intersect(find(IFidx==2),find(IFpost(:,2)<postPcut));
% IFcombIndPost=[IFindPost1;IFindPost2];
% % create IFdataT, IFinddataT, IFidxT, and IFnorpResponseT that include only cells with posterior probabilities higher than postPcut 
% % inddarta is the index to be used to recover the original indexing 
% IFdataT=IFdata;
% IFdataT(IFcombIndPost,:)=[];
% IFinddata=(1:size(IFdata,1))';
% IFinddataT=IFinddata;
% IFinddataT(IFcombIndPost,:)=[];
% IFidxT=IFidx;
% IFidxT(IFcombIndPost,:)=[];
% IFnorpResponseT=IFnorpResponse;
% IFnorpResponseT(:,IFcombIndPost)=[];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % cluster IF cells
% IFclust1=numel(find(IFidx==1));
% IFclust2=numel(find(IFidx==2));
% if IFclust2>IFclust1
%     IFmisclass=100*(IFclust2/(IFclust1+IFclust2));
% else
%     IFmisclass=100*(IFclust1/(IFclust1+IFclust2));
% end
% 
% % plots individual and mean responses of all IF retro cells DS cells clustered into 2
% % groups.
% % f15=figure;
% % set(f15,'position',[5 150 950 800]);
% % subplot(2,2,1);
% % plot(IFnorpResponseT(:,IFidxT==1),'k');
% % xlim([0 size(IFnorpResponseT,1)]);
% % subplot(2,2,2);
% % plot(IFnorpResponseT(:,IFidxT==2),'k');
% % xlim([0 size(IFnorpResponseT,1)]);
% % subplot(2,2,3);
% % shadedErrorBar(1:size(IFnorpResponseT,1),mean(IFnorpResponseT(:,IFidxT==1),2),std(IFnorpResponseT(:,IFidxT==1),[],2),'k');
% % xlim([0 size(IFnorpResponseT,1)]); ylim([0 1]);
% % title(['n = ',num2str(sum(IFidxT==1))]);
% % subplot(2,2,4);
% % shadedErrorBar(1:size(IFnorpResponseT,1),mean(IFnorpResponseT(:,IFidxT==2),2),std(IFnorpResponseT(:,IFidxT==2),[],2),'k');
% % xlim([0 size(IFnorpResponseT,1)]); ylim([0 1]);
% % title(['n = ',num2str(sum(IFidxT==2))]);
% % subtitle('IF retro cells');
% 
% f15=figure('Color',[1 1 1],'Renderer','painters');
% set(f15,'position',[5 150 950 400]);
% subplot(1,2,1,'Parent',f15,'YColor',[1 1 1],'XColor',[1 1 1]);
% shadedErrorBar(time,mean(IFnorpResponseT(:,IFidxT==1),2),std(IFnorpResponseT(:,IFidxT==1),[],2),'k');
% hold on
% plot(time,mean(IFnorpResponseT(:,IFidxT==1),2),'r','LineWidth',1);
% xlim([0 time(end)]); ylim([0 1]);
% ylabel('\Delta\itF/F','FontSize',16);
% xlabel('Time (sec)','FontSize',16);
% title(['\itn\rm = ',num2str(sum(IFidxT==1))],'FontSize',16);
% box off
% subplot(1,2,2,'Parent',f15,'YColor',[1 1 1],'XColor',[1 1 1]);
% shadedErrorBar(time,mean(IFnorpResponseT(:,IFidxT==2),2),std(IFnorpResponseT(:,IFidxT==2),[],2),'k');
% hold on
% plot(time,mean(IFnorpResponseT(:,IFidxT==2),2),'b','LineWidth',1);
% xlim([0 time(end)]); ylim([0 1]);
% xlabel('Time (sec)','FontSize',16);
% title(['\itn\rm = ',num2str(sum(IFidxT==2))],'FontSize',16);
% box off
% 
% 
%% cluster new data sets - SF retro DS cells (both ON and ONOFF were included)
% % extract responses from the following file:
% % pooledMapAlphaCorr_Calcium imaging data_40_DSI0.30_1.00_created_May_15_2015_15_45_SFDS
% filedir=cd;
% %namesSF=getFileList(filedir,'pooledMapAlphaCorr_Calcium imaging data_40_DSI0.30_1.00_created_Jul_16_2015_14_08_feature-based clustering_SF',0,'anywhere');
% namesSF=getFileList(filedir,'pooledMapAlphaCorr_Calcium imaging data_40_DSI0.30_1.00_created_Oct_08_2015_16_41_ONandONOFF_SF',0,'anywhere');
% SF=load(num2str(namesSF{1}));
% % standardize the length of SF responses
% % determine the minimum number of elements
% 
% % standardize the length of responses
% for nu=1:numel(SF.pooledmap.pResponse)
%     SFnumPresponseT{nu}=SF.pooledmap.pResponse{nu}(1:min(elements));
% end
% SFpResponse=reshape(cell2mat(SFnumPresponseT'),min(elements),[]);
% % normalize the responses
% for p=1:size(SFpResponse,2)
%     SFnorpResponse(:,p)=mat2gray(SFpResponse(:,p));
% end
% % smooth the reponses
% for p=1:size(SFpResponse,2)
%     SFpResponseS(:,p)=smooth(SFpResponse(:,p),15,'moving');
% end
% % normalize the smoothed responses
% SFpResponseNorS=(SFpResponseS-repmat(mean(SFpResponseS(10:40,:),1),size(SFpResponseS,1),1))./(repmat(max(SFpResponseS,[],1),size(SFpResponseS,1),1)-repmat(mean(SFpResponseS(10:40,:),1),size(SFpResponseS,1),1));
% 
% % switch between using the normalized responses or not
% switch NORMALRESstate
%     case 1
%         SFresponse=SFpResponseNorS;
%     case 0
%         SFresponse=SFpResponseS;
% end
% 
% 
% SFonset1=SF.pooledmap.onset1;
% 
% %% calculate parameters for clustering SF cells
% figure;
% SFpeakMinMag=SFresponse(130,:);     %data point number 130 is where the minima occur for both ON and ONOFF cells
% for p=1:size(SFresponse,2)
%     peakfinder(SFresponse(:,p),0.1,0.1,1,false); drawnow expose
%     [peakLocMax,peakMagMax]=peakfinder(SFresponse(:,p),0.1,0.1,1,false);
%     
%     SFpeakLocMax{p}=peakLocMax;
%     SFpeakMagMax{p}=peakMagMax;
%     
%     %%%% excludes from analysis every cell with a weak and noisy response (more than 4 detected peaks)
%     if size(SFpeakLocMax{p},1)>4
%         SFresponse(:,p)=0;
%     end
%     
%     %%% find OFF peak
%     offind=intersect(find(peakLocMax>145),find(peakLocMax<170));
%     if isempty(offind)
%         SFpeakLocOff(p)=NaN;
%     else
%         SFpeakLocOff(p)=peakLocMax(offind(1),1);
%     end
%     
%     %%% calculate slope of OFF peak
%     if ~isnan(SFpeakLocOff(p))
%         yOffData=SFresponse(SFpeakLocOff(p):SFpeakLocOff(p)+15,p);
%         % Set up fittype and options.
%         ft = fittype('poly1');
%         opts = fitoptions(ft);
%         % Fit model to data.
%         [fitresultOff, gofOff] = fit((1:size(yOffData,1))', yOffData(:,1), ft, opts );
%         SFslopeOff(p)=fitresultOff.p1;
%     else
%         SFslopeOff(p)=NaN;
%     end
%     
%     
%     %%% find ON peak
%     onind=find(peakLocMax>60);
%     if isempty(onind)
%         SFpeakLocOn(p)=NaN;
%     else
%         SFpeakLocOn(p)=peakLocMax(onind(1),1);
%     end
%     
%     %%% calculate slope of ON peak
%     if ~isnan(SFpeakLocOn(p))
%         yOnData=SFresponse(SFpeakLocOn(p):SFpeakLocOn(p)+slopePoints,p);
%         % Fit model to data.
%         [fitresultOn, gofOn] = fit((1:size(yOnData,1))', yOnData(:,1), ft, opts );
%         SFslopeOn(p)=fitresultOn.p1;
%     else
%         SFslopeOn(p)=NaN;
%     end
% end
% 
% SFonSumResponse=sum(SFresponse(60:130,:),1);
% SFoffSumResponse=sum(SFresponse(145:210,:),1);
% 
% SFpeakLocOn=time(SFpeakLocOn);      % convert DSpeakLocOn from number of data points to time
% SFpeakLocOn=SFpeakLocOn-SFonset1;    % adjust DSpeakLocOn relative to the onset of the stimulus
% 
% SFpeakLocOn=SFpeakLocOn;
% SFpeakLocOff=SFpeakLocOff';
% SFonSumResponse=SFonSumResponse';
% SFoffSumResponse=SFoffSumResponse';
% SFpeakMinMag=SFpeakMinMag';
% SFslopeOn=SFslopeOn';
% SFslopeOff=SFslopeOff';
% 
% %data=[SFpeakLocOn,SFpeakLocOff,SFonSumResponse,SFoffSumResponse,SFpeakMinMag,SFslopeOn,SFslopeOff];
% SFdata=[SFpeakLocOn,SFslopeOn];
% 
% 
% 
% SFidx = cluster(GM,SFdata);
% SFpost=posterior(GM,SFdata);        %calculates the posterior probability of cells
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;scatter(1:length(SFpost(SFidx==1,1)),SFpost(SFidx==1,1),'b')
% figure;scatter(1:length(SFpost(SFidx==2,2)),SFpost(SFidx==2,2),'r')
% SFindPost1=intersect(find(SFidx==1),find(SFpost(:,1)<postPcut));
% SFindPost2=intersect(find(SFidx==2),find(SFpost(:,2)<postPcut));
% SFcombIndPost=[SFindPost1;SFindPost2];
% % create SFdataT, SFinddataT, SFidxT, and SFnorpResponseT that include only cells with posterior probabilities higher than postPcut 
% % inddarta is the index to be used to recover the original indexing 
% SFdataT=SFdata;
% SFdataT(SFcombIndPost,:)=[];
% SFinddata=(1:size(SFdata,1))';
% SFinddataT=SFinddata;
% SFinddataT(SFcombIndPost,:)=[];
% SFidxT=SFidx;
% SFidxT(SFcombIndPost,:)=[];
% SFnorpResponseT=SFnorpResponse;
% SFnorpResponseT(:,SFcombIndPost)=[];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % cluster SF cells
% SFclust1=numel(find(SFidx==1));
% SFclust2=numel(find(SFidx==2));
% if SFclust2<SFclust1
%     SFmisclass=100*(SFclust2/(SFclust1+SFclust2));
% else
%     SFmisclass=100*(SFclust1/(SFclust1+SFclust2));
% end
% 
% % plots individual and mean responses of all SF-OT retro cells DS cells clustered into 2
% % groups.
% % f16=figure;
% % set(f16,'position',[5 150 950 800]);
% % subplot(2,2,1);
% % plot(SFnorpResponseT(:,SFidxT==1),'k');
% % xlim([0 size(SFnorpResponseT,1)]);
% % subplot(2,2,2);
% % plot(SFnorpResponseT(:,SFidxT==2),'k');
% % xlim([0 size(SFnorpResponseT,1)]);
% % subplot(2,2,3);
% % shadedErrorBar(1:size(SFnorpResponseT,1),mean(SFnorpResponseT(:,SFidxT==1),2),std(SFnorpResponseT(:,SFidxT==1),[],2),'k');
% % xlim([0 size(SFnorpResponseT,1)]); ylim([0 1]);
% % title(['n = ',num2str(sum(SFidxT==1))]);
% % subplot(2,2,4);
% % shadedErrorBar(1:size(SFnorpResponseT,1),mean(SFnorpResponseT(:,SFidxT==2),2),std(SFnorpResponseT(:,SFidxT==2),[],2),'k');
% % xlim([0 size(SFnorpResponseT,1)]); ylim([0 1]);
% % title(['n = ',num2str(sum(SFidxT==2))]);
% % subtitle('SF retro cells');
% 
% 
% 
% f16=figure('Color',[1 1 1],'Renderer','painters');
% set(f16,'position',[5 150 950 400]);
% subplot(1,2,1,'Parent',f16,'YColor',[1 1 1],'XColor',[1 1 1]);
% shadedErrorBar(time,mean(SFnorpResponseT(:,SFidxT==1),2),std(SFnorpResponseT(:,SFidxT==1),[],2),'k');
% hold on
% plot(time,mean(SFnorpResponseT(:,SFidxT==1),2),'r','LineWidth',1);
% xlim([0 time(end)]); ylim([0 1]);
% ylabel('\Delta\itF/F','FontSize',16);
% xlabel('Time (sec)','FontSize',16);
% title(['\itn\rm = ',num2str(sum(SFidxT==1))],'FontSize',16);
% box off
% subplot(1,2,2,'Parent',f16,'YColor',[1 1 1],'XColor',[1 1 1]);
% shadedErrorBar(time,mean(SFnorpResponseT(:,SFidxT==2),2),std(SFnorpResponseT(:,SFidxT==2),[],2),'k');
% hold on
% plot(time,mean(SFnorpResponseT(:,SFidxT==2),2),'b','LineWidth',1);
% xlim([0 time(end)]); ylim([0 1]);
% xlabel('Time (sec)','FontSize',16);
% title(['\itn\rm = ',num2str(sum(SFidxT==2))],'FontSize',16);
% box off
% 
%% cluster new data sets - MTN retro DS cells (both ON and ONOFF were included)
% % extract responses from the following file:
% % pooledMapAlphaCorr_Calcium imaging data_40_DSI0.30_1.00_created_May_15_2015_15_45_SFDS
% filedir=cd;
% %namesMTN=getFileList(filedir,'pooledMapAlphaCorr_Calcium imaging data_40_DSI0.30_1.00_created_Jul_16_2015_14_05_feature-based clustering_MTN',0,'anywhere');
% namesMTN=getFileList(filedir,'pooledMapAlphaCorr_Calcium imaging data_40_DSI0.30_1.00_created_Oct_08_2015_16_41__ONandONOFF_MTNd',0,'anywhere');
% MTN=load(num2str(namesMTN{1}));
% % standardize the length of MTN responses
% % determine the minimum number of elements
% 
% % standardize the length of responses
% for nu=1:numel(MTN.pooledmap.pResponse)
%     MTNnumPresponseT{nu}=MTN.pooledmap.pResponse{nu}(1:min(elements));
% end
% MTNpResponse=reshape(cell2mat(MTNnumPresponseT'),min(elements),[]);
% % normalize the responses
% for p=1:size(MTNpResponse,2)
%     MTNnorpResponse(:,p)=mat2gray(MTNpResponse(:,p));
% end
% % smooth the reponses
% for p=1:size(MTNpResponse,2)
%     MTNpResponseS(:,p)=smooth(MTNpResponse(:,p),15,'moving');
% end
% % normalize the smoothed responses
% MTNpResponseNorS=(MTNpResponseS-repmat(mean(MTNpResponseS(10:40,:),1),size(MTNpResponseS,1),1))./(repmat(max(MTNpResponseS,[],1),size(MTNpResponseS,1),1)-repmat(mean(MTNpResponseS(10:40,:),1),size(MTNpResponseS,1),1));
% 
% % switch between using the normalized responses or not
% switch NORMALRESstate
%     case 1
%         MTNresponse=MTNpResponseNorS;
%     case 0
%         MTNresponse=MTNpResponseS;
% end
% 
% 
% MTNonset1=MTN.pooledmap.onset1;
% 
% %% calculate parameters for clustering MTN cells
% figure;
% MTNpeakMinMag=MTNresponse(130,:);     %data point number 130 is where the minima occur for both ON and ONOFF cells
% for p=1:size(MTNresponse,2)
%     peakfinder(MTNresponse(:,p),0.1,0.1,1,false); drawnow expose
%     [peakLocMax,peakMagMax]=peakfinder(MTNresponse(:,p),0.1,0.1,1,false);
%     
%     MTNpeakLocMax{p}=peakLocMax;
%     MTNpeakMagMax{p}=peakMagMax;
%     
%     %%%% excludes from analysis every cell with a weak and noisy response (more than 4 detected peaks)
%     if size(MTNpeakLocMax{p},1)>4
%         MTNresponse(:,p)=0;
%     end
%     
%     %%% find OFF peak
%     offind=intersect(find(peakLocMax>145),find(peakLocMax<170));
%     if isempty(offind)
%         MTNpeakLocOff(p)=NaN;
%     else
%         MTNpeakLocOff(p)=peakLocMax(offind(1),1);
%     end
%     
%     %%% calculate slope of OFF peak
%     if ~isnan(MTNpeakLocOff(p))
%         yOffData=MTNresponse(MTNpeakLocOff(p):MTNpeakLocOff(p)+15,p);
%         % Set up fittype and options.
%         ft = fittype('poly1');
%         opts = fitoptions(ft);
%         % Fit model to data.
%         [fitresultOff, gofOff] = fit((1:size(yOffData,1))', yOffData(:,1), ft, opts );
%         MTNslopeOff(p)=fitresultOff.p1;
%     else
%         MTNslopeOff(p)=NaN;
%     end
%     
%     
%     %%% find ON peak
%     onind=find(peakLocMax>60);
%     if isempty(onind)
%         MTNpeakLocOn(p)=NaN;
%     else
%         MTNpeakLocOn(p)=peakLocMax(onind(1),1);
%     end
%     
%     %%% calculate slope of ON peak
%     if ~isnan(MTNpeakLocOn(p))
%         yOnData=MTNresponse(MTNpeakLocOn(p):MTNpeakLocOn(p)+slopePoints,p);
%         % Fit model to data.
%         [fitresultOn, gofOn] = fit((1:size(yOnData,1))', yOnData(:,1), ft, opts );
%         MTNslopeOn(p)=fitresultOn.p1;
%     else
%         MTNslopeOn(p)=NaN;
%     end
% end
% 
% MTNonSumResponse=sum(MTNresponse(60:130,:),1);
% MTNoffSumResponse=sum(MTNresponse(145:210,:),1);
% 
% MTNpeakLocOn=time(MTNpeakLocOn);      % convert DSpeakLocOn from number of data points to time
% MTNpeakLocOn=MTNpeakLocOn-MTNonset1;    % adjust DSpeakLocOn relative to the onset of the stimulus
% 
% MTNpeakLocOn=MTNpeakLocOn;
% MTNpeakLocOff=MTNpeakLocOff';
% MTNonSumResponse=MTNonSumResponse';
% MTNoffSumResponse=MTNoffSumResponse';
% MTNpeakMinMag=MTNpeakMinMag';
% MTNslopeOn=MTNslopeOn';
% MTNslopeOff=MTNslopeOff';
% 
% %data=[MTNpeakLocOn,MTNpeakLocOff,MTNonSumResponse,MTNoffSumResponse,MTNpeakMinMag,MTNslopeOn,MTNslopeOff];
% MTNdata=[MTNpeakLocOn,MTNslopeOn];
% 
% MTNidx = cluster(GM,MTNdata);
% MTNpost=posterior(GM,MTNdata);        %calculates the posterior probability of cells
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;scatter(1:length(MTNpost(MTNidx==1,1)),MTNpost(MTNidx==1,1),'b')
% figure;scatter(1:length(MTNpost(MTNidx==2,2)),MTNpost(MTNidx==2,2),'r')
% MTNindPost1=intersect(find(MTNidx==1),find(MTNpost(:,1)<postPcut));
% MTNindPost2=intersect(find(MTNidx==2),find(MTNpost(:,2)<postPcut));
% MTNcombIndPost=[MTNindPost1;MTNindPost2];
% % create MTNdataT, MTNinddataT, MTNidxT, and MTNnorpResponseT that include only cells with posterior probabilities higher than postPcut 
% % inddarta is the index to be used to recover the original indexing 
% MTNdataT=MTNdata;
% MTNdataT(MTNcombIndPost,:)=[];
% MTNinddata=(1:size(MTNdata,1))';
% MTNinddataT=MTNinddata;
% MTNinddataT(MTNcombIndPost,:)=[];
% MTNidxT=MTNidx;
% MTNidxT(MTNcombIndPost,:)=[];
% MTNnorpResponseT=MTNnorpResponse;
% MTNnorpResponseT(:,MTNcombIndPost)=[];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % cluster MTN cells
% MTNclust1=numel(find(MTNidx==1));
% MTNclust2=numel(find(MTNidx==2));
% if MTNclust2<MTNclust1
%     MTNmisclass=100*(MTNclust2/(MTNclust1+MTNclust2));
% else
%     MTNmisclass=100*(MTNclust1/(MTNclust1+MTNclust2));
% end
% 
% % plots individual and mean responses of all MTN retro cells DS cells clustered into 2
% % groups.
% %f17=figure;
% % set(f17,'position',[5 150 950 800]);
% % subplot(2,2,1);
% % plot(MTNnorpResponseT(:,MTNidxT==1),'k');
% % xlim([0 size(MTNnorpResponseT,1)]);
% % subplot(2,2,2);
% % plot(MTNnorpResponseT(:,MTNidxT==2),'k');
% % xlim([0 size(MTNnorpResponseT,1)]);
% % subplot(2,2,3);
% % shadedErrorBar(1:size(MTNnorpResponseT,1),mean(MTNnorpResponseT(:,MTNidxT==1),2),std(MTNnorpResponseT(:,MTNidxT==1),[],2),'k');
% % xlim([0 size(MTNnorpResponseT,1)]); ylim([0 1]);
% % title(['n = ',num2str(sum(MTNidxT==1))]);
% % subplot(2,2,4);
% % shadedErrorBar(1:size(MTNnorpResponseT,1),mean(MTNnorpResponseT(:,MTNidxT==2),2),std(MTNnorpResponseT(:,MTNidxT==2),[],2),'k');
% % xlim([0 size(MTNnorpResponseT,1)]); ylim([0 1]);
% % title(['n = ',num2str(sum(MTNidxT==2))]);
% % subtitle('MTN retro cells');
% 
% f17=figure('Color',[1 1 1],'Renderer','painters');
% set(f17,'position',[5 150 950 400]);
% subplot(1,2,1,'Parent',f17,'YColor',[1 1 1],'XColor',[1 1 1]);
% shadedErrorBar(time,mean(MTNnorpResponseT(:,MTNidxT==1),2),std(MTNnorpResponseT(:,MTNidxT==1),[],2),'k');
% hold on
% plot(time,mean(MTNnorpResponseT(:,MTNidxT==1),2),'r','LineWidth',1);
% xlim([0 time(end)]); ylim([0 1]);
% ylabel('\Delta\itF/F','FontSize',16);
% xlabel('Time (sec)','FontSize',16);
% title(['\itn\rm = ',num2str(sum(MTNidxT==1))],'FontSize',16);
% box off
% subplot(1,2,2,'Parent',f17,'YColor',[1 1 1],'XColor',[1 1 1]);
% shadedErrorBar(time,mean(MTNnorpResponseT(:,MTNidxT==2),2),std(MTNnorpResponseT(:,MTNidxT==2),[],2),'k');
% hold on
% plot(time,mean(MTNnorpResponseT(:,MTNidxT==2),2),'b','LineWidth',1);
% xlim([0 time(end)]); ylim([0 1]);
% xlabel('Time (sec)','FontSize',16);
% title(['\itn\rm = ',num2str(sum(MTNidxT==2))],'FontSize',16);
% box off
% 
% save('clustering workspace');
% 
% 
%% calculate the Gaussian mixture models for ON DS cells seperatly
% % with 2 to 10 clusters, 50 repetitions (every time it uses a different initial value), 
% % and determine the AIC, BIC, and negative log likelihood
% 
% % extract ON and ONOFF cells from the whole dataset
% if idx1>idx2
%     dataON=data(idx==2,:);
%     dataONOFF=data(idx==1,:);
%     DSnorpResponseON=DSnorpResponse(:,idx==2);
%     DSnorpResponseONOFF=DSnorpResponse(:,idx==1);
% else
%     dataON=data(idx==1,:);
%     dataONOFF=data(idx==2,:);
%     DSnorpResponseON=DSnorpResponse(:,idx==1);
%     DSnorpResponseONOFF=DSnorpResponse(:,idx==2);
% end
% 
% 
% options=statset('MaxIter',10000);
% for nK=2:1:10
% %GM=fitgmdist(data,nK,'RegularizationValue',0.001,'Replicates',50,'Options',options)
% GMon=fitgmdist(dataON,nK,'RegularizationValue',0.001,'Replicates',50,'Options',options)
% idxON(:,nK)=cluster(GMon,dataON);
% Pon{nK}=posterior(GMon,dataON);
% allLLon(nK)=GMon.NegativeLogLikelihood;
% end
% 
% % calculte bic and aic and deltaBIC and deltaAIC
% [aicON,bicON] = aicbic(-allLLon,1:10,size(dataON,1));
% for mb=2:10
% manualBICon(mb)=-2*-allLLon(mb)+mb*log(size(dataON,1));
% end
% 
% [Bon,Ion]=sort(bicON)
% for mb=1:9
% deltaBICon(mb)=Bon(1)-Bon(mb+1);
% end
% 
% [~,OptimalKon]=min(bicON(1:end))
% suboptimalKon=I(find(deltaBICon>-6)+1)     % additional models with deltaBIC smaller than 6 (that corresponds to a strong evidence against higher BIC)
% 
% figure;
% plot(2:10,bicON(2:10),'-or');
% figure;
% plot(2:10,manualBICon(2:10),'-ob');
% figure;
% plot(2:10,aicON(2:10),'-og');
% 
% 
% [BaicON,IaicON]=sort(aicON)
% for mb=1:9
% deltaAICon(mb)=BaicON(1)-BaicON(mb+1);
% end
% 
% [~,OptimalKaicON]=min(aicON(1:end))
% suboptimalKaicON=I(find(deltaAICon>-6)+1)     % additional models with deltaBIC smaller than 6 (that corresponds to a strong evidence against higher BIC)
% 
% 
% 
% OptimalKon=4;
% 
% % recalculate the Gaussian mixture models with the optimal number of
% % clusters (optimalK)
% %GM=fitgmdist(data,OptimalK,'RegularizationValue',0.001,'Replicates',50,'Options',options)
% GMon=fitgmdist(dataON,OptimalKon,'RegularizationValue',0.001,'Replicates',50,'Options',options)
% idxON=cluster(GMon,dataON);
% postON=posterior(GMon,dataON);        %calculates the posterior probability of cells
% 
% 
% % calculate the number of cells in each cluster
% idx1on=sum(idxON==1)
% idx2on=sum(idxON==2)
% idx3on=sum(idxON==3)
% idx4on=sum(idxON==4)
% idx5on=sum(idxON==5)
% 
% f11=figure;
% set(f11,'position',[5 150 950 800]);
% subplot(2,5,1);
% plot(DSnorpResponseON(:,idxON==1),'k');
% xlim([0 size(DSnorpResponseON,1)]);
% subplot(2,5,2);
% plot(DSnorpResponseON(:,idxON==2),'k');
% xlim([0 size(DSnorpResponseON,1)]);
% subplot(2,5,3);
% plot(DSnorpResponseON(:,idxON==3),'k');
% xlim([0 size(DSnorpResponseON,1)]);
% subplot(2,5,4);
% plot(DSnorpResponseON(:,idxON==4),'k');
% xlim([0 size(DSnorpResponseON,1)]);
% subplot(2,5,5);
% plot(DSnorpResponseON(:,idxON==5),'k');
% xlim([0 size(DSnorpResponseON,1)]);
% subplot(2,5,6);
% shadedErrorBar(1:size(DSnorpResponseON,1),mean(DSnorpResponseON(:,idxON==1),2),std(DSnorpResponseON(:,idxON==1),[],2),'k');
% xlim([0 size(DSnorpResponseON,1)]); ylim([0 1]);
% subplot(2,5,7);
% shadedErrorBar(1:size(DSnorpResponseON,1),mean(DSnorpResponseON(:,idxON==2),2),std(DSnorpResponseON(:,idxON==2),[],2),'k');
% xlim([0 size(DSnorpResponseON,1)]); ylim([0 1]);
% subplot(2,5,8);
% shadedErrorBar(1:size(DSnorpResponseON,1),mean(DSnorpResponseON(:,idxON==3),2),std(DSnorpResponseON(:,idxON==3),[],2),'k');
% xlim([0 size(DSnorpResponseON,1)]); ylim([0 1]);
% subplot(2,5,9);
% shadedErrorBar(1:size(DSnorpResponseON,1),mean(DSnorpResponseON(:,idxON==4),2),std(DSnorpResponseON(:,idxON==4),[],2),'k');
% xlim([0 size(DSnorpResponseON,1)]); ylim([0 1]);
% subplot(2,5,10);
% shadedErrorBar(1:size(DSnorpResponseON,1),mean(DSnorpResponseON(:,idxON==5),2),std(DSnorpResponseON(:,idxON==5),[],2),'k');
% xlim([0 size(DSnorpResponseON,1)]); ylim([0 1]);
% 
% 
% 
%% calculate Gaussian mixture models for ONOFF DS cells only
% options=statset('MaxIter',10000);
% for nK=2:1:10
% %GM=fitgmdist(data,nK,'RegularizationValue',0.001,'Replicates',50,'Options',options)
% GMonoff=fitgmdist(dataONOFF,nK,'RegularizationValue',0.001,'Replicates',50,'Options',options)
% idxONOFF(:,nK)=cluster(GMonoff,dataONOFF);
% Ponoff{nK}=posterior(GMonoff,dataONOFF);
% allLLonoff(nK)=GMonoff.NegativeLogLikelihood;
% end
% 
% % calculte bic and aic and deltaBIC and deltaAIC
% [aicONOFF,bicONOFF] = aicbic(-allLLonoff,1:10,size(dataONOFF,1));
% for mb=2:10
% manualBIConoff(mb)=-2*-allLLonoff(mb)+mb*log(size(dataONOFF,1));
% end
% 
% [Bonoff,Ionoff]=sort(bicONOFF)
% for mb=1:9
% deltaBIConoff(mb)=Bonoff(1)-Bonoff(mb+1);
% end
% 
% [~,OptimalKonoff]=min(bicONOFF(1:end))
% suboptimalKonoff=I(find(deltaBIConoff>-6)+1)     % additional models with deltaBIC smaller than 6 (that corresponds to a strong evidence against higher BIC)
% 
% figure;
% plot(2:10,bicONOFF(2:10),'-or');
% figure;
% plot(2:10,manualBIConoff(2:10),'-ob');
% figure;
% plot(2:10,aicONOFF(2:10),'-og');
% 
% 
% [BaicONOFF,IaicONOFF]=sort(aicONOFF)
% for mb=1:9
% deltaAIConoff(mb)=BaicONOFF(1)-BaicONOFF(mb+1);
% end
% 
% [~,OptimalKaicONOFF]=min(aicONOFF(1:end))
% suboptimalKaicONOFF=I(find(deltaAIConoff>-6)+1)     % additional models with deltaBIC smaller than 6 (that corresponds to a strong evidence against higher BIC)
% 
% 
% 
% OptimalKonoff=2;
% 
% % recalculate the Gaussian mixture models with the optimal number of
% % clusters (optimalK)
% %GM=fitgmdist(data,OptimalK,'RegularizationValue',0.001,'Replicates',50,'Options',options)
% GMonoff=fitgmdist(dataONOFF,OptimalKonoff,'RegularizationValue',0.001,'Replicates',50,'Options',options)
% idxONOFF=cluster(GMonoff,dataONOFF);
% postONOFF=posterior(GMonoff,dataONOFF);        %calculates the posterior probability of cells
% 
% 
% % calculate the number of cells in each cluster
% idx1onoff=sum(idxONOFF==1)
% idx2onoff=sum(idxONOFF==2)
% idx3onoff=sum(idxONOFF==3)
% idx4onoff=sum(idxONOFF==4)
% idx5onoff=sum(idxONOFF==5)
% 
% f11=figure;
% set(f11,'position',[5 150 950 800]);
% subplot(2,5,1);
% plot(DSnorpResponseONOFF(:,idxONOFF==1),'k');
% xlim([0 size(DSnorpResponseONOFF,1)]);
% subplot(2,5,2);
% plot(DSnorpResponseONOFF(:,idxONOFF==2),'k');
% xlim([0 size(DSnorpResponseONOFF,1)]);
% subplot(2,5,3);
% plot(DSnorpResponseONOFF(:,idxONOFF==3),'k');
% xlim([0 size(DSnorpResponseONOFF,1)]);
% subplot(2,5,4);
% plot(DSnorpResponseONOFF(:,idxONOFF==4),'k');
% xlim([0 size(DSnorpResponseONOFF,1)]);
% subplot(2,5,5);
% plot(DSnorpResponseONOFF(:,idxONOFF==5),'k');
% xlim([0 size(DSnorpResponseONOFF,1)]);
% subplot(2,5,6);
% shadedErrorBar(1:size(DSnorpResponseONOFF,1),mean(DSnorpResponseONOFF(:,idxONOFF==1),2),std(DSnorpResponseONOFF(:,idxONOFF==1),[],2),'k');
% xlim([0 size(DSnorpResponseONOFF,1)]); ylim([0 1]);
% subplot(2,5,7);
% shadedErrorBar(1:size(DSnorpResponseONOFF,1),mean(DSnorpResponseONOFF(:,idxONOFF==2),2),std(DSnorpResponseONOFF(:,idxONOFF==2),[],2),'k');
% xlim([0 size(DSnorpResponseONOFF,1)]); ylim([0 1]);
% subplot(2,5,8);
% shadedErrorBar(1:size(DSnorpResponseONOFF,1),mean(DSnorpResponseONOFF(:,idxONOFF==3),2),std(DSnorpResponseONOFF(:,idxONOFF==3),[],2),'k');
% xlim([0 size(DSnorpResponseONOFF,1)]); ylim([0 1]);
% subplot(2,5,9);
% shadedErrorBar(1:size(DSnorpResponseONOFF,1),mean(DSnorpResponseONOFF(:,idxONOFF==4),2),std(DSnorpResponseONOFF(:,idxONOFF==4),[],2),'k');
% xlim([0 size(DSnorpResponseONOFF,1)]); ylim([0 1]);
% subplot(2,5,10);
% shadedErrorBar(1:size(DSnorpResponseONOFF,1),mean(DSnorpResponseONOFF(:,idxONOFF==5),2),std(DSnorpResponseONOFF(:,idxONOFF==5),[],2),'k');
% xlim([0 size(DSnorpResponseONOFF,1)]); ylim([0 1]);
% 
% %% plot gaussian mixture graphs for ONOFF
% figure
% y=idxONOFF==1;
% h = gscatter(dataONOFF(:,1),dataONOFF(:,2),y);
% hold on
% xaxislim=xlim;
% yaxislim=ylim;
% 
% scatter(GMonoff.mu(1,1),GMonoff.mu(1,2),'ko')
% scatter(GMonoff.mu(2,1),GMonoff.mu(2,2),'ko')
% 
% % normalize sigma
% % maxsigma=max(GM.Sigma(:));
% % absGM.Sigma(:,:,1)=GM.Sigma(:,:,1)./maxsigma;
% % absGM.Sigma(:,:,2)=GM.Sigma(:,:,2)./maxsigma;
% 
% gauss1=gmdistribution(GMonoff.mu(1,:),GMonoff.Sigma(:,:,1));
% gauss2=gmdistribution(GMonoff.mu(2,:),GMonoff.Sigma(:,:,2));
% ezcontour(@(x1,x2)pdf(gauss1,[x1 x2]),xlim,ylim);
% ezcontour(@(x1,x2)pdf(gauss2,[x1 x2]),xlim,ylim);
% legend(h,'Model 0','Model1')
% hold off
% 
% %ezcontour(@(x1,x2)pdf(GM,[x1 x2]),xlim,ylim,7);
% 
% 
% figure;
% scatter(dataONOFF(idxONOFF==1,1),dataONOFF(idxONOFF==1,2),10,postONOFF(idxONOFF==1,1),'+')
% hold on
% scatter(dataONOFF(idxONOFF==2,1),dataONOFF(idxONOFF==2,2),10,postONOFF(idxONOFF==2,1),'o')
% hold off
% legend('Cluster 1','Cluster 2','Location','NW')
% clrmap = jet(80); colormap(clrmap(9:72,:))
% ylabel(colorbar,'Cluster 1 posterior probability')
% 
% figure;
% [~,order] = sort(postONOFF(:,1));
% plot(1:size(dataONOFF,1),postONOFF(order,1),'r-',1:size(dataONOFF,1),postONOFF(order,2),'b-');
% legend({'Cluster 1' 'Cluster 2'},'location','NW');
% ylabel('Cluster membership score');
% xlabel('Point ranking');
% hold on
% plot([0 2446],[0.95 0.95],'--k');
% plot([0 2446],[0.05 0.05],'--k');
% 

save('clustering summary');     % saves all variables in workspace and figures
end
