function [idxT,inddataT]=clusteringDS6(pResponse,onset1,NORMALRESstate,postPcut)

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
load('E:\Calcium imaging data\clusteringTest\mergedmeta_fov#11');
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

DSonSumResponse=sum(DSresponse(60:130,:),1);
DSoffSumResponse=sum(DSresponse(145:210,:),1);

DSpeakLocOn=time(DSpeakLocOn);      % convert DSpeakLocOn from number of data points to time
DSpeakLocOn=DSpeakLocOn-onset1;    % adjust DSpeakLocOn relative to the onset of the stimulus



% DSdiff=diff(DSresponse);
% [~,diffmax]=max(DSresponse(1:80,:));
% DSabsPeakLockOn=DSpeakLocOn-diffmax;

DSpeakLocOn=DSpeakLocOn;
DSpeakLocOff=DSpeakLocOff';
DSonSumResponse=DSonSumResponse';
DSoffSumResponse=DSoffSumResponse';
DSpeakMinMag=DSpeakMinMag';
DSslopeOn=DSslopeOn';
DSslopeOff=DSslopeOff';

%data=[DSpeakLocOn,DSpeakLocOff,DSonSumResponse,DSoffSumResponse,DSpeakMinMag,DSslopeOn,DSslopeOff];
data=[DSpeakLocOn,DSslopeOn];
%data=DSpeakLocOn;

figure;histogram(data(:,1),100);
xlabel('time to ON peak (data points)');

figure;histogram(data(:,2),100);
xlabel('slope following ON peak)');


%% calculate the Gaussian mixture models with 2 to 10 clusters, 50 repetitions (every time it uses a different initial value), 
% and determine the AIC, BIC, and negative log likelihood
options=statset('MaxIter',10000);
for nK=2:1:10
%GM=fitgmdist(data,nK,'RegularizationValue',0.001,'Replicates',50,'Options',options)
GM=fitgmdist(data,nK,'RegularizationValue',0.001,'Replicates',50,'Options',options)
idx(:,nK)=cluster(GM,data);
P{nK}=posterior(GM,data);
allLL(nK)=GM.NegativeLogLikelihood;
end

% calculte bic and aic and deltaBIC and deltaAIC
[aic,bic] = aicbic(-allLL,1:10,size(data,1));
for mb=2:10
manualBIC(mb)=-2*-allLL(mb)+mb*log(size(data,1));
end

[B,I]=sort(bic)
for mb=1:9
deltaBIC(mb)=B(1)-B(mb+1);
end

[~,OptimalK]=min(bic(1:end))
suboptimalK=I(find(deltaBIC>-6)+1);     % additional models with deltaBIC smaller than 6 (that corresponds to a strong evidence against higher BIC)

figure;
plot(2:10,bic(2:10),'-or');
figure;
plot(2:10,manualBIC(2:10),'-ob');


[Baic,Iaic]=sort(aic)
for mb=1:9
deltaAIC(mb)=Baic(1)-Baic(mb+1);
end

[~,OptimalKaic]=min(aic(1:end))
suboptimalKaic=I(find(deltaAIC>-6)+1);     % additional models with deltaBIC smaller than 6 (that corresponds to a strong evidence against higher BIC)



OptimalK=2;

% recalculate the Gaussian mixture models with the optimal number of
% clusters (optimalK)
%GM=fitgmdist(data,OptimalK,'RegularizationValue',0.001,'Replicates',50,'Options',options)
GM=fitgmdist(data,OptimalK,'RegularizationValue',0.001,'Replicates',50,'Options',options)
idx=cluster(GM,data);
post=posterior(GM,data);        %calculates the posterior probability of cells


% calculate the number of cells in each cluster
idx1=sum(idx==1)
idx2=sum(idx==2)
idx3=sum(idx==3)
idx4=sum(idx==4)
idx5=sum(idx==5)

f11=figure;
set(f11,'position',[5 150 950 800]);
subplot(2,5,1);
plot(DSnorpResponse(:,idx==1),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,5,2);
plot(DSnorpResponse(:,idx==2),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,5,3);
plot(DSnorpResponse(:,idx==3),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,5,4);
plot(DSnorpResponse(:,idx==4),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,5,5);
plot(DSnorpResponse(:,idx==5),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,5,6);
shadedErrorBar(1:size(DSnorpResponse,1),mean(DSnorpResponse(:,idx==1),2),std(DSnorpResponse(:,idx==1),[],2),'k');
xlim([0 size(DSnorpResponse,1)]); ylim([0 1]);
subplot(2,5,7);
shadedErrorBar(1:size(DSnorpResponse,1),mean(DSnorpResponse(:,idx==2),2),std(DSnorpResponse(:,idx==2),[],2),'k');
xlim([0 size(DSnorpResponse,1)]); ylim([0 1]);
subplot(2,5,8);
shadedErrorBar(1:size(DSnorpResponse,1),mean(DSnorpResponse(:,idx==3),2),std(DSnorpResponse(:,idx==3),[],2),'k');
xlim([0 size(DSnorpResponse,1)]); ylim([0 1]);
subplot(2,5,9);
shadedErrorBar(1:size(DSnorpResponse,1),mean(DSnorpResponse(:,idx==4),2),std(DSnorpResponse(:,idx==4),[],2),'k');
xlim([0 size(DSnorpResponse,1)]); ylim([0 1]);
subplot(2,5,10);
shadedErrorBar(1:size(DSnorpResponse,1),mean(DSnorpResponse(:,idx==5),2),std(DSnorpResponse(:,idx==5),[],2),'k');
xlim([0 size(DSnorpResponse,1)]); ylim([0 1]);

%% plot gaussian mixture graphs
figure
y=idx==1;
h = gscatter(data(:,1),data(:,2),y);
hold on
xaxislim=xlim;
yaxislim=ylim;

scatter(GM.mu(1,1),GM.mu(1,2),'ko')
scatter(GM.mu(2,1),GM.mu(2,2),'ko')

% normalize sigma
% maxsigma=max(GM.Sigma(:));
% absGM.Sigma(:,:,1)=GM.Sigma(:,:,1)./maxsigma;
% absGM.Sigma(:,:,2)=GM.Sigma(:,:,2)./maxsigma;

gauss1=gmdistribution(GM.mu(1,:),GM.Sigma(:,:,1));
gauss2=gmdistribution(GM.mu(2,:),GM.Sigma(:,:,2));
ezcontour(@(x1,x2)pdf(gauss1,[x1 x2]),xlim,ylim);
ezcontour(@(x1,x2)pdf(gauss2,[x1 x2]),xlim,ylim);
legend(h,'Model 0','Model1')
hold off

%ezcontour(@(x1,x2)pdf(GM,[x1 x2]),xlim,ylim,7);


figure;
scatter(data(idx==1,1),data(idx==1,2),10,post(idx==1,1),'+')
hold on
scatter(data(idx==2,1),data(idx==2,2),10,post(idx==2,1),'o')
hold off
legend('Cluster 1','Cluster 2','Location','NW')
clrmap = jet(80); colormap(clrmap(9:72,:))
ylabel(colorbar,'Cluster 1 posterior probability')

figure;
[~,order] = sort(post(:,1));
plot(1:size(data,1),post(order,1),'r-',1:size(data,1),post(order,2),'b-');
legend({'Cluster 1' 'Cluster 2'},'location','NW');
ylabel('Cluster membership score');
xlabel('Point ranking');
hold on
plot([0 2446],[0.95 0.95],'--k');
plot([0 2446],[0.05 0.05],'--k');



%% apply the posterior probability cutoff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots the posterior probability of cells belonging to each of the
% two first clusters
% uses postPcut as the cutoff for posterior probablity; cutoff of 0.5 means no restriction. 
figure;scatter(1:length(post(idx==1,1)),post(idx==1,1),'b')
figure;scatter(1:length(post(idx==2,2)),post(idx==2,2),'r')
% find the cells with posterior probability low than 0.6
% indPost1=find(post(idx==1,1)<0.6);
% indPost2=find(post(idx==2,2)<0.6);
indPost1=intersect(find(idx==1),find(post(:,1)<postPcut));
indPost2=intersect(find(idx==2),find(post(:,2)<postPcut));
combIndPost=[indPost1;indPost2];
% create dataT, inddataT, idxT, and DSnorpResponseT that include only cells with posterior probabilities higher than postPcut 
% inddarta is the 
dataT=data;
dataT(combIndPost,:)=[];
inddata=(1:size(data,1))';
inddataT=inddata;
inddataT(combIndPost,:)=[];
idxT=idx;
idxT(combIndPost,:)=[];
DSnorpResponseT=DSnorpResponse;
DSnorpResponseT(:,combIndPost)=[];
postT=post;
postT(combIndPost,:)=[];

figure;
scatter(dataT(idxT==1,1),dataT(idxT==1,2),10,postT(idxT==1,1),'+')
hold on
scatter(dataT(idxT==2,1),dataT(idxT==2,2),10,postT(idxT==2,1),'o')
hold off
legend('Cluster 1','Cluster 2','Location','NW')
clrmap = jet(80); colormap(clrmap(9:72,:))
ylabel(colorbar,'Cluster 1 posterior probability')


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
set(f12,'position',[5 150 950 400]);
subplot(1,2,1,'Parent',f12,'YColor',[1 1 1],'XColor',[1 1 1]);
shadedErrorBar(time,mean(DSnorpResponseT(:,idxT==1),2),std(DSnorpResponseT(:,idxT==1),[],2),'k');
hold on
plot(time,mean(DSnorpResponseT(:,idxT==1),2),'r','LineWidth',1);
xlim([0 time(end)]); ylim([0 1]);
ylabel('\Delta\itF/F','FontSize',16);
xlabel('Time (sec)','FontSize',16);
title(['\itn\rm = ',num2str(sum(idxT==1))],'FontSize',16);
box off
subplot(1,2,2,'Parent',f12,'YColor',[1 1 1],'XColor',[1 1 1]);
shadedErrorBar(time,mean(DSnorpResponseT(:,idxT==2),2),std(DSnorpResponseT(:,idxT==2),[],2),'k');
hold on
plot(time,mean(DSnorpResponseT(:,idxT==2),2),'b','LineWidth',1);
xlim([0 time(end)]); ylim([0 1]);
xlabel('Time (sec)','FontSize',16);
title(['\itn\rm = ',num2str(sum(idxT==2))],'FontSize',16);
box off


%% cluster new data sets - SNS cells (10 ONDS and 10 ONOFFDS)
% extract responses from the file: summary_All SNS cells
filedir=cd;
namesSNS=getFileList(filedir,'summary_All SNS cells',0,'anywhere');
load(num2str(namesSNS{1}));

SNSpResponseONDS=summary.pResponseONDS(1:min(elements),:);
SNSpResponseONOFFDS=summary.pResponseONOFFDS(1:min(elements),:);
% standardize the length of responses
SNSpResponseONDS(min(elements):size(SNSpResponseONDS,1),:)=[];
SNSpResponseONOFFDS(min(elements):size(SNSpResponseONOFFDS,1),:)=[];
% normalize the responses
for p=1:size(SNSpResponseONDS,2)
    SNSpResponseONDSnor(:,p)=mat2gray(SNSpResponseONDS(:,p));
end
for p=1:size(SNSpResponseONOFFDS,2)
    SNSpResponseONOFFDSnor(:,p)=mat2gray(SNSpResponseONOFFDS(:,p));
end
% smooth the responses
for p=1:size(SNSpResponseONDS,2)
    SNSpResponseONDSs(:,p)=smooth(SNSpResponseONDS(:,p),15,'moving');
end
for p=1:size(SNSpResponseONOFFDS,2)
    SNSpResponseONOFFDSs(:,p)=smooth(SNSpResponseONOFFDS(:,p),15,'moving');
end
% normalize the smoothed responses
SNSpResponseONDSnorS=(SNSpResponseONDSs-repmat(mean(SNSpResponseONDSs(10:40,:),1),size(SNSpResponseONDSs,1),1))./(repmat(max(SNSpResponseONDSs,[],1),size(SNSpResponseONDSs,1),1)-repmat(mean(SNSpResponseONDSs(10:40,:),1),size(SNSpResponseONDSs,1),1));
SNSpResponseONOFFDSnorS=(SNSpResponseONOFFDSs-repmat(mean(SNSpResponseONOFFDSs(10:40,:),1),size(SNSpResponseONOFFDSs,1),1))./(repmat(max(SNSpResponseONOFFDSs,[],1),size(SNSpResponseONOFFDSs,1),1)-repmat(mean(SNSpResponseONOFFDSs(10:40,:),1),size(SNSpResponseONOFFDSs,1),1));
SNSpResponseS=[SNSpResponseONDSs,SNSpResponseONOFFDSs];
SNSpResponsenorS=[SNSpResponseONDSnorS,SNSpResponseONOFFDSnorS];
SNSonset1=[summary.onset1ONDS,summary.onset1ONOFFDS];

% switch between using the normalized responses or not
switch NORMALRESstate
    case 1
        SNSresponse=[SNSpResponseONDSnorS,SNSpResponseONOFFDSnorS];
    case 0
        SNSresponse=[SNSpResponseONDSs,SNSpResponseONOFFDSs];
end


%% calculate parameters for clustering SNS cells
figure;
SNSpeakMinMag=SNSresponse(130,:);     %data point number 130 is where the minima occur for both ON and ONOFF cells
for p=1:size(SNSresponse,2)
    peakfinder(SNSresponse(:,p),0.1,0.1,1,false); drawnow expose
    [peakLocMax,peakMagMax]=peakfinder(SNSresponse(:,p),0.1,0.1,1,false);
    
    SNSpeakLocMax{p}=peakLocMax;
    SNSpeakMagMax{p}=peakMagMax;
    
    %%%% excludes from analysis every cell with a weak and noisy response (more than 4 detected peaks)
    if size(SNSpeakLocMax{p},1)>4
        SNSresponse(:,p)=0;
    end
    
    %%% find OFF peak
    offind=intersect(find(peakLocMax>145),find(peakLocMax<170));
    if isempty(offind)
        SNSpeakLocOff(p)=NaN;
    else
        SNSpeakLocOff(p)=peakLocMax(offind(1),1);
    end
    
    %%% calculate slope of OFF peak
    if ~isnan(SNSpeakLocOff(p))
        yOffData=SNSresponse(SNSpeakLocOff(p):SNSpeakLocOff(p)+15,p);
        % Set up fittype and options.
        ft = fittype('poly1');
        opts = fitoptions(ft);
        % Fit model to data.
        [fitresultOff, gofOff] = fit((1:size(yOffData,1))', yOffData(:,1), ft, opts );
        SNSslopeOff(p)=fitresultOff.p1;
    else
        SNSslopeOff(p)=NaN;
    end
    
    
    %%% find ON peak
    onind=find(peakLocMax>60);
    if isempty(onind)
        SNSpeakLocOn(p)=NaN;
    else
        SNSpeakLocOn(p)=peakLocMax(onind(1),1);
    end
    
    %%% calculate slope of ON peak
    if ~isnan(SNSpeakLocOn(p))
        yOnData=SNSresponse(SNSpeakLocOn(p):SNSpeakLocOn(p)+slopePoints,p);
        % Fit model to data.
        [fitresultOn, gofOn] = fit((1:size(yOnData,1))', yOnData(:,1), ft, opts );
        SNSslopeOn(p)=fitresultOn.p1;
    else
        SNSslopeOn(p)=NaN;
    end
end

SNSonSumResponse=sum(SNSresponse(60:130,:),1);
SNSoffSumResponse=sum(SNSresponse(145:210,:),1);

SNSpeakLocOn=time(SNSpeakLocOn);      % convert DSpeakLocOn from number of data points to time
SNSpeakLocOn=SNSpeakLocOn-SNSonset1';    % adjust DSpeakLocOn relative to the onset of the stimulus

SNSpeakLocOn=SNSpeakLocOn;
SNSpeakLocOff=SNSpeakLocOff';
SNSonSumResponse=SNSonSumResponse';
SNSoffSumResponse=SNSoffSumResponse';
SNSpeakMinMag=SNSpeakMinMag';
SNSslopeOn=SNSslopeOn';
SNSslopeOff=SNSslopeOff';

%data=[SNSpeakLocOn,SNSpeakLocOff,SNSonSumResponse,SNSoffSumResponse,SNSpeakMinMag,SNSslopeOn,SNSslopeOff];
SNSdata=[SNSpeakLocOn,SNSslopeOn];



SNSidx = cluster(GM,SNSdata);
SNSpost=posterior(GM,SNSdata);        %calculates the posterior probability of cells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;scatter(1:length(SNSpost(SNSidx==1,1)),SNSpost(SNSidx==1,1),'b')
figure;scatter(1:length(SNSpost(SNSidx==2,2)),SNSpost(SNSidx==2,2),'r')
SNSindPost1=intersect(find(SNSidx==1),find(SNSpost(:,1)<postPcut));
SNSindPost2=intersect(find(SNSidx==2),find(SNSpost(:,2)<postPcut));
SNScombIndPost=[SNSindPost1;SNSindPost2];
% create dataT, inddataT, SNSidxT, and DSnorpResponseT that include only cells with posterior probabilities higher than postPcut 
% inddarta is the index to be used to recover the original indexing 
SNSdataT=SNSdata;
SNSdataT(SNScombIndPost,:)=[];
SNSinddata=(1:size(SNSdata,1))';
SNSinddataT=SNSinddata;
SNSinddataT(SNScombIndPost,:)=[];
SNSidxT=SNSidx;
SNSidxT(SNScombIndPost,:)=[];
SNSpResponsenorST=SNSpResponsenorS;
SNSpResponsenorST(:,SNScombIndPost)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% cluster SNS cells
ONDSn=10;
ONOFFDSn=10;
SNSONclust1=0;
SNSONclust2=0;
SNSONOFFclust1=0;
SNSONOFFclust2=0;
for s=1:size(SNSinddataT,1)
    if SNSinddataT(s)<=ONDSn
        if SNSidxT(s)==1
            SNSONclust1=SNSONclust1+1;
        elseif SNSidxT(s)==2
            SNSONclust2=SNSONclust2+1;
        end
    end
    if SNSinddataT(s)>ONDSn
        if SNSidxT(s)==1
            SNSONOFFclust1=SNSONOFFclust1+1;
        elseif SNSidxT(s)==2
            SNSONOFFclust2=SNSONOFFclust2+1;
        end
    end
end

% if SNS2clust2<SNS2clust1
%     SNSmisclass=100*((SNS1clust1+SNS2clust2)/(SNS1clust1+SNS2clust2+SNS1clust2+SNS2clust1));
% else
%     SNSmisclass=100*((SNS1clust2+SNS2clust1)/(SNS1clust1+SNS2clust2+SNS1clust2+SNS2clust1));
% end

% plots the clustered SNS cells responses
% f13=figure;
% set(f13,'position',[5 150 950 800]);
% subplot(2,2,1);
% plot(SNSpResponsenorST(:,SNSidxT==1),'k');
% xlim([0 size(SNSpResponsenorST,1)]); ylim([0 1]);
% subplot(2,2,2);
% plot(SNSpResponsenorST(:,SNSidxT==2),'k');
% xlim([0 size(SNSpResponsenorST,1)]); ylim([0 1]);
subplot(2,2,3);
shadedErrorBar(1:size(SNSpResponsenorST,1),mean(SNSpResponsenorST(:,SNSidxT==1),2),std(SNSpResponsenorST(:,SNSidxT==1),[],2),'k');
xlim([0 size(SNSpResponsenorST,1)]); ylim([0 1]);
title(['n = ',num2str(sum(SNSidxT==1))]);
subplot(2,2,4);
shadedErrorBar(1:size(SNSpResponsenorST,1),mean(SNSpResponsenorST(:,SNSidxT==2),2),std(SNSpResponsenorST(:,SNSidxT==2),[],2),'k');
xlim([0 size(SNSpResponsenorST,1)]); ylim([0 1]);
title(['n = ',num2str(sum(SNSidxT==2))]);
% subtitle('SNS cells');

f13=figure('Color',[1 1 1],'Renderer','painters');
set(f13,'position',[5 150 950 400]);
subplot(1,2,1,'Parent',f13,'YColor',[1 1 1],'XColor',[1 1 1]);
shadedErrorBar(time(1:end-1,1),mean(SNSpResponsenorST(:,SNSidxT==1),2),std(SNSpResponsenorST(:,SNSidxT==1),[],2),'k');
hold on
plot(time(1:end-1,1),mean(SNSpResponsenorST(:,SNSidxT==1),2),'r','LineWidth',1);
xlim([0 time(end-1)]); ylim([0 1]);
ylabel('\Delta\itF/F','FontSize',16);
xlabel('Time (sec)','FontSize',16);
title(['\itn\rm = ',num2str(sum(SNSidxT==1))],'FontSize',16);
box off
subplot(1,2,2,'Parent',f13,'YColor',[1 1 1],'XColor',[1 1 1]);
shadedErrorBar(time(1:end-1,1),mean(SNSpResponsenorST(:,SNSidxT==2),2),std(SNSpResponsenorST(:,SNSidxT==2),[],2),'k');
hold on
plot(time(1:end-1,1),mean(SNSpResponsenorST(:,SNSidxT==2),2),'b','LineWidth',1);
xlim([0 time(end-1)]); ylim([0 1]);
xlabel('Time (sec)','FontSize',16);
title(['\itn\rm = ',num2str(sum(SNSidxT==2))],'FontSize',16);
box off



end
