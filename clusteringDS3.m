function [idxT,inddataT]=clusteringDS3(pResponse,NORMALRESstate,postPcut)

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

%% PCA analysis
[coeff,score,latent,~,explained] = pca(DSresponse);
% reconstructing the responses
s=zeros(size(score,1),size(coeff,1));
for u=1:size(coeff,1)
for v=1:size(coeff,2)
s(:,u)=s(:,u)+coeff(u,v).*score(:,v);
end
end

stdscore=std(score,[],1);

data=[stdscore(1).*coeff(:,1),stdscore(2).*coeff(:,2),stdscore(3).*coeff(:,3),stdscore(4).*coeff(:,4),stdscore(5).*coeff(:,5)];


%% calculate parameters for clustering
figure;
peakMinMag=DSresponse(130,:);     %data point number 130 is where the minima occur for both ON and ONOFF cells
for p=1:size(DSresponse,2)
    peakfinder(DSresponse(:,p),0.1,0.1,1,false); drawnow expose
    [peakLocMax,peakMagMax]=peakfinder(DSresponse(:,p),0.1,0.1,1,false);
    
    DSpeakLocMax{p}=peakLocMax;
    DSpeakMagMax{p}=peakMagMax;
    
    %%%% excludes from analysis every cell with a weak and noisy response (more than 4 detected peaks)
    if size(DSpeakLocMax{p},1)>4
        DSresponse(:,p)=0;
    end
    
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
        slopeOff(p)=fitresultOff.p1;
    else
        slopeOff(p)=NaN;
    end
    
    
    %%% find ON peak
    onind=find(peakLocMax>60);
    if isempty(onind)
        DSpeakLocOn(p)=NaN;
    else
        DSpeakLocOn(p)=peakLocMax(onind(1),1);
    end
    
    %%% calculate slope of ON peak
    if ~isnan(DSpeakLocOn(p))
        yOnData=DSresponse(DSpeakLocOn(p):DSpeakLocOn(p)+15,p);
        % Fit model to data.
        [fitresultOn, gofOn] = fit((1:size(yOnData,1))', yOnData(:,1), ft, opts );
        slopeOn(p)=fitresultOn.p1;
    else
        slopeOn(p)=NaN;
    end
end

onSumResponse=sum(DSresponse(60:130,:),1);
offSumResponse=sum(DSresponse(145:210,:),1);

DSpeakLocOn=DSpeakLocOn';
DSpeakLocOff=DSpeakLocOff';
onSumResponse=onSumResponse';
offSumResponse=offSumResponse';
peakMinMag=peakMinMag';
slopeOn=slopeOn';
slopeOff=slopeOff';

%data=[DSpeakLocOn,DSpeakLocOff,onSumResponse,offSumResponse,peakMinMag,slopeOn,slopeOff];
%data=[DSpeakLocOn,peakMinMag];
data=DSpeakLocOn;


%% calculate the Gaussian mixture models with 2 to 5 clusters, 50 repetitions (every time it uses a different initial value), 
% and determine the AIC, BIC, and negative log likelihood
options=statset('MaxIter',10000);
% for nK=2:1:5
% %GM=fitgmdist(data,nK,'RegularizationValue',0.001,'Replicates',50,'Options',options)
% GM=fitgmdist(data,nK,'RegularizationValue',0.001,'Replicates',50,'Options',options)
% idx(:,nK)=cluster(GM,data);
% P{nK}=posterior(GM,data);
% allAIC(nK)=GM.AIC;
% allBIC(nK)=GM.BIC;
% allLL(nK)=GM.NegativeLogLikelihood;
% end
% [~,OptimalK]=min(allBIC(2:end));

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

% plots individual and mean responses of all DS cells clustered into 2
% groups.
f12=figure;
set(f12,'position',[5 150 950 800]);
subplot(2,2,1);
plot(DSnorpResponseT(:,idxT==1),'k');
xlim([0 size(DSnorpResponseT,1)]);
subplot(2,2,2);
plot(DSnorpResponseT(:,idxT==2),'k');
xlim([0 size(DSnorpResponseT,1)]);
subplot(2,2,3);
shadedErrorBar(1:size(DSnorpResponseT,1),mean(DSnorpResponseT(:,idxT==1),2),std(DSnorpResponseT(:,idxT==1),[],2),'k');
xlim([0 size(DSnorpResponseT,1)]); ylim([0 1]);
title(['n = ',num2str(sum(idxT==1))]);
subplot(2,2,4);
shadedErrorBar(1:size(DSnorpResponseT,1),mean(DSnorpResponseT(:,idxT==2),2),std(DSnorpResponseT(:,idxT==2),[],2),'k');
xlim([0 size(DSnorpResponseT,1)]); ylim([0 1]);
title(['n = ',num2str(sum(idxT==2))]);
subtitle('all DS cells');


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

% switch between using the normalized responses or not
switch NORMALRESstate
    case 1
        SNSresponse=[SNSpResponseONDSnorS,SNSpResponseONOFFDSnorS];
    case 0
        SNSresponse=[SNSpResponseONDSs,SNSpResponseONOFFDSs];
end

% PCA analysis on SNS cells
[SNScoeff,SNSscore,SNSlatent,~,SNSexplained] = pca(SNSresponse);
% reconstructing the responses
SNSs=zeros(size(SNSscore,1),size(SNScoeff,1));
for u=1:size(SNScoeff,1)
for v=1:size(SNScoeff,2)
SNSs(:,u)=SNSs(:,u)+SNScoeff(u,v).*SNSscore(:,v);
end
end

SNSstdscore=std(SNSscore,[],1);

SNSdata=[SNSstdscore(1).*SNScoeff(:,1),SNSstdscore(2).*SNScoeff(:,2),SNSstdscore(3).*SNScoeff(:,3),SNSstdscore(4).*SNScoeff(:,4),SNSstdscore(5).*SNScoeff(:,5)];

% SNSdata=[SNSstdscore(1).*SNScoeff(:,1),SNSstdscore(2).*SNScoeff(:,2),SNSstdscore(3).*SNScoeff(:,3),SNSstdscore(4).*SNScoeff(:,4),SNSstdscore(5).*SNScoeff(:,5),...
%     SNSstdscore(6).*SNScoeff(:,6),SNSstdscore(7).*SNScoeff(:,7),SNSstdscore(8).*SNScoeff(:,8),SNSstdscore(9).*SNScoeff(:,9),SNSstdscore(10).*SNScoeff(:,10),...
%     SNSstdscore(11).*SNScoeff(:,11),SNSstdscore(12).*SNScoeff(:,12),SNSstdscore(13).*SNScoeff(:,13),SNSstdscore(14).*SNScoeff(:,14),SNSstdscore(15).*SNScoeff(:,15)];

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
f13=figure;
set(f13,'position',[5 150 950 800]);
subplot(2,2,1);
plot(SNSpResponsenorST(:,SNSidxT==1),'k');
xlim([0 size(SNSpResponsenorST,1)]); ylim([0 1]);
subplot(2,2,2);
plot(SNSpResponsenorST(:,SNSidxT==2),'k');
xlim([0 size(SNSpResponsenorST,1)]); ylim([0 1]);
subplot(2,2,3);
shadedErrorBar(1:size(SNSpResponsenorST,1),mean(SNSpResponsenorST(:,SNSidxT==1),2),std(SNSpResponsenorST(:,SNSidxT==1),[],2),'k');
xlim([0 size(SNSpResponsenorST,1)]); ylim([0 1]);
title(['n = ',num2str(sum(SNSidxT==1))]);
subplot(2,2,4);
shadedErrorBar(1:size(SNSpResponsenorST,1),mean(SNSpResponsenorST(:,SNSidxT==2),2),std(SNSpResponsenorST(:,SNSidxT==2),[],2),'k');
xlim([0 size(SNSpResponsenorST,1)]); ylim([0 1]);
title(['n = ',num2str(sum(SNSidxT==2))]);
subtitle('SNS cells');


%% cluster new data sets - CART ONOFFDS cells
% extract responses from the following file:
% pooledMapAlphaCorr_Calcium imaging data_40_DSI0.30_1.00_created_May_15_2015_14_47_CARTONOFFDS
filedir=cd;
namesCART=getFileList(filedir,'pooledMapAlphaCorr_Calcium imaging data_40_DSI0.30_1.00_created_Jul_15_2015_14_28_magCut1.5_unsupervised_CART',0,'anywhere');
CART=load(num2str(namesCART{1}));
% standardize the length of CART responses
% determine the minimum number of elements

% standardize the length of responses
for nu=1:numel(CART.pooledmap.pResponse)
    CARTnumPresponseT{nu}=CART.pooledmap.pResponse{nu}(1:min(elements));
end
CARTpResponse=reshape(cell2mat(CARTnumPresponseT'),min(elements),[]);
% normalize the responses
for p=1:size(CARTpResponse,2)
    CARTnorpResponse(:,p)=mat2gray(CARTpResponse(:,p));
end
% smooth the reponses
for p=1:size(CARTpResponse,2)
    CARTpResponseS(:,p)=smooth(CARTpResponse(:,p),15,'moving');
end
% normalize the smoothed responses
CARTpResponseNorS=(CARTpResponseS-repmat(mean(CARTpResponseS(10:40,:),1),size(CARTpResponseS,1),1))./(repmat(max(CARTpResponseS,[],1),size(CARTpResponseS,1),1)-repmat(mean(CARTpResponseS(10:40,:),1),size(CARTpResponseS,1),1));

% switch between using the normalized responses or not
switch NORMALRESstate
    case 1
        CARTresponse=CARTpResponseNorS;
    case 0
        CARTresponse=CARTpResponseS;
end

% PCA analysis on CART cells
[CARTcoeff,CARTscore,CARTlatent,~,CARTexplained] = pca(CARTresponse);
% reconstructing the responses
CARTs=zeros(size(CARTscore,1),size(CARTcoeff,1));
for u=1:size(CARTcoeff,1)
for v=1:size(CARTcoeff,2)
CARTs(:,u)=CARTs(:,u)+CARTcoeff(u,v).*CARTscore(:,v);
end
end

CARTstdscore=std(CARTscore,[],1);

CARTdata=[CARTstdscore(1).*CARTcoeff(:,1),CARTstdscore(2).*CARTcoeff(:,2),CARTstdscore(3).*CARTcoeff(:,3),CARTstdscore(4).*CARTcoeff(:,4),CARTstdscore(5).*CARTcoeff(:,5)];

% CARTdata=[CARTstdscore(1).*CARTcoeff(:,1),CARTstdscore(2).*CARTcoeff(:,2),CARTstdscore(3).*CARTcoeff(:,3),CARTstdscore(4).*CARTcoeff(:,4),CARTstdscore(5).*CARTcoeff(:,5),...
%     CARTstdscore(6).*CARTcoeff(:,6),CARTstdscore(7).*CARTcoeff(:,7),CARTstdscore(8).*CARTcoeff(:,8),CARTstdscore(9).*CARTcoeff(:,9),CARTstdscore(10).*CARTcoeff(:,10),...
%     CARTstdscore(11).*CARTcoeff(:,11),CARTstdscore(12).*CARTcoeff(:,12),CARTstdscore(13).*CARTcoeff(:,13),CARTstdscore(14).*CARTcoeff(:,14),CARTstdscore(15).*CARTcoeff(:,15)];

CARTidx = cluster(GM,CARTdata);
CARTpost=posterior(GM,CARTdata);        %calculates the posterior probability of cells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;scatter(1:length(CARTpost(CARTidx==1,1)),CARTpost(CARTidx==1,1),'b')
figure;scatter(1:length(CARTpost(CARTidx==2,2)),CARTpost(CARTidx==2,2),'r')
CARTindPost1=intersect(find(CARTidx==1),find(CARTpost(:,1)<postPcut));
CARTindPost2=intersect(find(CARTidx==2),find(CARTpost(:,2)<postPcut));
CARTcombIndPost=[CARTindPost1;CARTindPost2];
% create CARTdataT, CARTinddataT, CARTidxT, and CARTnorpResponseT that include only cells with posterior probabilities higher than postPcut 
% inddarta is the index to be used to recover the original indexing 
CARTdataT=CARTdata;
CARTdataT(CARTcombIndPost,:)=[];
CARTinddata=(1:size(CARTdata,1))';
CARTinddataT=CARTinddata;
CARTinddataT(CARTcombIndPost,:)=[];
CARTidxT=CARTidx;
CARTidxT(CARTcombIndPost,:)=[];
CARTnorpResponseT=CARTnorpResponse;
CARTnorpResponseT(:,CARTcombIndPost)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cluster CART cells
CARTclust1=numel(find(CARTidxT==1));
CARTclust2=numel(find(CARTidxT==2));
if CARTclust2<CARTclust1
    CARTmisclass=100*(CARTclust2/(CARTclust1+CARTclust2));
else
    CARTmisclass=100*(CARTclust1/(CARTclust1+CARTclust2));
end
    
% plots individual and mean responses of all CART positive DS cells clustered into 2
% groups.
f14=figure;
set(f14,'position',[5 150 950 800]);
subplot(2,2,1);
plot(CARTnorpResponseT(:,CARTidxT==1),'k');
xlim([0 size(CARTnorpResponseT,1)]);
subplot(2,2,2);
plot(CARTnorpResponseT(:,CARTidxT==2),'k');
xlim([0 size(CARTnorpResponseT,1)]);
subplot(2,2,3);
shadedErrorBar(1:size(CARTnorpResponseT,1),mean(CARTnorpResponseT(:,CARTidxT==1),2),std(CARTnorpResponseT(:,CARTidxT==1),[],2),'k');
xlim([0 size(CARTnorpResponseT,1)]); ylim([0 1]);
title(['n = ',num2str(sum(CARTidxT==1))]);
subplot(2,2,4);
shadedErrorBar(1:size(CARTnorpResponseT,1),mean(CARTnorpResponseT(:,CARTidxT==2),2),std(CARTnorpResponseT(:,CARTidxT==2),[],2),'k');
xlim([0 size(CARTnorpResponseT,1)]); ylim([0 1]);
title(['n = ',num2str(sum(CARTidxT==2))]);
subtitle('CART positive');

%% cluster new data sets - MTN and SF retro DS cells (both ON and ONOFF were included)
% extract responses from the following file:
% pooledMapAlphaCorr_Calcium imaging data_40_DSI0.30_1.00_created_May_15_2015_15_30_MTNSFDS
filedir=cd;
namesMTNSF=getFileList(filedir,'pooledMapAlphaCorr_Calcium imaging data_40_DSI0.30_1.00_created_Jul_15_2015_14_29_magCut1.5_unsupervised_MTNSF',0,'anywhere');
MTNSF=load(num2str(namesMTNSF{1}));
% standardize the length of MTNSF responses
% determine the minimum number of elements

% standardize the length of responses
for nu=1:numel(MTNSF.pooledmap.pResponse)
    MTNSFnumPresponseT{nu}=MTNSF.pooledmap.pResponse{nu}(1:min(elements));
end
MTNSFpResponse=reshape(cell2mat(MTNSFnumPresponseT'),min(elements),[]);
% normalize the responses
for p=1:size(MTNSFpResponse,2)
    MTNSFnorpResponse(:,p)=mat2gray(MTNSFpResponse(:,p));
end
% smooth the reponses
for p=1:size(MTNSFpResponse,2)
    MTNSFpResponseS(:,p)=smooth(MTNSFpResponse(:,p),15,'moving');
end
% normalize the smoothed responses
MTNSFpResponseNorS=(MTNSFpResponseS-repmat(mean(MTNSFpResponseS(10:40,:),1),size(MTNSFpResponseS,1),1))./(repmat(max(MTNSFpResponseS,[],1),size(MTNSFpResponseS,1),1)-repmat(mean(MTNSFpResponseS(10:40,:),1),size(MTNSFpResponseS,1),1));

% switch between using the normalized responses or not
switch NORMALRESstate
    case 1
        MTNSFresponse=MTNSFpResponseNorS;
    case 0
        MTNSFresponse=MTNSFpResponseS;
end

% PCA analysis on MTNSF cells
[MTNSFcoeff,MTNSFscore,MTNSFlatent,~,MTNSFexplained] = pca(MTNSFresponse);
% reconstructing the responses
MTNSFs=zeros(size(MTNSFscore,1),size(MTNSFcoeff,1));
for u=1:size(MTNSFcoeff,1)
for v=1:size(MTNSFcoeff,2)
MTNSFs(:,u)=MTNSFs(:,u)+MTNSFcoeff(u,v).*MTNSFscore(:,v);
end
end

MTNSFstdscore=std(MTNSFscore,[],1);

MTNSFdata=[MTNSFstdscore(1).*MTNSFcoeff(:,1),MTNSFstdscore(2).*MTNSFcoeff(:,2),MTNSFstdscore(3).*MTNSFcoeff(:,3),MTNSFstdscore(4).*MTNSFcoeff(:,4),MTNSFstdscore(5).*MTNSFcoeff(:,5)];

% MTNSFdata=[MTNSFstdscore(1).*MTNSFcoeff(:,1),MTNSFstdscore(2).*MTNSFcoeff(:,2),MTNSFstdscore(3).*MTNSFcoeff(:,3),MTNSFstdscore(4).*MTNSFcoeff(:,4),MTNSFstdscore(5).*MTNSFcoeff(:,5),...
%     MTNSFstdscore(6).*MTNSFcoeff(:,6),MTNSFstdscore(7).*MTNSFcoeff(:,7),MTNSFstdscore(8).*MTNSFcoeff(:,8),MTNSFstdscore(9).*MTNSFcoeff(:,9),MTNSFstdscore(10).*MTNSFcoeff(:,10),...
%     MTNSFstdscore(11).*MTNSFcoeff(:,11),MTNSFstdscore(12).*MTNSFcoeff(:,12),MTNSFstdscore(13).*MTNSFcoeff(:,13),MTNSFstdscore(14).*MTNSFcoeff(:,14),MTNSFstdscore(15).*MTNSFcoeff(:,15)];

MTNSFidx = cluster(GM,MTNSFdata);
MTNSFpost=posterior(GM,MTNSFdata);        %calculates the posterior probability of cells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;scatter(1:length(MTNSFpost(MTNSFidx==1,1)),MTNSFpost(MTNSFidx==1,1),'b')
figure;scatter(1:length(MTNSFpost(MTNSFidx==2,2)),MTNSFpost(MTNSFidx==2,2),'r')
MTNSFindPost1=intersect(find(MTNSFidx==1),find(MTNSFpost(:,1)<postPcut));
MTNSFindPost2=intersect(find(MTNSFidx==2),find(MTNSFpost(:,2)<postPcut));
MTNSFcombIndPost=[MTNSFindPost1;MTNSFindPost2];
% create MTNSFdataT, MTNSFinddataT, MTNSFidxT, and MTNSFnorpResponseT that include only cells with posterior probabilities higher than postPcut 
% inddarta is the index to be used to recover the original indexing 
MTNSFdataT=MTNSFdata;
MTNSFdataT(MTNSFcombIndPost,:)=[];
MTNSFinddata=(1:size(MTNSFdata,1))';
MTNSFinddataT=MTNSFinddata;
MTNSFinddataT(MTNSFcombIndPost,:)=[];
MTNSFidxT=MTNSFidx;
MTNSFidxT(MTNSFcombIndPost,:)=[];
MTNSFnorpResponseT=MTNSFnorpResponse;
MTNSFnorpResponseT(:,MTNSFcombIndPost)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% cluster MTNSF cells
MTNSFclust1=numel(find(MTNSFidx==1));
MTNSFclust2=numel(find(MTNSFidx==2));
if MTNSFclust2>MTNSFclust1
    MTNSFmisclass=100*(MTNSFclust2/(MTNSFclust1+MTNSFclust2));
else
    MTNSFmisclass=100*(MTNSFclust1/(MTNSFclust1+MTNSFclust2));
end

% plots individual and mean responses of all MTN and SF-OT retro cells DS cells clustered into 2
% groups.
f15=figure;
set(f15,'position',[5 150 950 800]);
subplot(2,2,1);
plot(MTNSFnorpResponseT(:,MTNSFidxT==1),'k');
xlim([0 size(MTNSFnorpResponseT,1)]);
subplot(2,2,2);
plot(MTNSFnorpResponseT(:,MTNSFidxT==2),'k');
xlim([0 size(MTNSFnorpResponseT,1)]);
subplot(2,2,3);
shadedErrorBar(1:size(MTNSFnorpResponseT,1),mean(MTNSFnorpResponseT(:,MTNSFidxT==1),2),std(MTNSFnorpResponseT(:,MTNSFidxT==1),[],2),'k');
xlim([0 size(MTNSFnorpResponseT,1)]); ylim([0 1]);
title(['n = ',num2str(sum(MTNSFidxT==1))]);
subplot(2,2,4);
shadedErrorBar(1:size(MTNSFnorpResponseT,1),mean(MTNSFnorpResponseT(:,MTNSFidxT==2),2),std(MTNSFnorpResponseT(:,MTNSFidxT==2),[],2),'k');
xlim([0 size(MTNSFnorpResponseT,1)]); ylim([0 1]);
title(['n = ',num2str(sum(MTNSFidxT==2))]);
subtitle('MTN and SF retro cells');


%% cluster new data sets - SF retro DS cells (both ON and ONOFF were included)
% extract responses from the following file:
% pooledMapAlphaCorr_Calcium imaging data_40_DSI0.30_1.00_created_May_15_2015_15_45_SFDS
filedir=cd;
namesSF=getFileList(filedir,'pooledMapAlphaCorr_Calcium imaging data_40_DSI0.30_1.00_created_Jul_15_2015_14_29__magCut1.5_unsupervisedSF',0,'anywhere');
SF=load(num2str(namesSF{1}));
% standardize the length of SF responses
% determine the minimum number of elements

% standardize the length of responses
for nu=1:numel(SF.pooledmap.pResponse)
    SFnumPresponseT{nu}=SF.pooledmap.pResponse{nu}(1:min(elements));
end
SFpResponse=reshape(cell2mat(SFnumPresponseT'),min(elements),[]);
% normalize the responses
for p=1:size(SFpResponse,2)
    SFnorpResponse(:,p)=mat2gray(SFpResponse(:,p));
end
% smooth the reponses
for p=1:size(SFpResponse,2)
    SFpResponseS(:,p)=smooth(SFpResponse(:,p),15,'moving');
end
% normalize the smoothed responses
SFpResponseNorS=(SFpResponseS-repmat(mean(SFpResponseS(10:40,:),1),size(SFpResponseS,1),1))./(repmat(max(SFpResponseS,[],1),size(SFpResponseS,1),1)-repmat(mean(SFpResponseS(10:40,:),1),size(SFpResponseS,1),1));

% switch between using the normalized responses or not
switch NORMALRESstate
    case 1
        SFresponse=SFpResponseNorS;
    case 0
        SFresponse=SFpResponseS;
end

% PCA analysis on SF cells
[SFcoeff,SFscore,SFlatent,~,SFexplained] = pca(SFresponse);
% reconstructing the responses
SFs=zeros(size(SFscore,1),size(SFcoeff,1));
for u=1:size(SFcoeff,1)
for v=1:size(SFcoeff,2)
SFs(:,u)=SFs(:,u)+SFcoeff(u,v).*SFscore(:,v);
end
end

SFstdscore=std(SFscore,[],1);

SFdata=[SFstdscore(1).*SFcoeff(:,1),SFstdscore(2).*SFcoeff(:,2),SFstdscore(3).*SFcoeff(:,3),SFstdscore(4).*SFcoeff(:,4),SFstdscore(5).*SFcoeff(:,5)];

% SFdata=[SFstdscore(1).*SFcoeff(:,1),SFstdscore(2).*SFcoeff(:,2),SFstdscore(3).*SFcoeff(:,3),SFstdscore(4).*SFcoeff(:,4),SFstdscore(5).*SFcoeff(:,5),...
%     SFstdscore(6).*SFcoeff(:,6),SFstdscore(7).*SFcoeff(:,7),SFstdscore(8).*SFcoeff(:,8),SFstdscore(9).*SFcoeff(:,9),SFstdscore(10).*SFcoeff(:,10),...
%     SFstdscore(11).*SFcoeff(:,11),SFstdscore(12).*SFcoeff(:,12),SFstdscore(13).*SFcoeff(:,13),SFstdscore(14).*SFcoeff(:,14),SFstdscore(15).*SFcoeff(:,15)];

SFidx = cluster(GM,SFdata);
SFpost=posterior(GM,SFdata);        %calculates the posterior probability of cells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;scatter(1:length(SFpost(SFidx==1,1)),SFpost(SFidx==1,1),'b')
figure;scatter(1:length(SFpost(SFidx==2,2)),SFpost(SFidx==2,2),'r')
SFindPost1=intersect(find(SFidx==1),find(SFpost(:,1)<postPcut));
SFindPost2=intersect(find(SFidx==2),find(SFpost(:,2)<postPcut));
SFcombIndPost=[SFindPost1;SFindPost2];
% create SFdataT, SFinddataT, SFidxT, and SFnorpResponseT that include only cells with posterior probabilities higher than postPcut 
% inddarta is the index to be used to recover the original indexing 
SFdataT=SFdata;
SFdataT(SFcombIndPost,:)=[];
SFinddata=(1:size(SFdata,1))';
SFinddataT=SFinddata;
SFinddataT(SFcombIndPost,:)=[];
SFidxT=SFidx;
SFidxT(SFcombIndPost,:)=[];
SFnorpResponseT=SFnorpResponse;
SFnorpResponseT(:,SFcombIndPost)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% cluster SF cells
SFclust1=numel(find(SFidx==1));
SFclust2=numel(find(SFidx==2));
if SFclust2<SFclust1
    SFmisclass=100*(SFclust2/(SFclust1+SFclust2));
else
    SFmisclass=100*(SFclust1/(SFclust1+SFclust2));
end

% plots individual and mean responses of all SF-OT retro cells DS cells clustered into 2
% groups.
f16=figure;
set(f16,'position',[5 150 950 800]);
subplot(2,2,1);
plot(SFnorpResponseT(:,SFidxT==1),'k');
xlim([0 size(SFnorpResponseT,1)]);
subplot(2,2,2);
plot(SFnorpResponseT(:,SFidxT==2),'k');
xlim([0 size(SFnorpResponseT,1)]);
subplot(2,2,3);
shadedErrorBar(1:size(SFnorpResponseT,1),mean(SFnorpResponseT(:,SFidxT==1),2),std(SFnorpResponseT(:,SFidxT==1),[],2),'k');
xlim([0 size(SFnorpResponseT,1)]); ylim([0 1]);
title(['n = ',num2str(sum(SFidxT==1))]);
subplot(2,2,4);
shadedErrorBar(1:size(SFnorpResponseT,1),mean(SFnorpResponseT(:,SFidxT==2),2),std(SFnorpResponseT(:,SFidxT==2),[],2),'k');
xlim([0 size(SFnorpResponseT,1)]); ylim([0 1]);
title(['n = ',num2str(sum(SFidxT==2))]);
subtitle('SF retro cells');


%% cluster new data sets - MTN retro DS cells (both ON and ONOFF were included)
% extract responses from the following file:
% pooledMapAlphaCorr_Calcium imaging data_40_DSI0.30_1.00_created_May_15_2015_15_45_SFDS
filedir=cd;
namesMTN=getFileList(filedir,'pooledMapAlphaCorr_Calcium imaging data_40_DSI0.30_1.00_created_Jul_15_2015_14_28__magCut1.5_unsupervised_MTN',0,'anywhere');
MTN=load(num2str(namesMTN{1}));
% standardize the length of MTN responses
% determine the minimum number of elements

% standardize the length of responses
for nu=1:numel(MTN.pooledmap.pResponse)
    MTNnumPresponseT{nu}=MTN.pooledmap.pResponse{nu}(1:min(elements));
end
MTNpResponse=reshape(cell2mat(MTNnumPresponseT'),min(elements),[]);
% normalize the responses
for p=1:size(MTNpResponse,2)
    MTNnorpResponse(:,p)=mat2gray(MTNpResponse(:,p));
end
% smooth the reponses
for p=1:size(MTNpResponse,2)
    MTNpResponseS(:,p)=smooth(MTNpResponse(:,p),15,'moving');
end
% normalize the smoothed responses
MTNpResponseNorS=(MTNpResponseS-repmat(mean(MTNpResponseS(10:40,:),1),size(MTNpResponseS,1),1))./(repmat(max(MTNpResponseS,[],1),size(MTNpResponseS,1),1)-repmat(mean(MTNpResponseS(10:40,:),1),size(MTNpResponseS,1),1));

% switch between using the normalized responses or not
switch NORMALRESstate
    case 1
        MTNresponse=MTNpResponseNorS;
    case 0
        MTNresponse=MTNpResponseS;
end

% PCA analysis on MTN cells
[MTNcoeff,MTNscore,MTNlatent,~,MTNexplained] = pca(MTNresponse);
% reconstructing the responses
MTNs=zeros(size(MTNscore,1),size(MTNcoeff,1));
for u=1:size(MTNcoeff,1)
for v=1:size(MTNcoeff,2)
MTNs(:,u)=MTNs(:,u)+MTNcoeff(u,v).*MTNscore(:,v);
end
end

MTNstdscore=std(MTNscore,[],1);

MTNdata=[MTNstdscore(1).*MTNcoeff(:,1),MTNstdscore(2).*MTNcoeff(:,2),MTNstdscore(3).*MTNcoeff(:,3),MTNstdscore(4).*MTNcoeff(:,4),MTNstdscore(5).*MTNcoeff(:,5)];

% MTNdata=[MTNstdscore(1).*MTNcoeff(:,1),MTNstdscore(2).*MTNcoeff(:,2),MTNstdscore(3).*MTNcoeff(:,3),MTNstdscore(4).*MTNcoeff(:,4),MTNstdscore(5).*MTNcoeff(:,5),...
%     MTNstdscore(6).*MTNcoeff(:,6),MTNstdscore(7).*MTNcoeff(:,7),MTNstdscore(8).*MTNcoeff(:,8),MTNstdscore(9).*MTNcoeff(:,9),MTNstdscore(10).*MTNcoeff(:,10),...
%     MTNstdscore(11).*MTNcoeff(:,11),MTNstdscore(12).*MTNcoeff(:,12),MTNstdscore(13).*MTNcoeff(:,13),MTNstdscore(14).*MTNcoeff(:,14),MTNstdscore(15).*MTNcoeff(:,15)];

MTNidx = cluster(GM,MTNdata);
MTNpost=posterior(GM,MTNdata);        %calculates the posterior probability of cells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;scatter(1:length(MTNpost(MTNidx==1,1)),MTNpost(MTNidx==1,1),'b')
figure;scatter(1:length(MTNpost(MTNidx==2,2)),MTNpost(MTNidx==2,2),'r')
MTNindPost1=intersect(find(MTNidx==1),find(MTNpost(:,1)<postPcut));
MTNindPost2=intersect(find(MTNidx==2),find(MTNpost(:,2)<postPcut));
MTNcombIndPost=[MTNindPost1;MTNindPost2];
% create MTNdataT, MTNinddataT, MTNidxT, and MTNnorpResponseT that include only cells with posterior probabilities higher than postPcut 
% inddarta is the index to be used to recover the original indexing 
MTNdataT=MTNdata;
MTNdataT(MTNcombIndPost,:)=[];
MTNinddata=(1:size(MTNdata,1))';
MTNinddataT=MTNinddata;
MTNinddataT(MTNcombIndPost,:)=[];
MTNidxT=MTNidx;
MTNidxT(MTNcombIndPost,:)=[];
MTNnorpResponseT=MTNnorpResponse;
MTNnorpResponseT(:,MTNcombIndPost)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cluster MTN cells
MTNclust1=numel(find(MTNidx==1));
MTNclust2=numel(find(MTNidx==2));
if MTNclust2<MTNclust1
    MTNmisclass=100*(MTNclust2/(MTNclust1+MTNclust2));
else
    MTNmisclass=100*(MTNclust1/(MTNclust1+MTNclust2));
end

% plots individual and mean responses of all MTN retro cells DS cells clustered into 2
% groups.
f17=figure;
set(f17,'position',[5 150 950 800]);
subplot(2,2,1);
plot(MTNnorpResponseT(:,MTNidxT==1),'k');
xlim([0 size(MTNnorpResponseT,1)]);
subplot(2,2,2);
plot(MTNnorpResponseT(:,MTNidxT==2),'k');
xlim([0 size(MTNnorpResponseT,1)]);
subplot(2,2,3);
shadedErrorBar(1:size(MTNnorpResponseT,1),mean(MTNnorpResponseT(:,MTNidxT==1),2),std(MTNnorpResponseT(:,MTNidxT==1),[],2),'k');
xlim([0 size(MTNnorpResponseT,1)]); ylim([0 1]);
title(['n = ',num2str(sum(MTNidxT==1))]);
subplot(2,2,4);
shadedErrorBar(1:size(MTNnorpResponseT,1),mean(MTNnorpResponseT(:,MTNidxT==2),2),std(MTNnorpResponseT(:,MTNidxT==2),[],2),'k');
xlim([0 size(MTNnorpResponseT,1)]); ylim([0 1]);
title(['n = ',num2str(sum(MTNidxT==2))]);
subtitle('MTN retro cells');


%%





















%% use Gaussian mixture models for clusterring
% repeat the evaluation of the number of clusters (every time it uses a different initial value) to determine the optimal number of clusters 
for nevaGM=2:1:20
evaGM=evalclusters(data,'gmdistribution','CalinskiHarabasz','KList',[1:6])
allOptimalK(nevaGM)=evaGM.OptimalK;
allOptimalY{nevaGM}=evaGM.OptimalY;
end
ind=find(allOptimalK==mode(allOptimalK));
finalOptimalK=allOptimalK(ind(1));
%finalOptimalY=allOptimalY{finalOptimalK(1)};

% recalculate the Gaussian mixture models with the optimal number of
% clusters (finalOptimalK) found above to calculate the AIC, BIC, and
% negative log likelihood
options = statset('MaxIter',1000);
for nK=2:1:10
GModel=fitgmdist(data,nK,'RegularizationValue',0.01,'Replicates',10,'Options',options);
idx(:,nK)=cluster(GModel,data);
P{nK}=posterior(GModel,data);
allAIC(nK)=GModel.AIC;
allBIC(nK)=GModel.BIC;
allLL(nK)=GModel.NegativeLogLikelihood;
end

idx = cluster(gm,Y);



evaCH=evalclusters(data,'kmeans','CalinskiHarabasz','KList',[1:6])
evaDB=evalclusters(data,'kmeans','DaviesBouldin','KList',[1:6])
evaG=evalclusters(data,'kmeans','gap','KList',[1:6])
evaS=evalclusters(data,'kmeans','silhouette','KList',[1:6])


evaCHg=evalclusters(data,'linkage','CalinskiHarabasz','KList',[1:6])
tree = linkage(data,'average');
figure()
dendrogram(tree)





idx2=kmeans(data,2,'Display','final','Replicates',5);
figure;
[silh2,h] = silhouette(data,idx2,'cityblock');
h = gca;
h.Children.EdgeColor = [.8 .8 1];
xlabel 'Silhouette Value';
ylabel 'Cluster';

idx3=kmeans(data,3,'Display','final','Replicates',5);
figure;
[silh3,h] = silhouette(data,idx3,'cityblock');
h = gca;
h.Children.EdgeColor = [.8 .8 1];
xlabel 'Silhouette Value';
ylabel 'Cluster';

idx4=kmeans(data,4,'Display','final','Replicates',5);
figure;
[silh4,h] = silhouette(data,idx4,'cityblock');
h = gca;
h.Children.EdgeColor = [.8 .8 1];
xlabel 'Silhouette Value';
ylabel 'Cluster';

idx5=kmeans(data,5,'Display','final','Replicates',5);
figure;
[silh5,h] = silhouette(data,idx5,'cityblock');
h = gca;
h.Children.EdgeColor = [.8 .8 1];
xlabel 'Silhouette Value';
ylabel 'Cluster';

cluster2=mean(silh2);
cluster3=mean(silh3);
cluster4=mean(silh4);
cluster5=mean(silh5);


idx=idx4;
%idx=evaGM.OptimalY;
%idx=evaCHg.OptimalY;

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
plot(mean(DSnorpResponse(:,idx==1),2),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,5,7);
plot(mean(DSnorpResponse(:,idx==2),2),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,5,8);
plot(mean(DSnorpResponse(:,idx==3),2),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,5,9);
plot(mean(DSnorpResponse(:,idx==4),2),'k');
xlim([0 size(DSnorpResponse,1)]);
subplot(2,5,10);
plot(mean(DSnorpResponse(:,idx==5),2),'k');
xlim([0 size(DSnorpResponse,1)]);


sum(idx==1)
sum(idx==2)
sum(idx==3)











%p
%%%%% CLUSTERING %%%%%%
% X=isDSnorpResponse';
% [idx,C] = kmeans(tmp,3);
%
% figure;
% plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',12)
% hold on
% plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12)
% plot(C(:,1),C(:,2),'kx',...
%      'MarkerSize',15,'LineWidth',3)
% legend('Cluster 1','Cluster 2','Centroids',...
%        'Location','NW')
% title 'Cluster Assignments and Centroids'
% hold off

end
