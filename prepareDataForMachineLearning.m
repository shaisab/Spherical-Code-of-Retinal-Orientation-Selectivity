
function []= prepareDataForMachineLearning()

%%% working directory should be 'All identified cells'


NORMALRESstate=1;

filedir=cd;
physFileList=getFileList(filedir,['physdata_DB'],0,'anywhere');         % used for plotting the onset and offset of stimulus

[filedir,go]=bbdirselector('select data folder',cd);
if ~go disp('no folder selected'); return; end

names=getFileList(filedir,'identifiedCell',0,'anywhere');

TimeONDS=zeros(300,300);
pResponseONDS=zeros(300,300);
TimeONOFFDS=zeros(300,300);
pResponseONOFFDS=zeros(300,300);

if ~isempty(names)
    for ncell=1:size(names,2)
        curname=cell2mat(names(ncell));
        inds=find(curname=='_');
        if strcmp(curname(inds(end)+1:length(curname)-4),'ONDS')
            ONDSres=load(num2str(names{ncell}));
            
            dif=size(TimeONDS,1)-size(ONDSres.Time,1);
            if dif>=0
                TimeONDS(size(ONDSres.Time,1)+1:end,:)=[];
                pResponseONDS(size(ONDSres.pResponse,1)+1:end,:)=[]; 
            elseif dif<0
                ONDSres.Time(size(TimeONDS,1)+1:end,:)=[];
                ONDSres.pResponse(size(pResponseONDS,1)+1:end,:)=[];
            end
            
            TimeONDS(:,ncell)=ONDSres.Time;
            pResponseONDS(:,ncell)=ONDSres.pResponse;   
            IDONDS{ncell}=curname;
            
        elseif strcmp(curname(inds(end)+1:length(curname)-4),'ONOFFDS')
            ONOFFDSres=load(num2str(names{ncell}));
            
            dif=size(TimeONOFFDS,1)-size(ONOFFDSres.Time,1);
            if dif>=0
                TimeONOFFDS(size(ONOFFDSres.Time,1)+1:end,:)=[];
                pResponseONOFFDS(size(ONOFFDSres.pResponse,1)+1:end,:)=[];
            elseif dif<0
                ONOFFDSres.Time(size(TimeONOFFDS,1)+1:end,:)=[];
                ONOFFDSres.pResponse(size(pResponseONOFFDS,1)+1:end,:)=[];
            end
            
            TimeONOFFDS(:,ncell)=ONOFFDSres.Time;
            pResponseONOFFDS(:,ncell)=ONOFFDSres.pResponse;  
            IDONOFFDS{ncell}=curname;
        end
        
    end
end

TimeONDS(:,~any(TimeONDS,1))=[];  %columns
pResponseONDS(:,~any(pResponseONDS,1))=[];  %columns
IDONDS=IDONDS(~cellfun('isempty',IDONDS));  
TimeONOFFDS(:,~any(TimeONOFFDS,1))=[];  %columns
pResponseONOFFDS(:,~any(pResponseONOFFDS,1))=[];  %columns
IDONOFFDS=IDONOFFDS(~cellfun('isempty',IDONOFFDS));  


response=[repmat(1,1,size(pResponseONDS,2)),repmat(2,1,size(pResponseONOFFDS,2))];
predictor=[pResponseONDS,pResponseONOFFDS];


% normalize the predictors
for p=1:size(predictor,2)
    norPredictor(:,p)=mat2gray(predictor(:,p));
end

% smooth the reponses
for p=1:size(predictor,2)
    predictorS(:,p)=smooth(predictor(:,p),15,'moving');
end

% normalize the smoothed predictors
predictorNorS=(predictorS-repmat(mean(predictorS(10:40,:),1),size(predictorS,1),1))./(repmat(max(predictorS,[],1),size(predictorS,1),1)-repmat(mean(predictorS(10:40,:),1),size(predictorS,1),1));

% switch between using the normalized predictors or not
switch NORMALRESstate
    case 1
        DSpredictor=predictorNorS;
    case 0
        DSpredictor=predictorS;
end


% PCA analysis
[coeff,score,latent,~,explained] = pca(DSpredictor);
% reconstructing the responses
s=zeros(size(score,1),size(coeff,1));
for u=1:size(coeff,1)
for v=1:size(coeff,2)
s(:,u)=s(:,u)+coeff(u,v).*score(:,v);
end
end

stdscore=std(score,[],1);


%data=[coeff(:,1),coeff(:,2),coeff(:,3),coeff(:,4),coeff(:,5)];
%data=[coeff(:,1),coeff(:,2),coeff(:,3),coeff(:,4),coeff(:,5),coeff(:,6),coeff(:,7),coeff(:,8),coeff(:,9),coeff(:,10),coeff(:,11),coeff(:,12),coeff(:,13),coeff(:,14),coeff(:,15)];
%data=[stdscore(1).*coeff(:,1),stdscore(2).*coeff(:,2),stdscore(3).*coeff(:,3),stdscore(4).*coeff(:,4),stdscore(5).*coeff(:,5)];

data=table(stdscore(1).*coeff(:,1),stdscore(2).*coeff(:,2),stdscore(3).*coeff(:,3),stdscore(4).*coeff(:,4),stdscore(5).*coeff(:,5),stdscore(6).*coeff(:,6),...
    stdscore(7).*coeff(:,7),stdscore(8).*coeff(:,8),stdscore(9).*coeff(:,9),stdscore(10).*coeff(:,10),stdscore(11).*coeff(:,11),...
    stdscore(12).*coeff(:,12),stdscore(13).*coeff(:,13),stdscore(14).*coeff(:,14),stdscore(15).*coeff(:,15),response');

save('trainset.mat','data');

end