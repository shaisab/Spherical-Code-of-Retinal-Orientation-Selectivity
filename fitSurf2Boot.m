function []= fitSurf2Boot(channelID,resolution,mapType)

warning('off','all');

%% import data
% channelID is equivalent to polarity2

% the single maps and pooledmap file should be located in a folder other
% than the boot folder

% import rotation single channel maps
[FileName,PathName] = uigetfile('*.mat','Select allPossibleFields files of single-channel rotation maps','MultiSelect','on');
% sort file order by ascending alpha
for k=1:size(FileName,2)
    indUnder=find(FileName{k}=='_');
    order(k)=str2num(FileName{k}(indUnder(6)+1:indUnder(7)-1));     % for individual retinas
    %order(k)=str2num(FileName{k}(indUnder(4)+1:indUnder(5)-1));     % for all retinas pooled
    %order(k)=str2num(FileName{k}(indUnder(3)+1:indUnder(4)-1));     % use when the text '_Calcium imaging data' was omitted from file name
end
[~,I]=sort(order);
FileNameSorted=FileName(I);

for i=1:size(FileName,2)
    load(FileNameSorted{i});
    rotField{i}=generateMatchindMatrix(allPossibleFields);
end

% import translation single channel maps
[FileName,PathName] = uigetfile('*.mat','Select allPossibleFields files of single-channel translation or rotation maps','MultiSelect','on');
% sort file order by ascending alpha
for k=1:size(FileName,2)
    indUnder=find(FileName{k}=='_');
    order(k)=str2num(FileName{k}(indUnder(6)+1:indUnder(7)-1));     % for individual retinas
    %order(k)=str2num(FileName{k}(indUnder(4)+1:indUnder(5)-1));     % for all retinas pooled
end
[~,I]=sort(order);
FileNameSorted=FileName(I);

for i=1:size(FileName,2)
    load(FileNameSorted{i});
    transField{i}=generateMatchindMatrix(allPossibleFields);
end

% import real translation and rotation files
filedir=cd;
names=getFileList(filedir,'.mat',0,'anywhere');     % detect all .mat files
Right=[];
type={'translation','rotation'};
for t=1:2
    typestr=type{t};
    for i=1:size(names,2)
        if strfind(names{i},typestr)
            if strfind(names{i},'_nboot')
                load(num2str(names{i}));
                real.(typestr){i}=generateMatchindMatrix(allPossibleFields);
            end
        end
    end
end

real.translation=real.translation(find(~cellfun(@isempty,real.translation)));
real.rotation=real.rotation(find(~cellfun(@isempty,real.rotation)));

save('tmp.mat');

%% perform fit
% four rotation channels based on real singularities

        guess=[0.5,0.5];
        lowerBound=[0,0];
        upperBound=[1,1];
        
  for nboot=1:size(real.rotation,2)      
        C1=reshape(rotField{1},[],1);
        C2=reshape(rotField{2},[],1);
        Cr=reshape(real.rotation{1,nboot},[],1);
        
        regressor=[C1,C2];
        
        options=optimset('MaxFunEvals',1e+9,'MaxIter',1000,'Display','off');
        [x,resnorm,residual,exitflag,output]=lsqnonlin(@(guess)myfunSurfDynamic(guess,regressor,Cr),guess,lowerBound,upperBound,options);
        
        % multiple linear regression
        [b,bint,r,rint,stats]=regress(Cr,[regressor]);
        
        % multiple regrression without a constant from my dimentionality of color
        % vision work based on Maloney's paper.
        X=regressor;
        Y=Cr;
        B=X\Y;          % this the same as: B=inv(X'*X)*X'*Y;
        SSR=B'*X'*Y;
        SST=Y'*Y;
        R2rot(nboot)=SSR/SST;
        R2rot_adjusted=1-((1-R2rot(nboot)).*((length(Cr)-1)/(length(Cr)-2)));
  end       
        
      guess=[0.5,0.5];
        lowerBound=[0,0];
        upperBound=[1,1];
        
  for nboot=1:size(real.translation,2)      
        C1=reshape(transField{1},[],1);
        C2=reshape(transField{2},[],1);
        Cr=reshape(real.translation{1,nboot},[],1);
        
        regressor=[C1,C2];
        
        options=optimset('MaxFunEvals',1e+9,'MaxIter',1000,'Display','off');
        [x,resnorm,residual,exitflag,output]=lsqnonlin(@(guess)myfunSurfDynamic(guess,regressor,Cr),guess,lowerBound,upperBound,options);
        
        % multiple linear regression
        [b,bint,r,rint,stats]=regress(Cr,[regressor]);
        
        % multiple regrression without a constant from my dimentionality of color
        % vision work based on Maloney's paper.
        X=regressor;
        Y=Cr;
        B=X\Y;          % this the same as: B=inv(X'*X)*X'*Y;
        SSR=B'*X'*Y;
        SST=Y'*Y;
        R2trans(nboot)=SSR/SST;
        R2trans_adjusted=1-((1-R2trans(nboot)).*((length(Cr)-1)/(length(Cr)-2)));
        
  end       

%% caclulate means and confidence intervals

% ONDS
% realR2Trans=0.95;
% realR2Rot=0.87;

%Hb9
% realR2Trans=0.78;
% realR2Rot=0.22;

% %%%%%% insert R2 for ON sus
% realR2Trans=0.84;
% realR2Rot=0.88;

% %%%%%% insert R2 for ON Trans
% realR2Rot=0.82;
% realR2Trans=0.88;
% 
% %%%%%% insert R2 for ON OFF Sustained
% realR2Rot=0.85;
% realR2Trans=0.87;

%%%%%% insert R2 for ON OFF Trans
realR2Rot=0.89;
realR2Trans=0.92;

meanRot=mean(R2rot_adjusted)
meanTrans=mean(R2trans_adjusted)

diffRot=realR2Rot-meanRot;
diffTrans=realR2Trans-meanTrans;

R2rot=R2rot_adjusted+repmat(diffRot,1,size(R2rot,2));
R2trans=R2trans_adjusted+repmat(diffTrans,1,size(R2trans,2));

ciRot=quantile(R2rot,[0.025 0.975])
ciTrans=quantile(R2trans,[0.025 0.975])

%% save file and figure
filedir=cd;
ind=strfind(filedir,'\');
filename=filedir;
filename(ind(end):length(filename))=[];
filename(1:ind(end-1))=[];

currenttime=datestr(now,'mmm_dd_yyyy_HH_MM');
save(['fitBoot_',num2str(filename),'_DS_created_',currenttime]);

end

