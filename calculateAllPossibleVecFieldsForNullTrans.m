function calculateAllPossibleVecFieldsForNullTrans(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,numIter,polarity2)

load(FileName);     % loads pooledmap file

polarity=1;
discPoints=40;
vizAngle=180;
interpType=2;       % find closest grid point

%% loads all six single-channel maps
switch polarity2
    case 1
        [FileName,PathName] = uigetfile('*.mat','Select six allPossibleFields single-channel maps','MultiSelect','on');
        for i=1:size(FileName,2)
            if strfind(FileName{i},'_ASC')
                ASC=load(FileName{i});
            elseif strfind(FileName{i},'_PSC')
                PSC=load(FileName{i});
            elseif strfind(FileName{i},'_LSC')
                LSC=load(FileName{i});
            else
                real=load(FileName{i});
            end
        end
        
    case 2
        [FileName,PathName] = uigetfile('*.mat','Select six allPossibleFields single-channel maps','MultiSelect','on');
        for i=1:size(FileName,2)
            if strfind(FileName{i},'_ASC')
                ASC=load(FileName{i});
            elseif strfind(FileName{i},'_PSC')
                PSC=load(FileName{i});
            elseif strfind(FileName{i},'_LSC')
                LSC=load(FileName{i});
            elseif strfind(FileName{i},'_inverseASC')
                inverseASC=load(FileName{i});
            elseif strfind(FileName{i},'_inversePSC')
                inversePSC=load(FileName{i});
            elseif strfind(FileName{i},'_inverseLSC')
                inverseLSC=load(FileName{i});
            elseif strfind(FileName{i},'_real')
                real=load(FileName{i});
            end
        end
end

%% code for generating matchind matrices and plotting contour plots
c=10;
cutoff=(['c',num2str(c)]);
k=1;
h=1;
for alpha=0:resolution:360
    for beta=0:resolution:180
        mat180ASC(k,h)=ASC.allPossibleFields.matchind180(k,h).(cutoff);
        mat90ASC(k,h)=ASC.allPossibleFields.matchind90(k,h).(cutoff);
        mat180PSC(k,h)=PSC.allPossibleFields.matchind180(k,h).(cutoff);
        mat90PSC(k,h)=PSC.allPossibleFields.matchind90(k,h).(cutoff);
        mat180LSC(k,h)=LSC.allPossibleFields.matchind180(k,h).(cutoff);
        mat90LSC(k,h)=LSC.allPossibleFields.matchind90(k,h).(cutoff);
%         mat180real(k,h)=real.allPossibleFields.matchind180(k,h).(cutoff);
%         mat90real(k,h)=real.allPossibleFields.matchind90(k,h).(cutoff);
        try mat180iASC(k,h)=inverseASC.allPossibleFields.matchind180(k,h).(cutoff); end
        try mat90iASC(k,h)=inverseASC.allPossibleFields.matchind90(k,h).(cutoff); end
        try mat180iPSC(k,h)=inversePSC.allPossibleFields.matchind180(k,h).(cutoff); end
        try mat90iPSC(k,h)=inversePSC.allPossibleFields.matchind90(k,h).(cutoff); end
        try mat180iLSC(k,h)=inverseLSC.allPossibleFields.matchind180(k,h).(cutoff); end
        try mat90iLSC(k,h)=inverseLSC.allPossibleFields.matchind90(k,h).(cutoff); end
        k=k+1;
    end
    k=1;
    h=h+1;
end
%% generate alpha-beta maps for random vectors 

if DATAMODELstate
   data=pooledmap.clusterdata{1,CLUSTNstate}; 
elseif ~DATAMODELstate
   data=pooledmap;
end


for i=1:numIter
    clear allPossibleFields
    allPossibleFields.field={};
    k=1;
    h=1;
    a=deg2rad(360*rand(size(data.XYZvecComp.VecXcompLSC,1),1));
    a=a';
for alpha=0:resolution:360
for beta=0:resolution:180
    fieldType=['alpha_',num2str(alpha),'_','beta_',num2str(beta),'_','iteration_',num2str(i)];
    allPossibleFields=VecComparisonPooled4(data,a,deg2rad(alpha),deg2rad(beta),polarity,fieldType,vizAngle,allPossibleFields,k,h,interpType,histType);
    k=k+1;
end
k=1;
h=h+1;
fieldType
end


% currenttime=datestr(now,'mmm_dd_yyyy_HH_MM');
% switch ALPHACORRstate
%     case 1
%         switch histType
%             case 10
%                     save(['allPossibleFieldsAlphaCorr_',num2str(filename),'_rand_',num2str(i),num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');
%         end
%         
%     case 0
%         switch histType  
%             case 10
%                     save(['allPossibleFields_',num2str(filename),'_rand_',num2str(i),num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');
%         end
% end



%%% code for generating matchind matrices
c=10;
cutoff=(['c',num2str(c)]);
k=1;
h=1;
for alpha=0:resolution:360
for beta=0:resolution:180
    null.randmat(i).mat180(k,h)=allPossibleFields.matchind180(k,h).(cutoff);
    %mat90(k,h)=allPossibleFields.matchind90(k,h).(cutoff);
    k=k+1;
end
k=1;
h=h+1;
end

mat180rand= null.randmat(i).mat180;      % extract the current map to be used for the fit below

%% process mat180 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C1=reshape(mat180ASC,[],1);
C2=reshape(mat180PSC,[],1);
C3=reshape(mat180LSC,[],1);
Cr=reshape(mat180rand,[],1);
try Ci1=reshape(mat180iASC,[],1); end
try Ci2=reshape(mat180iPSC,[],1); end
try Ci3=reshape(mat180iLSC,[],1); end

switch polarity2
    case 1
% non linear optimization
guess=[0.5,0.5,0.5];
options=optimset('MaxFunEvals',1e+9,'MaxIter',1000,'Display','iter');
[x180,resnorm180,residual180,exitflag180,output180]=lsqnonlin(@(guess) myfunSurf(guess,C1,C2,C3,Cr),guess,[0,0,0],[1,1,1],options)

% multiple linear regression
[b180,bint180,r180,rint180,stats180]=regress(Cr,[C1,C2,C3]);
        
% multiple regrression without a constant from my dimentionality of color
% vision work based on Maloney's paper.
X=[C1,C2,C3];
Y=Cr;
B=X\Y;          % this the same as: B=inv(X'*X)*X'*Y;
SSR=B'*X'*Y;
SST=Y'*Y;
R2_180=SSR/SST;

% reconstruct the contour using the estimated coefficients
Cfit=x180(1).*C1+x180(2).*C2+x180(3).*C3;
%Cfit=b180(1).*C1+b180(2).*C2+b180(3).*C3;
CfitRs180=reshape(Cfit,size(mat180ASC,1),[]);

% normalize canal coefficients
xn180=x180./sum(x180);

    case 2
        
        % non linear optimization
guess=[0.5,0.5,0.5,0.5,0.5,0.5];
options=optimset('MaxFunEvals',1e+9,'MaxIter',1000,'Display','iter');
[x180,resnorm180,residual180,exitflag180,output180]=lsqnonlin(@(guess) myfunSurf6c(guess,C1,C2,C3,Ci1,Ci2,Ci3,Cr),guess,[0,0,0,0,0,0],[1,1,1,1,1,1],options)

% multiple linear regression
[b180,bint180,r180,rint180,stats180]=regress(Cr,[C1,C2,C3,Ci1,Ci2,Ci3]);
        
% multiple regrression without a constant from my dimentionality of color
% vision work based on Maloney's paper.
X=[C1,C2,C3,Ci1,Ci2,Ci3];
Y=Cr;
B=X\Y;          % this the same as: B=inv(X'*X)*X'*Y;
SSR=B'*X'*Y;
SST=Y'*Y;
R2_180=SSR/SST;

% reconstruct the contour using the estimated coefficients
Cfit=x180(1).*C1+x180(2).*C2+x180(3).*C3+x180(4).*Ci1+x180(5).*Ci2+x180(6).*Ci3;
%Cfit=b180(1).*C1+b180(2).*C2+b180(3).*C3;
CfitRs180=reshape(Cfit,size(mat180ASC,1),[]);

% normalize canal coefficients
xn180=x180./sum(x180);

end
        
null.coeff(i,:)=xn180;
null.R2_180(i,1)=R2_180;
null.b180(i,:)=b180';
null.x180(i,:)=x180;
null.FileName{i,1}=FileName;
%fitRandSummary{i}=ws2struct();

end         % loop over all iterations

%% save all data files

filedir=cd;
ind=strfind(filedir,'\');
filename=filedir;
filename(ind(end):length(filename))=[];
filename(1:ind(end-1))=[];

currenttime=datestr(now,'mmm_dd_yyyy_HH_MM');
switch ALPHACORRstate
    case 1
        switch histType
            case 10
                    save(['nullAlphaCorr_',num2str(filename),num2str(discPoints),'_DS_created_',currenttime,'.mat'],'null');
                    %save(['fitRandSummaryAlphaCorr_',num2str(filename),num2str(discPoints),'_DS_created_',currenttime,'.mat'],'fitRandSummary');
        end
        
    case 0
        switch histType  
            case 10
                    save(['null_',num2str(filename),num2str(discPoints),'_DS_created_',currenttime,'.mat'],'null');
                    %save(['fitRandSummary_',num2str(filename),num2str(discPoints),'_DS_created_',currenttime,'.mat'],'fitRandSummary');
        end
end

end

