function calculateVecFieldStepwise(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,numIter,polarity2)

load(FileName);     % loads pooledmap file

polarity=1;
discPoints=40;
vizAngle=180;
interpType=2;       % find closest grid point

%% loads all six single-channel maps
  
[FileName,PathName] = uigetfile('*.mat','Select real allPossibleFields map','MultiSelect','on');
real=load(FileName);
      
%% code for generating matchind matrices and plotting contour plots
c=10;
cutoff=(['c',num2str(c)]);
k=1;
h=1;
for alpha=0:resolution:360
    for beta=0:resolution:180 
        mat180real(k,h)=real.allPossibleFields.matchind180(k,h).(cutoff);
        mat90real(k,h)=real.allPossibleFields.matchind90(k,h).(cutoff);        
        k=k+1;
    end
    k=1;
    h=h+1;
end
%% select data type

if DATAMODELstate
   data=pooledmap.clusterdata{1,CLUSTNstate}; 
elseif ~DATAMODELstate
   data=pooledmap;
end

%% find max of mat180real

mat180res=mat180real;       % for the first step, set mat180res to mat180real 

for i=1:numIter
    
[M,I]=max(mat180res(:));
[Irow,Icol]=ind2sub(size(mat180real),I);
alpharange=0:resolution:360;
betarange=0:resolution:180;
alphamax=alpharange(Icol);
betamax=betarange(Irow);

%% calculate the vector field for alphamax and betamax, and place it in maxVec.X, maxVec.Y, and maxVec.Z

    clear allPossibleFields
    allPossibleFields.field={};
    k=1;
    h=1;
    a=deg2rad(360*rand(size(data.XYZvecComp.VecXcompLSC,1),1));     % for compatability with VecComparisonPooled4
    a=a';
    maxVec=[];  % for compatability with VecComparisonPooled4
    histType=5;
alpha=alphamax;
beta=betamax;
fieldType=['alpha_',num2str(alpha),'_','beta_',num2str(beta)];
allPossibleFields=VecComparisonPooled4(data,a,maxVec,deg2rad(alpha),deg2rad(beta),polarity,fieldType,vizAngle,allPossibleFields,k,h,interpType,histType);
fieldType

figure;
quiver3(pooledmap.XYZvecComp.X,pooledmap.XYZvecComp.Y,pooledmap.XYZvecComp.Z,allPossibleFields.VecXcomp{:,1},allPossibleFields.VecYcomp{:,1},allPossibleFields.VecZcomp{:,1});

maxVec.X=allPossibleFields.VecXcomp{:,1};
maxVec.Y=allPossibleFields.VecYcomp{:,1};
maxVec.Z=allPossibleFields.VecZcomp{:,1};

%%

    clear allPossibleFields
    allPossibleFields.field={};
    k=1;
    h=1;
    a=deg2rad(360*rand(size(data.XYZvecComp.VecXcompLSC,1),1));         % random angles between 0 and 360 degrees
    a=a';
    histType=20;
for alpha=0:resolution:360
for beta=0:resolution:180
    fieldType=['alpha_',num2str(alpha),'_','beta_',num2str(beta),'_','iteration_',num2str(i)];
    allPossibleFields=VecComparisonPooled4(data,a,maxVec,deg2rad(alpha),deg2rad(beta),polarity,fieldType,vizAngle,allPossibleFields,k,h,interpType,histType);
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




%% code for generating matchind matrices
c=10;
cutoff=(['c',num2str(c)]);
k=1;
h=1;
for alpha=0:resolution:360
for beta=0:resolution:180
    step.maxmat(i).mat180(k,h)=allPossibleFields.matchind180(k,h).(cutoff); 
    k=k+1;
end
k=1;
h=h+1;
end

mat180max= step.maxmat(i).mat180;      % extract the current map to be used for the fit below

% plot real and max heat maps
f1=figure;
alpha=0:resolution:360;
beta=0:resolution:180;
subplot(2,2,1);
contourf(alpha,beta,mat180real);
subplot(2,2,2);
contourf(alpha,beta,mat180max);


%% fit  regressor and Cr  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C1=reshape(mat180max,[],1);
Cr=reshape(mat180real,[],1);
if i==1
regressor=C1;
else
regressor=[regressor,C1];
end

% non linear optimization
guess=repmat(0.5,1,i);

options=optimset('MaxFunEvals',1e+9,'MaxIter',1000,'Display','iter');
[x180,resnorm180,residual180,exitflag180,output180]=lsqnonlin(@(guess) myfunSurfDynamic(guess,regressor,Cr),guess,0,1,options)  

% multiple linear regression
[b180,bint180,r180,rint180,stats180]=regress(Cr,regressor);
        
% multiple regrression without a constant from my dimentionality of color
% vision work based on Maloney's paper.
X=regressor;
Y=Cr;
B=X\Y;          % this the same as: B=inv(X'*X)*X'*Y;
SSR=B'*X'*Y;
SST=Y'*Y;
R2_180=SSR/SST;

% reconstruct the contour using the estimated coefficients
Cfit=0;
for k=1:size(x180,2)
Cfit=Cfit+x180(k).*regressor(:,k);
end
CfitRs180=reshape(Cfit,size(mat180real,1),[]);

figure;
contourf(alpha,beta,CfitRs180);

mat180res=mat180real-CfitRs180;       % set the now real matrix to the residual matrix after the fit

figure;
contourf(alpha,beta,mat180res);

% [M,~]=min(mat180real(:));
% Mrep=repmat(M,size(mat180real,1),size(mat180real,2));
% mat180real=mat180real-Mrep; 

% figure;
% contourf(alpha,beta,mat180real);        % adjust mat180real to be non-negative



%%
        
step.alphamax(i,1)=alphamax;
step.betamax(i,1)=betamax;
step.R2_180(i,1)=R2_180;
step.b180{i,:}=b180';
step.x180{i,:}=x180;
step.FileName{i,1}=FileName;
step.CfitRs180{i,1}=CfitRs180;
step.mat180res{i,1}=mat180res;
step.regressor{i,1}=regressor;
if i==1
    step.mat180real{i,1}=mat180real;
end
step.maxVec=maxVec;
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
            case 20
                    save(['stepAlphaCorr_',num2str(filename),num2str(discPoints),'_DS_created_',currenttime,'.mat'],'step');
                    %save(['fitRandSummaryAlphaCorr_',num2str(filename),num2str(discPoints),'_DS_created_',currenttime,'.mat'],'fitRandSummary');
        end
        
    case 0
        switch histType  
            case 20
                    save(['step_',num2str(filename),num2str(discPoints),'_DS_created_',currenttime,'.mat'],'step');
                    %save(['fitRandSummary_',num2str(filename),num2str(discPoints),'_DS_created_',currenttime,'.mat'],'fitRandSummary');
        end
end

end

