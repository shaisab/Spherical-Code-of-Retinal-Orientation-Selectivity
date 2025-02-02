function calculateVecFieldStepwiseWithNoise(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,numIter,polarity2,addNoise,fieldType)

load(FileName);     % loads pooledmap file

polarity=1;
discPoints=40;
vizAngle=180;
interpType=2;       % find closest grid point

%% loads real alpha-beta map
  
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
   
switch fieldType
    case 1
        [M,I]=max(mat180res(:));
        [Irow,Icol]=ind2sub(size(mat180real),I);
    case 0
        [M,I]=max(mat180res(10,:));
        Icol=I;
        Irow=floor(size(mat180real,1)/2)+1;     % corresponds to beta=90 degrees, the margin of the retina. equals 10 in the case of 10 degrees resolution.
end


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
fieldID=['alpha_',num2str(alpha),'_','beta_',num2str(beta)];

switch fieldType
    case 1
allPossibleFields=VecComparisonPooled4(data,a,maxVec,deg2rad(alpha),deg2rad(beta),polarity,fieldID,vizAngle,allPossibleFields,k,h,interpType,histType);
    case 0
allPossibleFields=VecComparisonPooled4Trans(data,a,maxVec,deg2rad(alpha),deg2rad(beta),polarity,fieldID,vizAngle,allPossibleFields,k,h,interpType,histType);
end
fieldID

figure;
quiver3(data.XYZvecComp.X,data.XYZvecComp.Y,data.XYZvecComp.Z,allPossibleFields.VecXcomp{:,1},allPossibleFields.VecYcomp{:,1},allPossibleFields.VecZcomp{:,1});

maxVec.X=allPossibleFields.VecXcomp{:,1};
maxVec.Y=allPossibleFields.VecYcomp{:,1};
maxVec.Z=allPossibleFields.VecZcomp{:,1};

%% add noise to the direction of vectors in the field

phys=dlmread('phys');       % get R from the phys file
R=phys(1);

numerator=-data.XYZvecComp.Y.*maxVec.X+data.XYZvecComp.X.*maxVec.Y;
denominator=sqrt(2.*data.XYZvecComp.Z.*R-data.XYZvecComp.Z.^2);

%convert the x,y,z components to angles in degrees between 0 and 360
ang=acosd(numerator./denominator);
%ang=acosd((-data.XYZvecComp.Y.*maxVec.X+data.XYZvecComp.X.*maxVec.Y)./sqrt(data.XYZvecComp.X.^2+data.XYZvecComp.Y.^2));



cond=maxVec.X.*((data.XYZvecComp.X.*(R-data.XYZvecComp.Z))./sqrt(data.XYZvecComp.X.^2+data.XYZvecComp.Y.^2))+...
    maxVec.Y.*((data.XYZvecComp.Y.*(R-data.XYZvecComp.Z))./sqrt(data.XYZvecComp.X.^2+data.XYZvecComp.Y.^2))+...
        maxVec.Z.*sqrt((data.XYZvecComp.X.^2+data.XYZvecComp.Y.^2)./R.^2);

for j=1:size(cond,1)
    if cond(j)>=0
        ang(j)=ang(j);
    elseif cond(j)<0
        ang(j)=360-ang(j);
    end
end

% add noise to the angles
angNoise=addNoise.*randn(size(ang,1),1)+ang;        % this is the angle with noise; addNoise is the standard deviation of noise to add in degrees 
angNoise=deg2rad(angNoise);                         % convert angNoise to back to radians so that I can calculate the vector components

% convert the angle back to the x,y,z components
maxVecNoise.X=(1./(sqrt(data.XYZvecComp.X.^2+data.XYZvecComp.Y.^2))).*((((R-data.XYZvecComp.Z).*sin(angNoise).*data.XYZvecComp.X)./R)-cos(angNoise).*data.XYZvecComp.Y);
maxVecNoise.Y=(1./(sqrt(data.XYZvecComp.X.^2+data.XYZvecComp.Y.^2))).*((((R-data.XYZvecComp.Z).*sin(angNoise).*data.XYZvecComp.Y)./R)+cos(angNoise).*data.XYZvecComp.X);
maxVecNoise.Z=(sqrt(2.*R.*data.XYZvecComp.Z-data.XYZvecComp.Z.^2).*sin(angNoise))./R;
       
% normalize Vector Fields
temp1=(maxVecNoise.X.^2+maxVecNoise.Y.^2+maxVecNoise.Z.^2).^(1/2);
maxVecNoise.X=maxVecNoise.X./temp1;
maxVecNoise.Y=maxVecNoise.Y./temp1;
maxVecNoise.Z=maxVecNoise.Z./temp1;

%% generate alpha-beta map with the observed field set to the max field with noise

    clear allPossibleFields
    allPossibleFields.field={};
    k=1;
    h=1;
    a=deg2rad(360*rand(size(data.XYZvecComp.VecXcompLSC,1),1));         % random angles between 0 and 360 degrees
    a=a';
    histType=20;
for alpha=0:resolution:360
for beta=0:resolution:180
    fieldID=['alpha_',num2str(alpha),'_','beta_',num2str(beta),'_','iteration_',num2str(i)];
    
    switch fieldType
    case 1
    allPossibleFields=VecComparisonPooled4(data,a,maxVecNoise,deg2rad(alpha),deg2rad(beta),polarity,fieldID,vizAngle,allPossibleFields,k,h,interpType,histType);
    case 0
        allPossibleFields=VecComparisonPooled4Trans(data,a,maxVecNoise,deg2rad(alpha),deg2rad(beta),polarity,fieldID,vizAngle,allPossibleFields,k,h,interpType,histType);
    end   
    
    k=k+1;
end
k=1;
h=h+1;
fieldID
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

mat180max=100*mat2gray(mat180max);          % normalize mat180max to range between 0 and 100


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

% adjust the size of regressor every step
if i==1
regressor=C1;
else
regressor=[regressor,C1];
end

% adjust the size of lowerBound and upperBound every step
% if i==1
% upperBound=M/100;
% lowerBound=0;
% else
% upperBound=[upperBound,M/100];
% lowerBound=[lowerBound,0];
% end


% non linear optimization
% adjust the size of guess, upperBound, and lowerBound every step
guess=repmat(0.5,1,i);
upperBound=repmat(1,1,i);
lowerBound=repmat(0,1,i);

options=optimset('MaxFunEvals',1e+9,'MaxIter',1000,'Display','iter');
[x180,resnorm180,residual180,exitflag180,output180]=lsqnonlin(@(guess) myfunSurfDynamic(guess,regressor,Cr),guess,lowerBound,upperBound,options)  

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
step.maxVecNoise=maxVecNoise;
%fitRandSummary{i}=ws2struct();

end         % loop over all iterations

%% plot real map with locations of alphamax and betamax

f10=figure;
alpha=0:resolution:360;
beta=0:resolution:180;
contourf(alpha,beta,step.mat180real{1,1});
hold on
stepnum=(1:size(step.alphamax,1))';
scatter(step.alphamax,step.betamax,'r');
text(step.alphamax,step.betamax,[repmat('   ',size(step.alphamax,1),1),num2str(stepnum)],'HorizontalAlignment','left','Color','r');

%% save all data files and summary real map with locations of alphamax and betamax

filedir=cd;
ind=strfind(filedir,'\');
filename=filedir;
filename(ind(end):length(filename))=[];
filename(1:ind(end-1))=[];

currenttime=datestr(now,'mmm_dd_yyyy_HH_MM');

switch ALPHACORRstate
    case 1
        switch fieldType
            case 1
                    save(['stepAlphaCorr_',num2str(filename),'_rotation_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'step'); 
                    savefig(f10,['stepAlphaCorr_',num2str(filename),'_rotation_',num2str(discPoints),'_DS_created_',currenttime]);
            case 0
                    save(['stepAlphaCorr_',num2str(filename),'_translation_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'step');
                    savefig(f10,['stepAlphaCorr_',num2str(filename),'_translation_',num2str(discPoints),'_DS_created_',currenttime]);
        end
        
    case 0
        switch fieldType 
            case 1
                    save(['step_',num2str(filename),'_rotation_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'step');
                    savefig(f10,['step_',num2str(filename),'_rotation_',num2str(discPoints),'_DS_created_',currenttime]);
            case 0
                    save(['step_',num2str(filename),'_translation_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'step'); 
                    savefig(f10,['step_',num2str(filename),'_translation_',num2str(discPoints),'_DS_created_',currenttime]);
        end
end

end

