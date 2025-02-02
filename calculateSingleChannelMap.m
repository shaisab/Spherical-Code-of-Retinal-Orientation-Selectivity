function calculateSingleChannelMap(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,polarity2,addNoise,mapType,fieldName,alphamax,betamax,GSstate)

load(FileName);     % loads pooledmap file

polarity=1;
discPoints=40;
vizAngle=180;
interpType=2;       % find closest grid point

%% select data type

if DATAMODELstate
   data=pooledmap.clusterdata{1,CLUSTNstate}; 
elseif ~DATAMODELstate
   data=pooledmap;
end

%% calculate the vector field for alphamax and betamax, and place it in maxVec.X, maxVec.Y, and maxVec.Z

    clear allPossibleFields
    allPossibleFields.field={};
    k=1;
    h=1;
    a=[];         % for compatability with VecComparisonPooled4
    maxVec.X=[];  % for compatability with VecComparisonPooled4
    maxVec.Y=[];  % for compatability with VecComparisonPooled4
    maxVec.Z=[];  % for compatability with VecComparisonPooled4
    
fieldID=['alpha_',num2str(alphamax),'_','beta_',num2str(betamax)];

switch mapType
    case 1
allPossibleFields=VecComparisonPooled4(data,a,maxVec,deg2rad(alphamax),deg2rad(betamax),polarity,fieldID,vizAngle,allPossibleFields,k,h,interpType,histType);
    case 0
allPossibleFields=VecComparisonPooled4Trans(data,a,maxVec,deg2rad(alphamax),deg2rad(betamax),polarity,fieldID,vizAngle,allPossibleFields,k,h,interpType,histType);
    case 2
allPossibleFields=VecComparisonPooled4Trans(data,a,maxVec,deg2rad(alphamax),deg2rad(betamax),polarity,fieldID,vizAngle,allPossibleFields,k,h,interpType,histType);
end
fieldID

figure;
quiver3(data.XYZvecComp.X,data.XYZvecComp.Y,data.XYZvecComp.Z,allPossibleFields.VecXcomp{:,1},allPossibleFields.VecYcomp{:,1},allPossibleFields.VecZcomp{:,1});

[StU,StV,StVecU,StVecV]=XYZtoStandard(data.XYZvecComp.X,data.XYZvecComp.Y,data.XYZvecComp.Z,allPossibleFields.VecXcomp{:,1},allPossibleFields.VecYcomp{:,1},allPossibleFields.VecZcomp{:,1},discPoints)

quiver(StU,StV,StVecU',StVecV','k');

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
    a=[];
    histType=20;
for alpha=120:resolution:300 
% for alpha=0:resolution:360
for beta=0:resolution:90
% for beta=0:resolution:180
    
%     if  beta==0
%         beta=180;
%     else
%         beta=mod(180-beta,180);
%     end
%     
%     [alpha,beta]=rotateNormalToExtraPersonalSpace3('Right',alpha,beta);
    
    fieldID=['alpha_',num2str(alpha),'_','beta_',num2str(beta),'_','iteration_'];
    
    switch mapType
    case 1
    allPossibleFields=VecComparisonPooled4(data,a,maxVecNoise,deg2rad(alpha),deg2rad(beta),polarity,fieldID,vizAngle,allPossibleFields,k,h,interpType,histType);
    case 0
    allPossibleFields=VecComparisonPooled4Trans(data,a,maxVecNoise,deg2rad(alpha),deg2rad(beta),polarity,fieldID,vizAngle,allPossibleFields,k,h,interpType,histType);
    case 2
    allPossibleFields=VecComparisonPooled4Trans(data,a,maxVecNoise,deg2rad(alpha),deg2rad(beta),polarity,fieldID,vizAngle,allPossibleFields,k,h,interpType,histType);

    end   
    
    k=k+1;
end
k=1;
h=h+1;
fieldID
end

allPossibleFields.metadata.alphamax=alphamax;
allPossibleFields.metadata.betamax=betamax;
allPossibleFields.metadata.mapType=mapType;
allPossibleFields.metadata.ALPHACORRstate=ALPHACORRstate;
allPossibleFields.metadata.DATAMODELstate=DATAMODELstate;
allPossibleFields.metadata.CLUSTNstate=CLUSTNstate;
allPossibleFields.metadata.histType=histType;
allPossibleFields.metadata.FileName=FileName;
allPossibleFields.metadata.resolution=resolution;

allPossibleFields.maxVec=maxVecNoise;

allPossibleFields.St=StU;
allPossibleFields.St=StV;
allPossibleFields.St=StVecU;
allPossibleFields.St=StVecV;

%% code for generating matchind matrices
c=10;
cutoff=(['c',num2str(c)]);
k=1;
h=1;
for alpha=120:resolution:300 
% for alpha=0:resolution:360
for beta=0:resolution:90 
% for beta=0:resolution:180
    single.mat180(k,h)=allPossibleFields.matchind180(k,h).(cutoff);
    single.matMix(k,h)=allPossibleFields.matchindMix(k,h).(cutoff);
    try single.matOSI(k,h)=allPossibleFields.matchindOSI(k,h).(cutoff); end
    k=k+1;
end
k=1;
h=h+1;
end

mat180single=single.mat180;      % extract the current map to be used for the fit below
matMixsingle=single.matMix;      % extract the current map to be used for the fit below
%mat180single=single.matOSI;      % extract the current map to be used for the fit below
%% convert alpha beta map to azimuth elevation map for global space and plot it
if GSstate==1
elements=floor(size(mat180single,2)/2);
d=circshift(mat180single,[0 elements]);   %shift azimuth 0 to middle of graph
d(:,elements+1)=[];   % remove the duplicate column at zero
d=[d(:,end),d];   % add the last column to the start
allPossibleFields.mat180single=d;
f1=figure;
hold on
% alpha=-180:resolution:180;
alpha=-90:resolution:90;
% beta=-90:resolution:90;
beta=-45:resolution:45;
contourf(alpha,beta,d);  
% xlim([-180 180]);ylim([-90 90]);
xlim([-90 90]);ylim([-45 45]);

% a(a>180)=a(a>180)-360;
% a=-a;
% plot([-180 180],[0 0],'m','LineWidth',1);
% plot([0 0],[-90 90],'m','LineWidth',1);
% %c=colorbar('northoutside');
% grid on
% ax = gca;
% ax.XTick = -180:45:180;
% ax.YTick = -90:45:90;
% %ax.XTickLabel = {'N','45','D','135','T','225','V','315','N'};
% ax.FontSize=12;
% xlabel('Azimuth (deg.)','FontSize',16);
% ax.YTickLabel = {'','','','',''};
% 
% ylabh=ylabel('Elevation (deg.)','FontSize',16);
% %set(ylabh, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
% ax.YTickLabel = {'nadir','45', 'horizon','135','zenith'};
% 
% c.Label.String='Goodness of fit (%)';
% c.FontSize=14;
% 
% hold off;

%% plot meanOutput map
dOutput=circshift(allPossibleFields.meanOutput,[0 elements]);   %shift azimuth 0 to middle of graph
dOutput(:,elements+1)=[];
dOutput=[dOutput(:,end),dOutput];
dOutput=fliplr(dOutput);
allPossibleFields.meanOutput=dOutput;
f2=figure;
hold on
% alpha=-180:resolution:180;
alpha=-90:resolution:90;
% beta=-90:resolution:90;
beta=-45:resolution:45;
contourf(alpha,beta,dOutput); 
% xlim([-180 180]);ylim([-90 90]);
xlim([-90 90]);ylim([-45 45]);

% figure
% alpha=-180:resolution:180;
% beta=-90:resolution:90;
% contourf(alpha,beta,allPossibleFields.meanOutput); 
% xlim([-180 180]);ylim([-90 90]);
% 
% hold on
% 
% plot([-180 180],[0 0],'m','LineWidth',1);
% plot([0 0],[-90 90],'m','LineWidth',1);
% c=colorbar('northoutside');
% grid on
% ax = gca;
% ax.XTick = -180:45:180;
% ax.YTick = -90:45:90;
% ax.XTickLabel = {'N','45','D','135','T','225','V','315','N'};
% ax.FontSize=12;
% xlabel('Azimuth (deg.)','FontSize',16);
% ax.YTickLabel = {'','','','',''};
% 
% ylabh=ylabel('Elevation (deg.)','FontSize',16);
% set(ylabh, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
% ax.YTickLabel = {'nadir','45', 'horizon','135','zenith'};
% 
% c.Label.String='Goodness of fit (%)';
% c.FontSize=14;
% 
% hold off;

else
%% plot single map
f1=figure;
alpha=120:resolution:300;
% alpha=0:resolution:360;
beta=0:resolution:90;
% beta=0:resolution:180;
if mapType==2
    contourf(alpha,beta,matMixsingle);
else
    contourf(alpha,beta,mat180single);
end
% hold on;
% scatter(alphamax,betamax,'kx','LineWidth',3);
% hold off;
%% plot meanOutput map
f2=figure;
alpha=120:resolution:300; 
% alpha=0:resolution:360;
beta=0:resolution:90; 
% beta=0:resolution:180;
try contourf(alpha,beta,allPossibleFields.meanOutput); end
end
%% save single map mat file
filedir=cd;
ind=strfind(filedir,'\');
filename=filedir;
filename(ind(end):length(filename))=[];
filename(1:ind(end-1))=[];
currenttime=datestr(now,'mmm_dd_yyyy_HH_MM');
switch mapType
case 0          %translation
switch ALPHACORRstate
    case 1
        save(['allPossibleFieldsAlphaCorr_',num2str(filename),'_translation_','alpha_',num2str(alphamax),'_beta_',num2str(betamax),'_noise_',num2str(addNoise),'_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');
    case 0
        save(['allPossibleFields_',num2str(filename),'_translation_','alpha_',num2str(alphamax),'_beta_',num2str(betamax),'_noise_',num2str(addNoise),'_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');
end

case 2          % mixed orthogonal
switch ALPHACORRstate
    case 1
        save(['allPossibleFieldsAlphaCorr_',num2str(filename),'_mixed orthogonal_','alpha_',num2str(alphamax),'_beta_',num2str(betamax),'_noise_',num2str(addNoise),'_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');
    case 0
        save(['allPossibleFields_',num2str(filename),'_mixed orthogonal_','alpha_',num2str(alphamax),'_beta_',num2str(betamax),'_noise_',num2str(addNoise),'_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');
end

case 1          % rotation
switch ALPHACORRstate
    case 1
        if isempty(fieldName) 
            save(['allPossibleFieldsAlphaCorr_',num2str(filename),'_rotation_','alpha_',num2str(alphamax),'_beta_',num2str(betamax),'_noise_',num2str(addNoise),'_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');
        else
            save(['allPossibleFieldsAlphaCorr_',num2str(filename),'_rotation_',num2str(fieldName),'_noise_',num2str(addNoise),'_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');
        end
    case 0
        save(['allPossibleFields_',num2str(filename),'_rotation_',num2str(fieldName),'_noise_',num2str(addNoise),'_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');
end
end
%% save single map figure
switch mapType
case 0          %translation
switch ALPHACORRstate
    case 1
        savefig(f1,['allPossibleFieldsAlphaCorr_',num2str(filename),'_translation_','alpha_',num2str(alphamax),'_beta_',num2str(betamax),'_noise_',num2str(addNoise),'_cutoff_',num2str(c),'_DS_created_',currenttime]);
        try savefig(f2,['allPossibleFieldsAlphaCorr_',num2str(filename),'_translation_','mean output_','alpha_',num2str(alphamax),'_beta_',num2str(betamax),'_noise_',num2str(addNoise),'_cutoff_',num2str(c),'_DS_created_',currenttime]);end
    case 0
        savefig(f1,['allPossibleFields_',num2str(filename),'_translation_','alpha_',num2str(alphamax),'_beta_',num2str(betamax),'_noise_',num2str(addNoise),'_cutoff_',num2str(c),'_DS_created_',currenttime]);
end

case 2          % mixed orthogonal
switch ALPHACORRstate
    case 1
        savefig(f1,['allPossibleFieldsAlphaCorr_',num2str(filename),'_mixed orthogonal_','alpha_',num2str(alphamax),'_beta_',num2str(betamax),'_noise_',num2str(addNoise),'_cutoff_',num2str(c),'_DS_created_',currenttime]);
        try savefig(f2,['allPossibleFieldsAlphaCorr_',num2str(filename),'_mixed orthogonal_','mean output_','alpha_',num2str(alphamax),'_beta_',num2str(betamax),'_noise_',num2str(addNoise),'_cutoff_',num2str(c),'_DS_created_',currenttime]);end
    case 0
        savefig(f1,['allPossibleFields_',num2str(filename),'_mixed orthogonal_','alpha_',num2str(alphamax),'_beta_',num2str(betamax),'_noise_',num2str(addNoise),'_cutoff_',num2str(c),'_DS_created_',currenttime]);
end


    case 1      % rotation
 switch ALPHACORRstate
    case 1
        if isempty(fieldName) 
        savefig(f1,['allPossibleFieldsAlphaCorr_',num2str(filename),'_rotation_','alpha_',num2str(alphamax),'_beta_',num2str(betamax),'_noise_',num2str(addNoise),'_cutoff_',num2str(c),'_DS_created_',currenttime]);
        try savefig(f2,['allPossibleFieldsAlphaCorr_',num2str(filename),'_translation_','mean output_','alpha_',num2str(alphamax),'_beta_',num2str(betamax),'_noise_',num2str(addNoise),'_cutoff_',num2str(c),'_DS_created_',currenttime]);end
        else
        savefig(f1,['allPossibleFieldsAlphaCorr_',num2str(filename),'_rotation_',num2str(fieldName),'_noise_',num2str(addNoise),'_cutoff_',num2str(c),'_DS_created_',currenttime]);
        try savefig(f2,['allPossibleFieldsAlphaCorr_',num2str(filename),'_translation_','mean output_','alpha_',num2str(alphamax),'_beta_',num2str(betamax),'_noise_',num2str(addNoise),'_cutoff_',num2str(c),'_DS_created_',currenttime]);end
        end
    case 0
        savefig(f1,['allPossibleFields_',num2str(filename),'_rotation_','mean output_',num2str(fieldName),'_noise_',num2str(addNoise),'_cutoff_',num2str(c),'_DS_created_',currenttime]);
 end
end

end

