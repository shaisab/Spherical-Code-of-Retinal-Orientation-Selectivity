function calculateAllPossibleVecFields(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,fieldType,DISTstate,locDist,widthDist,ABstate,GSstate)

%% load data
load(FileName);     % loads pooledmap file OR combinedMap file OR single-channel map


if DATAMODELstate==1
   data=pooledmap.clusterdata{1,CLUSTNstate}; 
elseif DATAMODELstate==0
   data=pooledmap;
elseif DATAMODELstate==2
   data=combinedMap;
elseif DATAMODELstate==3
   data=allPossibleFields;
end

clear allPossibleFields

%% generate alpha-beta map
polarity=1;
discPoints=40;
vizAngle=180;
interpType=2;       % find closest grid point
allPossibleFields.field={};
k=1;
h=1;
a=[];       % for compatability with VecComparisonPooled4
maxVec=[];  % for compatability with VecComparisonPooled4

%figure;
%for alpha=0:resolution:360 
for alpha=120:resolution:300
%for beta=0:resolution:180
for beta=0:resolution:90
    fieldID=['alpha_',num2str(alpha),'_','beta_',num2str(beta)];
    
    switch fieldType
        case 1  % rotation
            allPossibleFields=VecComparisonPooled4(data,a,maxVec,deg2rad(alpha),deg2rad(beta),polarity,fieldID,vizAngle,allPossibleFields,k,h,interpType,histType,DISTstate,locDist,widthDist);
        case 0  % translation
            allPossibleFields=VecComparisonPooled4Trans(data,a,maxVec,deg2rad(alpha),deg2rad(beta),polarity,fieldID,vizAngle,allPossibleFields,k,h,interpType,histType);
        case 2  % mixed orthoginal (rotation, translation)
            allPossibleFields=VecComparisonPooled4Trans(data,a,maxVec,deg2rad(alpha),deg2rad(beta),polarity,fieldID,vizAngle,allPossibleFields,k,h,interpType,histType);     
    end
    k=k+1;
end
k=1;
h=h+1;
fieldID
end


%% generate matchind matrix
c=10;
cutoff=(['c',num2str(c)]);
k=1;
h=1;
%for alpha=0:resolution:360 
for alpha=120:resolution:300
%for beta=0:resolution:180 
for beta=0:resolution:90
    mat180(k,h)=allPossibleFields.matchind180(k,h).(cutoff);
    matMix(k,h)=allPossibleFields.matchindMix(k,h).(cutoff);
    try matOSI(k,h)=allPossibleFields.matchindOSI(k,h).(cutoff); end
    try mat90(k,h)=allPossibleFields.matchind90(k,h).(cutoff); end
    k=k+1;
end
k=1;
h=h+1;
end

%% find max singularities in alpha beta real translation map
if fieldType==0
mat180tmp=mat180;
% mat180tmp=mat90;
% mat180tmp=matOSI;

% UNCOMMENT FOR ALL RETINAS POOLED
for i=1:2
    if i==1
        [~,I]=max(mat180tmp(:));
        [Irow,Icol]=ind2sub(size(mat180tmp),I);
    else
        if i==1
            [~,I]=max(mat180tmp(:));
            [Irow,Icol]=ind2sub(size(mat180tmp),I);
        else
            if Icol<=(70/resolution)-1 && Icol>0
                mat180tmp(:,Icol-Icol+1:Icol+(70/resolution))=0;
                mat180tmp(:,end-((70/resolution)-Icol-1):end)=0;
                [~,I]=max(mat180tmp(:));
                [Irow,Icol]=ind2sub(size(mat180tmp),I);
            elseif Icol>=size(mat180tmp,2)-(70/resolution)+2
                mat180tmp(:,Icol-(70/resolution):Icol)=0;
                [~,I]=max(mat180tmp(:));
                [Irow,Icol]=ind2sub(size(mat180tmp),I);
            elseif Icol>(70/resolution) && Icol<size(mat180tmp,2)-(70/resolution)+2
                mat180tmp(:,Icol-(70/resolution):Icol+(70/resolution))=0;
                [~,I]=max(mat180tmp(:));
                [Irow,Icol]=ind2sub(size(mat180tmp),I);
            end
        end
    end
    Icolmax(i)=Icol;
    Irowmax(i)=Irow;
end

% fix for individual retinas
% Icolmax=[1,2,3,4];
% Irowmax=[1,2,3,4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Icolmax,J]=sort(Icolmax);      % sort Icolmax to assend along alpha
Irowmax=Irowmax(J);             % sort Irowmax to assend along alpha

alpha=120:resolution:300;
%alpha=0:resolution:360;
beta=0:resolution:90;
%beta=0:resolution:180;
allPossibleFields.alphasing=alpha(Icolmax);
allPossibleFields.betasing=beta(Irowmax);
allPossibleFields.Icolmax=Icolmax;
allPossibleFields.Irowmax=Irowmax;
end

%% find max singularities in alpha beta real mixed orthogonal map
if fieldType==2
mat180tmp=matMix;

% UNCOMMENT FOR ALL RETINAS POOLED
for i=1:1
    if i==1
        [~,I]=max(mat180tmp(:));
        [Irow,Icol]=ind2sub(size(mat180tmp),I);
    else
        if i==1
            [~,I]=max(mat180tmp(:));
            [Irow,Icol]=ind2sub(size(mat180tmp),I);
        else
            if Icol<=(70/resolution)-1 && Icol>0
                mat180tmp(:,Icol-Icol+1:Icol+(70/resolution))=0;
                mat180tmp(:,end-((70/resolution)-Icol-1):end)=0;
                [~,I]=max(mat180tmp(:));
                [Irow,Icol]=ind2sub(size(mat180tmp),I);
            elseif Icol>=size(mat180tmp,2)-(70/resolution)+2
                mat180tmp(:,Icol-(70/resolution):Icol)=0;
                [~,I]=max(mat180tmp(:));
                [Irow,Icol]=ind2sub(size(mat180tmp),I);
            elseif Icol>(70/resolution) && Icol<size(mat180tmp,2)-(70/resolution)+2
                mat180tmp(:,Icol-(70/resolution):Icol+(70/resolution))=0;
                [~,I]=max(mat180tmp(:));
                [Irow,Icol]=ind2sub(size(mat180tmp),I);
            end
        end
    end
    Icolmax(i)=Icol;
    Irowmax(i)=Irow;
end

% fix for individual retinas
% Icolmax=[1,2,3,4];
% Irowmax=[1,2,3,4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Icolmax,J]=sort(Icolmax);      % sort Icolmax to assend along alpha
Irowmax=Irowmax(J);             % sort Irowmax to assend along alpha

alpha=120:resolution:300;
%alpha=0:resolution:360;
beta=0:resolution:90;
%beta=0:resolution:180;
allPossibleFields.alphasing=alpha(Icolmax);
allPossibleFields.betasing=beta(Irowmax);
allPossibleFields.Icolmax=Icolmax;
allPossibleFields.Irowmax=Irowmax;
end


%% find max singularities in alpha beta real rotation map
if fieldType==1
mat180tmp=mat180;
% mat180tmp=matOSI;
for i=1:2
    if i==1
        [~,I]=max(mat180tmp(:));
        [Irow,Icol]=ind2sub(size(mat180tmp),I);
    else   
       try mat180tmp(Irow-(40/resolution):Irow+(40/resolution),Icol-(40/resolution):Icol+(40/resolution))=0; end 
       try mat180tmp(Irow-(30/resolution):Irow+(30/resolution),Icol-(30/resolution):Icol+(30/resolution))=0; end 
       try mat180tmp(Irow-(20/resolution):Irow+(20/resolution),Icol-(20/resolution):Icol+(20/resolution))=0; end 
       try mat180tmp(Irow-(10/resolution):Irow+(10/resolution),Icol-(10/resolution):Icol+(10/resolution))=0; end 
           % set to zero all elements that are 40 degrees away from the maximum points in both rows and columns
       [~,I]=max(mat180tmp(:));
       [Irow,Icol]=ind2sub(size(mat180tmp),I);
       
%        Icol=Icolmax(i-1)+(90/resolution);         % forces the alpha
%         % seperation between hotspots to be 90 degrees
%         if Icol>size(mat180,2)
%             Icol=Icol-size(mat180,2);
%         end
%         [~,Irow]=max(mat180(:,Icol)) ;    
    end
Icolmax(i)=Icol;
Irowmax(i)=Irow;
end
alpha=120:resolution:300;
%alpha=0:resolution:360;
beta=0:resolution:90;
%beta=0:resolution:180;
allPossibleFields.alphasing=alpha(Icolmax);
allPossibleFields.betasing=beta(Irowmax);
allPossibleFields.Icolmax=Icolmax;
allPossibleFields.Irowmax=Irowmax;
end
%% save allPossibleFields file
filedir=cd;
ind=strfind(filedir,'\');
filename=filedir;
filename(ind(end):length(filename))=[];
filename(1:ind(end-1))=[];

currenttime=datestr(now,'mmm_dd_yyyy_HH_MM');
switch histType
    
case 5    
switch ALPHACORRstate
    case 1
                switch fieldType
                    case 1
                        save(['allPossibleFieldsAlphaCorr_',num2str(filename),'_real_rotation_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');
                    case 0
                        save(['allPossibleFieldsAlphaCorr_',num2str(filename),'_real_translation_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');
                    case 2
                        save(['allPossibleFieldsAlphaCorr_',num2str(filename),'_real_mixed orthogonal_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');

                end   
    case 0
                switch fieldType
                    case 1
                        save(['allPossibleFields_',num2str(filename),'_real_rotation_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');
                    case 0
                        save(['allPossibleFields_',num2str(filename),'_real_translation_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');
                    case 2
                        save(['allPossibleFields_',num2str(filename),'_real_mixed orthogonal_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');

                end
end


case 6    
switch ALPHACORRstate
    case 1
                switch fieldType
                    case 1
                        save(['allPossibleFieldsAlphaCorr_',num2str(filename),'_fake_combinedMap_rotation_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');
                    case 0
                        save(['allPossibleFieldsAlphaCorr_',num2str(filename),'_fake_combinedMap_translation_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');
                    case 2
                        save(['allPossibleFieldsAlphaCorr_',num2str(filename),'_fake_combinedMap_mixed orthogonal_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');
end   
    case 0
                switch fieldType
                    case 1
                        save(['allPossibleFields_',num2str(filename),'_fake_combinedMap_rotation_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');
                    case 0
                        save(['allPossibleFields_',num2str(filename),'_fake_combinedMap_translation_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');
                    case 2
                        save(['allPossibleFields_',num2str(filename),'_fake_combinedMap_mixed orthogonal_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');

                end
end


case 7    
switch ALPHACORRstate
    case 1
                switch fieldType
                    case 1
                        save(['allPossibleFieldsAlphaCorr_',num2str(filename),'_single channel probed with_rotation_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');
                    case 0
                        save(['allPossibleFieldsAlphaCorr_',num2str(filename),'_single channel probed with_translation_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');
                    case 2
                        save(['allPossibleFieldsAlphaCorr_',num2str(filename),'_single channel probed with_mixed orthogonal_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');

                end   
    case 0
                switch fieldType
                    case 1
                        save(['allPossibleFields_',num2str(filename),'_single channel probed with_rotation_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');
                    case 0
                        save(['allPossibleFields_',num2str(filename),'_single channel probed with_translation_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');
                    case 2
                        save(['allPossibleFields_',num2str(filename),'_single channel probed with_mixed orthogonal_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'allPossibleFields');

                end
end

end


%% plot contours with singularities

switch fieldType
case 1          % rotation
    
%% calculate and plot rotation (case 1) or translation (case 2) singularities

%normalASC=[0.596,0.766,0.241];          %ASC
normalASC=[0.646079040832055;0.758000000000000;-0.080292421793915];
%normalASC=[0.604,0.758,0.243]; % ASC symmetric
[alphaASC,betaASC,~]=ABGlobalVisual2RetinaRight(normalASC);
betaiASC=90-abs(90-betaASC);
if alphaASC<=180
    alphaiASC=180+alphaASC;
else alphaASC>180 
    alphaiASC=alphaASC-180;
end

%normalPSC=[0.471,-0.766,0.438];         %PSC
normalPSC=[0.618330720986842;-0.760000000000000;0.193486742398264];
%normalPSC=[0.447,-0.76,0.469];    % PSC symmetric
[alphaPSC,betaPSC,~]=ABGlobalVisual2RetinaRight(normalPSC);
betaiPSC=90-abs(90-betaPSC);
if alphaPSC<=180
    alphaiPSC=180+alphaPSC;
else alphaPSC>180 
    alphaiPSC=alphaPSC-180;
end

%normalLSC=[0.421,-0.060,-0.905];        %LSC
normalLSC=[-0.036837075032666;-0.085000000000000;-0.992592580016110];
%normalLSC=[0.449,-0.085,-0.886];   % LSC symmetric
[alphaLSC,betaLSC,~]=ABGlobalVisual2RetinaRight(normalLSC);
betaiLSC=90-abs(90-betaLSC);
if alphaLSC<=180
    alphaiLSC=180+alphaLSC;
else alphaLSC>180 
    alphaiLSC=alphaLSC-180;
end


%% plot alpha-beta or azimuth-elevation map

% convert alpha beta map to azimuth elevation map for global space and plot it
if GSstate==1       % global space
    
    if fieldType==2 % mixed orthogonal
        elements=floor(size(matMix,2)/2);
    else
        elements=floor(size(mat180,2)/2);
    end

% plot alphabeta map
if fieldType==2 % mixed orthogonal
    d=circshift(matMix,[0 elements]);   %shift azimuth 0 to middle of graph
else
    d=circshift(mat180,[0 elements]);   %shift azimuth 0 to middle of graph
end

d(:,elements+1)=[];   % remove the duplicate column at zero
d=[d(:,end),d];   % add the last column to the start
f2=figure;
hold on
alpha=-180:resolution:180;
beta=-90:resolution:90;
contourf(alpha,beta,d,10);  
xlim([-180 180]);ylim([-90 90]);

% plot meanOutput map
% dOutput=circshift(allPossibleFields.meanOutput,[0 elements]);   %shift azimuth 0 to middle of graph
% dOutput(:,elements+1)=[];   % remove the duplicate column at zero
% dOutput=[dOutput(:,end),dOutput];   % add the last column to the start
% f3=figure;
% alpha=-180:resolution:180;
% beta=-90:resolution:90;
% contourf(alpha,beta,dOutput);
% xlim([-180 180]);ylim([-90 90]);

else
f2=figure;
alpha=120:resolution:300;
%alpha=0:resolution:360;
beta=0:resolution:90;
%beta=0:resolution:180;
%subplot(1,2,1);

if fieldType==2 % mixed orthogonal
    contourf(alpha,beta,matMix,10);
else
%     contourf(alpha,beta,matOSI,10);
    contourf(alpha,beta,mat180,10);
end
hold on;
if ABstate==2
scatter(alphaASC,betaASC,60,'bx','LineWidth',2);
scatter(alphaiASC,betaiASC,50,'bo','LineWidth',2);
scatter(alphaPSC,betaPSC,60,'gx','LineWidth',2);
scatter(alphaiPSC,betaiPSC,50,'go','LineWidth',2);
scatter(alphaLSC,betaLSC,60,'rx','LineWidth',2);
scatter(alphaiLSC,betaiLSC,50,'ro','LineWidth',2);
elseif ABstate==1
scatter(allPossibleFields.alphasing,allPossibleFields.betasing,'kx','LineWidth',2);    
end

% hold on;
% set(gcf, 'color', [1 1 1]);
% plot([0 360],[90 90],'--m','LineWidth',3);
% 
% %grid on
% ax = gca;
% ax.XTick = 0:45:360;
% ax.YTick = 0:45:180;
% ax.XTickLabel = {'N','45','D','135','T','225','V','315','N'};
% ax.FontSize=16;
% xlabel('Direction (deg.)','FontSize',20);
% ax.YTickLabel = {'','','','',''};
% 
% ylabh=ylabel('Eccentricity (deg.)','FontSize',20);
% set(ylabh, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
% ax.YTickLabel = {'anti-optic disc','135', 'margin','45','optic disc'};
% 
% c=colorbar('eastoutside');
% c.Label.String='Goodness of fit (%)';
% c.FontSize=14;


% plot meanOutput map
% f3=figure;
% alpha=120:resolution:300;
% beta=0:resolution:90;
% try contourf(alpha,beta,allPossibleFields.meanOutput); end

end

    case 0          % translation
        
f2=figure;
alpha=120:resolution:300;
%alpha=0:resolution:360;
beta=0:resolution:90;
% beta=0:resolution:180;
%subplot(1,2,1);
contourf(alpha,beta,mat180,10);
% contourf(alpha,beta,mat90,10);
% contourf(alpha,beta,matOSI,10);
hold on;
scatter(allPossibleFields.alphasing,allPossibleFields.betasing,'kx','LineWidth',3);


% hold on;
% set(gcf, 'color', [1 1 1]);
% plot([0 360],[90 90],'--m','LineWidth',3);
% 
% %grid on
% ax = gca;
% ax.XTick = 0:45:360;
% ax.YTick = 0:45:180;
% ax.XTickLabel = {'N','45','D','135','T','225','V','315','N'};
% ax.FontSize=12;
% xlabel('Direction (deg.)','FontSize',16);
% ax.YTickLabel = {'','','','',''};
% 
% ylabh=ylabel('Eccentricity (deg.)','FontSize',16);
% set(ylabh, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
% ax.YTickLabel = {'anti-optic disc','45', 'margin','135','optic disc'};

% c=colorbar('northoutside');
% c.Label.String='Goodness of fit (%)';
% c.FontSize=14;

save('tmp');
      

  case 2          % mixed orthogonal
        
f2=figure;
alpha=120:resolution:300;
% alpha=0:resolution:360;
beta=0:resolution:90;
% beta=0:resolution:180;
%subplot(1,2,1);
contourf(alpha,beta,matMix,10);
hold on;
scatter(allPossibleFields.alphasing,allPossibleFields.betasing,'kx','LineWidth',3);


% hold on;
% set(gcf, 'color', [1 1 1]);
% plot([0 360],[90 90],'--m','LineWidth',3);
% 
% %grid on
% ax = gca;
% ax.XTick = 0:45:360;
% ax.YTick = 0:45:180;
% ax.XTickLabel = {'N','45','D','135','T','225','V','315','N'};
% ax.FontSize=12;
% xlabel('Direction (deg.)','FontSize',16);
% ax.YTickLabel = {'','','','',''};
% 
% ylabh=ylabel('Eccentricity (deg.)','FontSize',16);
% set(ylabh, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
% ax.YTickLabel = {'anti-optic disc','45', 'margin','135','optic disc'};

% c=colorbar('northoutside');
% c.Label.String='Goodness of fit (%)';
% c.FontSize=14;

save('tmp');
end

%% save figure
currenttime=datestr(now,'mmm_dd_yyyy_HH_MM');
switch histType
    
case 5    
switch ALPHACORRstate
    case 1
                switch fieldType
                    case 1
                        savefig(f2,['allPossibleFieldsAlphaCorr_',num2str(filename),'_real_rotation_','_cutoff_',num2str(c),'_OS_created_',currenttime]);
                    case 0
                        saveas(f2,['allPossibleFieldsAlphaCorr_',num2str(filename),'_real_translation_','_cutoff_',num2str(c),'_OS_created_',currenttime],'fig');
                    case 2
                        saveas(f2,['allPossibleFieldsAlphaCorr_',num2str(filename),'_real_mixed orthogonal_','_cutoff_',num2str(c),'_OS_created_',currenttime],'fig');

                end   
    case 0
                switch fieldType
                    case 1
                        savefig(f2,['allPossibleFields_',num2str(filename),'_real_rotation_','_cutoff_',num2str(c),'_OS_created_',currenttime]);
                    case 0
                        savefig(f2,['allPossibleFields_',num2str(filename),'_real_translation_','_cutoff_',num2str(c),'_OS_created_',currenttime]);
                    case 2
                        savefig(f2,['allPossibleFields_',num2str(filename),'_real_mixed orthogonal_','_cutoff_',num2str(c),'_OS_created_',currenttime]);

                end
end


case 6    
switch ALPHACORRstate
    case 1
                switch fieldType
                    case 1
                        savefig(f2,['allPossibleFieldsAlphaCorr_',num2str(filename),'_fake_combinedMap_rotation_','_cutoff_',num2str(c),'_OS_created_',currenttime]);
                    case 0
                        savefig(f2,['allPossibleFieldsAlphaCorr_',num2str(filename),'_fake_combinedMap_translation_','_cutoff_',num2str(c),'_OS_created_',currenttime]);
                    case 2
                        savefig(f2,['allPossibleFieldsAlphaCorr_',num2str(filename),'_fake_combinedMap_mixed orthogonal_','_cutoff_',num2str(c),'_OS_created_',currenttime]);

                end   
    case 0
                switch fieldType
                    case 1
                        savefig(f2,['allPossibleFields_',num2str(filename),'_fake_combinedMap_rotation_','_cutoff_',num2str(c),'_OS_created_',currenttime]);
                    case 0
                        savefig(f2,['allPossibleFields_',num2str(filename),'_fake_combinedMap_translation_','_cutoff_',num2str(c),'_OS_created_',currenttime]);
                    case 2
                        savefig(f2,['allPossibleFields_',num2str(filename),'_fake_combinedMap_mixed orthogonal_','_cutoff_',num2str(c),'_OS_created_',currenttime]);

                end
end


case 7    
switch ALPHACORRstate
    case 1
                switch fieldType
                    case 1
                        savefig(f2,['allPossibleFieldsAlphaCorr_',num2str(filename),'_single channel probed with_rotation_','_cutoff_',num2str(c),'_OS_created_',currenttime]);
                        savefig(f3,['allPossibleFieldsAlphaCorr_',num2str(filename),'_single channel probed with_rotation_','mean output_','_OS_created_',currenttime]);
                    case 0
                        savefig(f2,['allPossibleFieldsAlphaCorr_',num2str(filename),'_single channel probed with_translation_','_cutoff_',num2str(c),'_OS_created_',currenttime]);
                    case 2
                        savefig(f2,['allPossibleFieldsAlphaCorr_',num2str(filename),'_single channel probed with_mixed orthogonal_','_cutoff_',num2str(c),'_OS_created_',currenttime]);

                end   
    case 0
                switch fieldType
                    case 1
                        savefig(f2,['allPossibleFields_',num2str(filename),'_single channel probed with_rotation_','_cutoff_',num2str(c),'_OS_created_',currenttime]);
                    case 0
                        savefig(f2,['allPossibleFields_',num2str(filename),'_single channel probed with_translation_','_cutoff_',num2str(c),'_OS_created_',currenttime]);
                    case 2
                        savefig(f2,['allPossibleFields_',num2str(filename),'_single channel probed with_mixed orthogonal_','_cutoff_',num2str(c),'_OS_created_',currenttime]);

                end
end



end

end

