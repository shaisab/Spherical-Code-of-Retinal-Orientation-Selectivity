function calculateAllPossibleVecFieldsBoot(ALPHACORRstate,DATAMODELstate,CLUSTNstate,histType,FileName,resolution,fieldType,DISTstate,locDist,widthDist,ABstate,GSstate)

%% load data
load(FileName);     % loads pooledmap file OR combinedMap file OR single-channel map
data=pooledmap;

%%%%%% bootstrap %%%%%%

for nboot=1:1000
    nboot
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
for alpha=120:resolution:300
    for beta=0:resolution:90
        fieldID=['alpha_',num2str(alpha),'_','beta_',num2str(beta)];
        
        switch fieldType
            case 1
                allPossibleFields=VecComparisonPooled4Boot2(data,a,maxVec,deg2rad(alpha),deg2rad(beta),polarity,fieldID,vizAngle,allPossibleFields,k,h,interpType,histType,DISTstate,locDist,widthDist);
            case 0
                allPossibleFields=VecComparisonPooled4TransBoot2(data,a,maxVec,deg2rad(alpha),deg2rad(beta),polarity,fieldID,vizAngle,allPossibleFields,k,h,interpType,histType);
            case 2
                allPossibleFields=VecComparisonPooled4TransBoot2(data,a,maxVec,deg2rad(alpha),deg2rad(beta),polarity,fieldID,vizAngle,allPossibleFields,k,h,interpType,histType);

        end
        k=k+1;
    end
    k=1;
    h=h+1;
    %fieldID
end


%% generate matchind matrix
c=10;
cutoff=(['c',num2str(c)]);
k=1;
h=1;
for alpha=120:resolution:300
    for beta=0:resolution:90
        mat180(k,h)=allPossibleFields.matchind180(k,h).(cutoff);
        matMix(k,h)=allPossibleFields.matchindMix(k,h).(cutoff);
        k=k+1;
    end
    k=1;
    h=h+1;
end

%% find max singularities in alpha beta real translation map
if fieldType==0
    mat180tmp=mat180;
    
    % fix for individual retinas
    Icolmax=[1,2];
    Irowmax=[1,2];
    
    
    [Icolmax,J]=sort(Icolmax);      % sort Icolmax to assend along alpha
    Irowmax=Irowmax(J);             % sort Irowmax to assend along alpha
    
    alpha=120:resolution:300;
    beta=0:resolution:90;
    allPossibleFields.alphasing=alpha(Icolmax);
    allPossibleFields.betasing=beta(Irowmax);
    allPossibleFields.Icolmax=Icolmax;
    allPossibleFields.Irowmax=Irowmax;
end


%% find max singularities in alpha beta real mixed orthogonal map
if fieldType==2
    mat180tmp=matMix;
    
    % fix for individual retinas
    Icolmax=1;
    Irowmax=1;
    
    
    [Icolmax,J]=sort(Icolmax);      % sort Icolmax to assend along alpha
    Irowmax=Irowmax(J);             % sort Irowmax to assend along alpha
    
    alpha=120:resolution:300;
    beta=0:resolution:90;
    allPossibleFields.alphasing=alpha(Icolmax);
    allPossibleFields.betasing=beta(Irowmax);
    allPossibleFields.Icolmax=Icolmax;
    allPossibleFields.Irowmax=Irowmax;
end


%% find max singularities in alpha beta real rotation map
if fieldType==1
    mat180tmp=mat180;
    
    % fix for individual retinas
    Icolmax=[1,2];
    Irowmax=[1,2];
    
    [Icolmax,J]=sort(Icolmax);      % sort Icolmax to assend along alpha
    Irowmax=Irowmax(J); 
    
    alpha=120:resolution:300;
    beta=0:resolution:90;
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

switch fieldType
    case 1
        save(['allPossibleFieldsAlphaCorr_',num2str(filename),'_real_rotation_',num2str(discPoints),'_DS_created_',currenttime,'_nboot_',num2str(nboot),'.mat'],'allPossibleFields');
    case 0
        save(['allPossibleFieldsAlphaCorr_',num2str(filename),'_real_translation_',num2str(discPoints),'_DS_created_',currenttime,'_nboot_',num2str(nboot),'.mat'],'allPossibleFields');
    case 2
        save(['allPossibleFieldsAlphaCorr_',num2str(filename),'_real_mixed orthogonal_',num2str(discPoints),'_DS_created_',currenttime,'_nboot_',num2str(nboot),'.mat'],'allPossibleFields');
end

%% plot contours with singularities

% switch fieldType
%     case 1          % rotation
%         
%         f2=figure;
%         alpha=0:resolution:360;
%         beta=0:resolution:180;
%         contourf(alpha,beta,mat180,10);
%         hold on;
%         scatter(allPossibleFields.alphasing,allPossibleFields.betasing,'kx','LineWidth',2);
%         
%     case 0          % translation
%         
%         f2=figure;
%         alpha=0:resolution:360;
%         beta=0:resolution:180;
%         contourf(alpha,beta,mat180,10);
%         hold on;
%         scatter(allPossibleFields.alphasing,allPossibleFields.betasing,'kx','LineWidth',3);
%         hold off;
% end

end
end

