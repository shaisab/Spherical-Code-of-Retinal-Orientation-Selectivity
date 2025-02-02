function combineSingleMaps(polarity2,mapType,DISTstate,locDist,widthDist,ABstate)

%% import pooledmap file
[FileName,PathName] = uigetfile('*.mat','Select a pooledMapAlphaCorr file');
load(FileName);     % loads pooledmap file


%% import nessecary single channel maps
switch mapType  
    case 1
        switch ABstate
            case 2  
                switch polarity2
                    case 2      %import single-channel maps for 6 rotation fields       
                        [FileName,PathName] = uigetfile('*.mat','Select allPossibleFields files of the six single-chennel rotation maps','MultiSelect','on');
                        for i=1:size(FileName,2)
                            if strfind(FileName{i},'_ASC')
                                ASC=load(FileName{i});
                            elseif strfind(FileName{i},'_PSC')
                                PSC=load(FileName{i});
                            elseif strfind(FileName{i},'_LSC')
                                LSC=load(FileName{i});
                            elseif strfind(FileName{i},'_iASC')
                                iASC=load(FileName{i});
                            elseif strfind(FileName{i},'_iPSC')
                                iPSC=load(FileName{i});
                            elseif strfind(FileName{i},'_iLSC')
                                iLSC=load(FileName{i});
                            end
                        end
                        
                    case 3      %import single-channel maps for 4 ASC and LSC rotation fields
                        
                        [FileName,PathName] = uigetfile('*.mat','Select allPossibleFields files of the four single-chennel ASC and LSC rotation maps','MultiSelect','on');
                        for i=1:size(FileName,2)
                            if strfind(FileName{i},'_ASC')
                                ASC=load(FileName{i});
                            elseif strfind(FileName{i},'_LSC')
                                LSC=load(FileName{i});
                            elseif strfind(FileName{i},'_iASC')
                                iASC=load(FileName{i});
                            elseif strfind(FileName{i},'_iLSC')
                                iLSC=load(FileName{i});
                            end
                        end
                end
                     
            case 1
                [FileName,PathName] = uigetfile('*.mat','Select allPossibleFields files of the four single-channel rotation maps based on real singularities','MultiSelect','on');
                for i=1:size(FileName,2)
                    transField{i}=load(FileName{i});
                end
        end
        
    case 0      %import single-channel maps for 4 ONOFF translational fields
        
        [FileName,PathName] = uigetfile('*.mat','Select allPossibleFields files of the four single-chennel ONOFF maps','MultiSelect','on');
        for i=1:size(FileName,2)
            transField{i}=load(FileName{i});
        end
end

%% generate combined map


switch mapType
    
    case 0       % 25% from each of the four translation vector fields
        
        %f=randi(4,[size(transField{1}.allPossibleFields.VecX,1),1]);
        f=randi(2,[size(transField{1}.allPossibleFields.VecX,1),1]);
        
%         MX=[transField{1}.allPossibleFields.VecX,transField{2}.allPossibleFields.VecX,transField{3}.allPossibleFields.VecX,transField{4}.allPossibleFields.VecX];
%         MY=[transField{1}.allPossibleFields.VecY,transField{2}.allPossibleFields.VecY,transField{3}.allPossibleFields.VecY,transField{4}.allPossibleFields.VecY];
%         MZ=[transField{1}.allPossibleFields.VecZ,transField{2}.allPossibleFields.VecZ,transField{3}.allPossibleFields.VecZ,transField{4}.allPossibleFields.VecZ];
        
        MX=[transField{1}.allPossibleFields.VecX,transField{2}.allPossibleFields.VecX];
        MY=[transField{1}.allPossibleFields.VecY,transField{2}.allPossibleFields.VecY];
        MZ=[transField{1}.allPossibleFields.VecZ,transField{2}.allPossibleFields.VecZ];
        
        for k=1:size(MX,1)
            VecX(k,1)=MX(k,f(k));
            VecY(k,1)=MY(k,f(k));
            VecZ(k,1)=MZ(k,f(k));
        end
        
        figure;quiver3(pooledmap.XYZvecComp.X,pooledmap.XYZvecComp.Y,pooledmap.XYZvecComp.Z,VecX,VecY,VecZ)
        
        combinedMap.X=pooledmap.XYZvecComp.X;
        combinedMap.Y=pooledmap.XYZvecComp.Y;
        combinedMap.Z=pooledmap.XYZvecComp.Z;
        combinedMap.VecX=VecX;
        combinedMap.VecY=VecY;
        combinedMap.VecZ=VecZ;
        
    case 1      % rotation
        
        if DISTstate==0         % generate combined rotation map WITHOUT distance function
            
            switch polarity2
                
                case 2       % 16.7% from each of the six rotation vector fields
                    
                    f=randi(6,[size(ASC.allPossibleFields.VecX,1),1]);
                    
                    MX=[ASC.allPossibleFields.VecX,PSC.allPossibleFields.VecX,LSC.allPossibleFields.VecX,iASC.allPossibleFields.VecX,iPSC.allPossibleFields.VecX,iLSC.allPossibleFields.VecX];
                    MY=[ASC.allPossibleFields.VecY,PSC.allPossibleFields.VecY,LSC.allPossibleFields.VecY,iASC.allPossibleFields.VecY,iPSC.allPossibleFields.VecY,iLSC.allPossibleFields.VecY];
                    MZ=[ASC.allPossibleFields.VecZ,PSC.allPossibleFields.VecZ,LSC.allPossibleFields.VecZ,iASC.allPossibleFields.VecZ,iPSC.allPossibleFields.VecZ,iLSC.allPossibleFields.VecZ];
                    
                    for k=1:size(MX,1)
                        VecX(k,1)=MX(k,f(k));
                        VecY(k,1)=MY(k,f(k));
                        VecZ(k,1)=MZ(k,f(k));
                    end
                    
                    figure;quiver3(pooledmap.XYZvecComp.X,pooledmap.XYZvecComp.Y,pooledmap.XYZvecComp.Z,VecX,VecY,VecZ)
                    
                    combinedMap.X=pooledmap.XYZvecComp.X;
                    combinedMap.Y=pooledmap.XYZvecComp.Y;
                    combinedMap.Z=pooledmap.XYZvecComp.Z;
                    combinedMap.VecX=VecX;
                    combinedMap.VecY=VecY;
                    combinedMap.VecZ=VecZ;
                    
                    
                case 3      % 25% from each of the four ASC and LSC rotation vector fields
                    
                    f=randi(4,[size(ASC.allPossibleFields.VecX,1),1]);
                    
                    MX=[ASC.allPossibleFields.VecX,LSC.allPossibleFields.VecX,iASC.allPossibleFields.VecX,iLSC.allPossibleFields.VecX];
                    MY=[ASC.allPossibleFields.VecY,LSC.allPossibleFields.VecY,iASC.allPossibleFields.VecY,iLSC.allPossibleFields.VecY];
                    MZ=[ASC.allPossibleFields.VecZ,LSC.allPossibleFields.VecZ,iASC.allPossibleFields.VecZ,iLSC.allPossibleFields.VecZ];
                    
                    for k=1:size(MX,1)
                        VecX(k,1)=MX(k,f(k));
                        VecY(k,1)=MY(k,f(k));
                        VecZ(k,1)=MZ(k,f(k));
                    end
                    
                    figure;quiver3(pooledmap.XYZvecComp.X,pooledmap.XYZvecComp.Y,pooledmap.XYZvecComp.Z,VecX,VecY,VecZ)
                    
                    combinedMap.X=pooledmap.XYZvecComp.X;
                    combinedMap.Y=pooledmap.XYZvecComp.Y;
                    combinedMap.Z=pooledmap.XYZvecComp.Z;
                    combinedMap.VecX=VecX;
                    combinedMap.VecY=VecY;
                    combinedMap.VecZ=VecZ;
            end
            
            
            
            
        elseif DISTstate==1         % generate combined rotation map WITH distance function
            
            
            switch polarity2
                
                case 2       % 16.7% from each of the six rotation vector fields
                    
                    %%% APLLY DISTANCE FUNCTION ---  NEED TO REPEAT THIS
                    %%% FOR EACH CANAL
                    phys=dlmread('phys');
                    R=phys(1);
                    
                    alphaCanal=deg2rad(alphaCanal);
                    betaCanal=deg2rad(betaCanal);
                    scoreDist=distanceFunction(alphaCanal,betaCanal,pooledmap.XYZvecComp.X,pooledmap.XYZvecComp.Y,pooledmap.XYZvecComp.Z,R,locDist,widthDist);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    
                    f=randi(6,[size(ASC.allPossibleFields.VecX,1),1]);
                    
                    MX=[ASC.allPossibleFields.VecX,PSC.allPossibleFields.VecX,LSC.allPossibleFields.VecX,iASC.allPossibleFields.VecX,iPSC.allPossibleFields.VecX,iLSC.allPossibleFields.VecX];
                    MY=[ASC.allPossibleFields.VecY,PSC.allPossibleFields.VecY,LSC.allPossibleFields.VecY,iASC.allPossibleFields.VecY,iPSC.allPossibleFields.VecY,iLSC.allPossibleFields.VecY];
                    MZ=[ASC.allPossibleFields.VecZ,PSC.allPossibleFields.VecZ,LSC.allPossibleFields.VecZ,iASC.allPossibleFields.VecZ,iPSC.allPossibleFields.VecZ,iLSC.allPossibleFields.VecZ];
                    
                    for k=1:size(MX,1)
                        VecX(k,1)=MX(k,f(k));
                        VecY(k,1)=MY(k,f(k));
                        VecZ(k,1)=MZ(k,f(k));
                    end
                    
                    figure;quiver3(pooledmap.XYZvecComp.X,pooledmap.XYZvecComp.Y,pooledmap.XYZvecComp.Z,VecX,VecY,VecZ)
                    
                    combinedMap.X=pooledmap.XYZvecComp.X;
                    combinedMap.Y=pooledmap.XYZvecComp.Y;
                    combinedMap.Z=pooledmap.XYZvecComp.Z;
                    combinedMap.VecX=VecX;
                    combinedMap.VecY=VecY;
                    combinedMap.VecZ=VecZ;
                    
                    
                case 3      % 25% from each of the four ASC and LSC rotation vector fields
                    
                    f=randi(4,[size(ASC.allPossibleFields.VecX,1),1]);
                    
                    MX=[ASC.allPossibleFields.VecX,LSC.allPossibleFields.VecX,iASC.allPossibleFields.VecX,iLSC.allPossibleFields.VecX];
                    MY=[ASC.allPossibleFields.VecY,LSC.allPossibleFields.VecY,iASC.allPossibleFields.VecY,iLSC.allPossibleFields.VecY];
                    MZ=[ASC.allPossibleFields.VecZ,LSC.allPossibleFields.VecZ,iASC.allPossibleFields.VecZ,iLSC.allPossibleFields.VecZ];
                    
                    for k=1:size(MX,1)
                        VecX(k,1)=MX(k,f(k));
                        VecY(k,1)=MY(k,f(k));
                        VecZ(k,1)=MZ(k,f(k));
                    end
                    
                    figure;quiver3(pooledmap.XYZvecComp.X,pooledmap.XYZvecComp.Y,pooledmap.XYZvecComp.Z,VecX,VecY,VecZ)
                    
                    combinedMap.X=pooledmap.XYZvecComp.X;
                    combinedMap.Y=pooledmap.XYZvecComp.Y;
                    combinedMap.Z=pooledmap.XYZvecComp.Z;
                    combinedMap.VecX=VecX;
                    combinedMap.VecY=VecY;
                    combinedMap.VecZ=VecZ;
            end
        end
        
end


%% save single map
filedir=cd;
ind=strfind(filedir,'\');
filename=filedir;
filename(ind(end):length(filename))=[];
filename(1:ind(end-1))=[];
currenttime=datestr(now,'mmm_dd_yyyy_HH_MM');

discPoints=40;

switch mapType
    case 0
        save(['combinedMap_',num2str(filename),'_translation_4channels_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'combinedMap');
    case 1
        switch polarity2
            case 2      %save combined map for 6 rotation fields
                save(['combinedMap_',num2str(filename),'_rotation_6channels_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'combinedMap');
            case 3      %save combined map for ASC and LSC rotation fields
                save(['combinedMap_',num2str(filename),'_rotation_ASCLSC_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'combinedMap');
        end
end

end
