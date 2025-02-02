%% import single-channel maps for 4 ONOFF translational fields 
[FileName,PathName] = uigetfile('*.mat','Select allPossibleFields files of the four single-chennel ONOFF maps','MultiSelect','on');
        for i=1:size(FileName,2)
            if strfind(FileName{i},'_alpha_0_beta_90')
                ONOFF0=load(FileName{i});
            elseif strfind(FileName{i},'_alpha_90_beta_90')
                ONOFF90=load(FileName{i});
            elseif strfind(FileName{i},'_alpha_180_beta_90')
                ONOFF180=load(FileName{i});
            elseif strfind(FileName{i},'_alpha_270_beta_90')
                ONOFF270=load(FileName{i});
            end
        end
        
[FileName,PathName] = uigetfile('*.mat','Select a pooledMapAlphaCorr file');
load(FileName);     % loads pooledmap file

%% 25% from each of the four ONOFF translation vector fields
f=randi(4,[size(ONOFF0.allPossibleFields.maxVec.X,1),1]);


MX=[ONOFF0.allPossibleFields.VecX,ONOFF90.allPossibleFields.VecX,ONOFF180.allPossibleFields.VecX,ONOFF270.allPossibleFields.VecX];
MY=[ONOFF0.allPossibleFields.VecY,ONOFF90.allPossibleFields.VecY,ONOFF180.allPossibleFields.VecY,ONOFF270.allPossibleFields.VecY];
MZ=[ONOFF0.allPossibleFields.VecZ,ONOFF90.allPossibleFields.VecZ,ONOFF180.allPossibleFields.VecZ,ONOFF270.allPossibleFields.VecZ];

for k=1:size(MX,1)
VecX(k,1)=MX(k,f(k));
VecY(k,1)=MY(k,f(k));
VecZ(k,1)=MZ(k,f(k));
end

figure;quiver3(pooledmap.XYZvecComp.X,pooledmap.XYZvecComp.Y,pooledmap.XYZvecComp.Z,VecX,VecY,VecZ)

maxVec.X=VecX;
maxVec.Y=VecY;
maxVec.Z=VecZ;


%% save single map
filedir=cd;
ind=strfind(filedir,'\');
filename=filedir;
filename(ind(end):length(filename))=[];
filename(1:ind(end-1))=[];
currenttime=datestr(now,'mmm_dd_yyyy_HH_MM');

discPoints=40;

save(['maxVec_',num2str(filename),'_translation_25percentAll4ONOFF_',num2str(discPoints),'_DS_created_',currenttime,'.mat'],'maxVec');
   
       

