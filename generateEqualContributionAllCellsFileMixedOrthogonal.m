function []= generateEqualContributionAllCellsFileMixedOrthogonal()
%function []= generateEqualContributionPooledmap(channelID,resolution,mapType)

discPoints=40;

%% Load files
[FileName,PathName] = uigetfile('*.mat','Select allPossibleFields files of single-channel translation or rotation maps','MultiSelect','on');

% sort file order by ascending alpha
if iscell(FileName)
    for k=1:size(FileName,2)
        indUnder=find(FileName{k}=='_');
        order(k)=str2num(FileName{k}(indUnder(4)+1:indUnder(5)-1));
    end
    [~,I]=sort(order);
    FileNameSorted=FileName(I);
    for i=1:size(FileName,2)
        transField{i}=load(FileNameSorted{i});
    end
else
    transField=load(FileName);
end

%% extract maxVec for each channel
X=[];
Y=[];
Z=[];
VecX=[];
VecY=[];
VecZ=[];
if iscell(FileName)
    for i=1:size(transField,2)
        X=[X;transField{1,i}.allPossibleFields.X];
        Y=[Y;transField{1,i}.allPossibleFields.Y];
        Z=[Z;transField{1,i}.allPossibleFields.Z];
        VecX=[VecX;transField{1,i}.allPossibleFields.VecX];
        VecY=[VecY;transField{1,i}.allPossibleFields.VecY];
        VecZ=[VecZ;transField{1,i}.allPossibleFields.VecZ];
    end
else
    X=transField.allPossibleFields.X;
    Y=transField.allPossibleFields.Y;
    Z=transField.allPossibleFields.Z;
    VecX=transField.allPossibleFields.VecX;
    VecY=transField.allPossibleFields.VecY;
    VecZ=transField.allPossibleFields.VecZ;
end

% % orthogonal vectors
% VecXOtemp=(VecZ.*Y-VecY.*Z)./(VecY.*X-VecX.*Y);
% 
% VecYOtemp=(VecX.*Z-VecZ.*X)./(VecY.*X-VecX.*Y);
% 
% VecXO=VecXOtemp./(1+VecXOtemp.^2+VecYOtemp.^2).^(1/2);
% 
% % VecXO=VecX;
% 
% VecYO=VecYOtemp./(1+VecXOtemp.^2+VecYOtemp.^2).^(1/2);
% 
% % VecYO=VecY;
% 
% VecZO=1./(1+VecXOtemp.^2+VecYOtemp.^2).^(1/2);

% VecZO=VecZ;

% denom=((Y.*VecZ-Z.*VecY).^2+(X.*VecZ-Z.*VecX).^2+(X.*VecY-Y.*VecX).^2).^(1/2);
% 
% VecXO=(Y.*VecZ-Z.*VecY)./denom;
% 
% VecYO=-1.*(X.*VecX-Z.*VecX)./denom;
% 
% VecZO=(X.*VecY-Y.*VecX)./denom;

phys=dlmread('physAlphaCorr');
R=phys(1);

VecXO=(VecY.*(Z-R)-VecZ.*Y)./R;

VecYO=(VecZ.*X-VecX.*(Z-R))./R;

VecZO=(VecX.*Y-VecY.*X)./R;

% figure
% quiver3(X,Y,Z,VecXO,VecYO,VecZO)
% hold on
% quiver3(X,Y,Z,VecX,VecY,VecZ)
% 

X=[X;X];
Y=[Y;Y];
Z=[Z;Z];
VecX=[VecX;VecXO];
VecY=[VecY;VecYO];
VecZ=[VecZ;VecZO];

% convert (x,y,z) coordinates to standard retina coordinates
[USt,VSt,VecUSt,VecVSt]=XYZtoStandard(X,Y,Z,VecX,VecY,VecZ,discPoints);

denom_2=((VecUSt).^2+(VecVSt).^2).^(1/2);

NormVecUSt=VecUSt./denom_2;
NormVecVSt=VecVSt./denom_2;

allCells.UVStvecComp.USt=[USt;USt];
allCells.UVStvecComp.VSt=[VSt;VSt];
allCells.UVStvecComp.VecUSt=[NormVecUSt';-1.*NormVecUSt'];
allCells.UVStvecComp.VecVSt=[NormVecVSt';-1.*NormVecVSt'];


Opposite_VecX=-1.*VecX;
Opposite_VecY=-1.*VecY;
Opposite_VecZ=-1.*VecZ;

Opposite_VecXO=-1.*VecXO;
Opposite_VecYO=-1.*VecYO;
Opposite_VecZO=-1.*VecZO;

Opposite_VecX=[Opposite_VecX;Opposite_VecXO];
Opposite_VecY=[Opposite_VecY;Opposite_VecYO];
Opposite_VecZ=[Opposite_VecZ;Opposite_VecZO];

allCells.XYZvecComp.X=[X;X];
allCells.XYZvecComp.Y=[Y;Y];
allCells.XYZvecComp.Z=[Z;Z];
allCells.XYZvecComp.VecX=[VecX;Opposite_VecX];
allCells.XYZvecComp.VecY=[VecY;Opposite_VecY];
allCells.XYZvecComp.VecZ=[VecZ;Opposite_VecZ];

allCells.cellID=num2cell(ones(2*size(X,1),1));         % for comapatibility with AllCellsAnalysis
allCells.injectionSite=num2cell(ones(2*size(X,1),1));  % for comapatibility with AllCellsAnalysis
allCells.pResponse=num2cell(ones(2*size(X,1),1));      % for comapatibility with AllCellsAnalysis
allCells.OSI=ones(2*size(X,1),1);                      % for comapatibility with AllCellsAnalysis
allCells.postmax=ones(2*size(X,1),1);                  % for comapatibility with AllCellsAnalysis
allCells.ONorONOFForOFF=ones(2*size(X,1),1);           % for comapatibility with AllCellsAnalysis
allCells.isCART=ones(2*size(X,1),1);                   % for comapatibility with AllCellsAnalysis
allCells.isRBPMS=ones(2*size(X,1),1);                  % for comapatibility with AllCellsAnalysis
allCells.isRetro=ones(2*size(X,1),1);                  % for comapatibility with AllCellsAnalysis
allCells.onset1=ones(2*size(X,1),1);                   % for comapatibility with AllCellsAnalysis
allCells.symmetry_ratio=ones(2*size(X,1),1);           % for comapatibility with AllCellsAnalysis
allCells.sumCorrRep2=ones(2*size(X,1),1);              % for comapatibility with AllCellsAnalysis
allCells.poriSz=ones(2*size(X,1),1);                   % for comapatibility with AllCellsAnalysis
allCells.oriR2=ones(2*size(X,1),1);                    % for comapatibility with AllCellsAnalysis
allCells.original_rmse=ones(2*size(X,1),1);            % for comapatibility with AllCellsAnalysis
allCells.shuffle_rmse=ones(2*size(X,1),1);             % for comapatibility with AllCellsAnalysis
allCells.theta=ones(2*size(X,1),9);                    % for comapatibility with AllCellsAnalysis 
allCells.meanr=ones(2*size(X,1),8);                    % for comapatibility with AllCellsAnalysis
allCells.meanrs=ones(2*size(X,1),8);                   % for comapatibility with AllCellsAnalysis 

allCells.XYZvecComp.VecXcompASC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.XYZvecComp.VecYcompASC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.XYZvecComp.VecZcompASC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.XYZvecComp.angleXYZASC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.XYZvecComp.VecXcompPSC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.XYZvecComp.VecYcompPSC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.XYZvecComp.VecZcompPSC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.XYZvecComp.angleXYZPSC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.XYZvecComp.VecXcompLSC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.XYZvecComp.VecYcompLSC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.XYZvecComp.VecZcompLSC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.XYZvecComp.angleXYZLSC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis

allCells.UVvecComp.U=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.UVvecComp.V=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.UVvecComp.VecU=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.UVvecComp.VecV=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.UVvecComp.VecUcompASC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.UVvecComp.VecVcompASC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.UVvecComp.angleUVASC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.UVvecComp.VecUcompPSC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.UVvecComp.VecVcompPSC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.UVvecComp.angleUVPSC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.UVvecComp.VecUcompLSC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.UVvecComp.VecVcompLSC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.UVvecComp.angleUVLSC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis

allCells.UVStvecComp.VecUStcompASC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.UVStvecComp.VecVStcompASC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.UVStvecComp.angleUVStASC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.UVStvecComp.VecUStcompPSC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.UVStvecComp.VecVStcompPSC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.UVStvecComp.angleUVStPSC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.UVStvecComp.VecUStcompLSC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.UVStvecComp.VecVStcompLSC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis
allCells.UVStvecComp.angleUVStLSC=ones(2*size(X,1),1);  % for comapatibility with AllCellsAnalysis

%% save allCells file

discPoints=40;

filedir=cd;
ind=strfind(filedir,'\');
filename=filedir;
filename(ind(end):length(filename))=[];
filename(1:ind(end-1))=[];
currenttime=datestr(now,'mmm_dd_yyyy_HH_MM');
save(['allCellsAlphaCorr_',num2str(filename),'_',num2str(discPoints),'_equal channel contribution and uniform distribution','_created_',currenttime,'.mat'],'allCells');

end 