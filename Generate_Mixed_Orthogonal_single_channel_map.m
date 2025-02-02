function []=Generate_Mixed_Orthogonal_single_channel_map()
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

VecXOtemp=(VecZ.*Y-VecY.*Z)./(VecY.*X-VecX.*Y);

VecYOtemp=(VecX.*Z-VecZ.*X)./(VecY.*X-VecX.*Y); 

VecX=VecXOtemp./(1+VecXOtemp.^2+VecYOtemp.^2).^(1/2);

VecY=VecYOtemp./(1+VecXOtemp.^2+VecYOtemp.^2).^(1/2);

VecZ=1./(1+VecXOtemp.^2+VecYOtemp.^2).^(1/2);

transField.allPossibleFields.VecX=VecX;
transField.allPossibleFields.VecY=VecY;
transField.allPossibleFields.VecZ=VecZ;


%% save allCells file

filedir=cd;
ind=strfind(filedir,'\');
filename=filedir;
filename(ind(end):length(filename))=[];
filename(1:ind(end-1))=[];
currenttime=datestr(now,'mmm_dd_yyyy_HH_MM');
allPossibleFields=transField.allPossibleFields;
save([num2str(FileName),'_Orthogonal','.mat'],'allPossibleFields');
end