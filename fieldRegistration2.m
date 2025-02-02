
function []= fieldRegistration2()

%% select files of the same fov to be registered
FolderName=uipickfiles;

%% sort file order by ascending ND and identify the middle file
for k=1:size(FolderName,2)
    names=getFileList(FolderName{k},'logFile',0,'anywhere');
    log{k}=load(num2str(names{1}));
    startTime(k)=datenum(log{k}.startTime);
end
clear log
[~,I]=sort(startTime);
FolderNameSorted=FolderName(I);
clear startTime
fixedID=find(I==ceil(median(I)));      % use the file that was acquired in in the middle

% load multi-tiff files and calculate projection
figure;
for j=1:size(FolderNameSorted,2)
    FileTif{j}=char(getFileList(FolderNameSorted{j},'.tif',0,'endonly'));
    
    warning('off','MATLAB:imagesci:tiffmexutils:libtiffErrorAsWarning');
    warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');
    
    InfoImage=imfinfo(FileTif{j});
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;
    NumberImages=length(InfoImage);
    stack=zeros(nImage,mImage,NumberImages,'uint16');
    
    TifLink = Tiff(FileTif{j}, 'r');
    for i=1:NumberImages
        TifLink.setDirectory(i);
        stack(:,:,i)=TifLink.read();
    end
    if sum(sum(stack(:,:,end)))==0
        stack(:,:,end)=[];
    end
    TifLink.close();
    Tstacks{j}=stack;
    
    proj{j}=uint16(mean(stack,3));
    subplot(4,3,j);
    imagesc(proj{1,j});
    axis equal;
    xlim([1 mImage]); ylim([1 nImage]);
end

% load projection images
f1=figure;
set(f1, 'color', [1 1 1]);
set(f1,'position',[150 150 1600 700]);
fixed=proj{fixedID};
for j=1:size(FolderNameSorted,2)
moving=proj{j};
subplot(1,3,1);
h=imshowpair(fixed,moving);

% calculate transformation for projection images
tformEstimate=imregcorr(moving,fixed);

% apply transformation on the moving projection image
Rfixed=imref2d(size(fixed));
movingReg=imwarp(moving,tformEstimate,'OutputView',Rfixed);
subplot(1,3,2);
imshowpair(fixed,movingReg,'falsecolor');

% apply transformation on the moving Tstack
stackReg=zeros(nImage,mImage,size(Tstacks{1,j},3));
Rfixed=imref2d(size(fixed));
for n=1:size(Tstacks{1,j},3)
    movingStack=Tstacks{1,j}(:,:,n);
    stackReg(:,:,n)=imwarp(movingStack,tformEstimate,'OutputView',Rfixed);  
end
TstackReg{j}=stackReg;

% calculate and plot the projection following registration
projReg{j}=uint16(mean(stackReg,3));
imagesc(projReg{1,j});
axis equal;
xlim([1 mImage]); ylim([1 nImage]);
subplot(1,3,3);
imshowpair(fixed, projReg{1,j},'Scaling','joint');


% save file
outputFileName=FileTif{j};
delete(outputFileName);     % delete the unregistered original tif file
img=uint16(TstackReg{j});
for K=1:length(img(1, 1, :))
   imwrite(img(:, :, K),outputFileName,'WriteMode','append','Compression','none');
end


ind=strfind(outputFileName,'\');
outputFileName1=([outputFileName(1:3),'\',outputFileName(4:27),'\',outputFileName(28:39),'\',outputFileName(40:47),'\',outputFileName(48:70),'\',outputFileName(71:end)]);
a=(['path=[',outputFileName,']']);
Miji(false);        % start Fiji
MIJ.run('Open...',{a});
MIJ.run('Green');
MIJ.run('Save');
MIJ.run('Close');

end

save('tmp.mat');

MIJ.exit;

%% register gr images of fields

Miji(false);        % start Fiji
% convert gr .oib files into .tif files 
filedir=cd;
grnames=getFileList(filedir,'.oib',0,'endonly');
if ~isempty(grnames)
    for g=1:length(grnames)
        thename=shortfile(grnames{g});
        a=(['path=[',filedir,'\',thename,']']);
        c=([a(1:end-4),'tif',']']);
        MIJ.run('Open...', {a});
        MIJ.run('Tiff...',{c});
        MIJ.run('Save', {c});
        MIJ.run('Close');
    end

% extract green and red channels into G and R
grnames=getFileList(filedir,'.tif',0,'endonly');
for g=1:length(grnames)    
    InfoImage=imfinfo(grnames{g});
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;
    NumberImages=length(InfoImage);
    grstack=zeros(nImage,mImage,NumberImages,'uint16');
    
    TifLink = Tiff(grnames{g}, 'r');
    for i=1:NumberImages
        TifLink.setDirectory(i);
        grstack(:,:,i)=TifLink.read();
    end
    TifLink.close();
    GRstack{g}=grstack;
    
    G{g}=GRstack{g}(:,:,1);
    R{g}=GRstack{g}(:,:,2);
    figure;
    subplot(1,2,1);
    imagesc(G{1,g});
    axis equal;
    xlim([1 mImage]); ylim([1 nImage]);
    subplot(1,2,2);
    imagesc(R{1,g});
    axis equal;
    xlim([1 mImage]); ylim([1 nImage]);
end

    
for g=1:length(grnames)   
moving=G{g};
% calculate transformation for green images
tformEstimate=imregcorr(moving,fixed);

% apply transformation on the moving (green) and red projection images
Rfixed=imref2d(size(fixed));
GmovingReg{g}=imwarp(moving,tformEstimate,'OutputView',Rfixed);
RmovingReg{g}=imwarp(R{g},tformEstimate,'OutputView',Rfixed);  

% save green and red files following registration
greenFileName=[grnames{g}(1:end-7),'_g','.tif'];
redFileName=[grnames{g}(1:end-7),'_r','.tif'];
imwrite(GmovingReg{g},greenFileName,'Compression','none');
imwrite(RmovingReg{g},redFileName,'Compression','none');
end

% apply green LUT on green image channel
% indg=strfind(greenFileName,'\');
% greenFileName=([greenFileName(1:indg(1)),'\',greenFileName(indg(1)+1:indg(2)),'\',greenFileName(indg(2)+1:indg(3)),'\',greenFileName(indg(3)+1:indg(4)),'\']);
a=(['path=[',greenFileName,']']);
Miji(false);        % start Fiji
MIJ.run('Open...',{a});
MIJ.run('Green');
MIJ.run('8-bit');
MIJ.run('Save');
MIJ.run('Close');

% apply red LUT on red image channel
% indr=strfind(redFileName,'\');
% redFileName=([redFileName(1:indr(1)),'\',redFileName(indr(1)+1:indr(2)),'\',redFileName(indr(2)+1:indr(3)),'\',redFileName(indr(3)+1:indr(4)),'\']);
a=(['path=[',redFileName,']']);
Miji(false);        % start Fiji
MIJ.run('Open...',{a});
MIJ.run('Red');
MIJ.run('8-bit');
MIJ.run('Save');
MIJ.run('Close');



end        
MIJ.exit;
end
