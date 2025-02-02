
function []= registerTandZstacks()

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
    TifLink.close();
    Tstacks{j}=stack;
    
    proj{j}=uint16(mean(stack,3));
    subplot(3,3,j);
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

% calculate transformation projection images
tformEstimate=imregcorr(moving,fixed);

% apply transformation on the moving projection image
Rfixed=imref2d(size(fixed));
movingReg=imwarp(moving,tformEstimate,'OutputView',Rfixed);
subplot(1,3,2);
imshowpair(fixed,movingReg,'falsecolor');

% apply transformation on the moving Tstack
stackReg=zeros(nImage,mImage,NumberImages);
Rfixed=imref2d(size(fixed));
for n=1:NumberImages
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
outputFileName=([outputFileName(1:3),'\',outputFileName(4:27),'\',outputFileName(28:39),'\',outputFileName(40:47),'\',outputFileName(48:70),'\',outputFileName(71:end)]);
a=(['path=[',outputFileName,']']);
Miji(false);        % start Fiji
MIJ.run('Open...',{a});
MIJ.run('Green');
MIJ.run('Save');
MIJ.run('Close');

end

MIJ.exit;
end
