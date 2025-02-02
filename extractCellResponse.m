
% Should be called as: extractCellResponse(13) OR extractCellResponse([13,14])
% Should be in the Tstacks folder, where the mergeddata files are.

function []= extractCellResponse(cellnumber)

thisdir = dir;
str = {thisdir.name};
[s,v] = listdlg('PromptString','Select mergedmeta file :',...
    'SelectionMode','multiple',...
    'ListString',str);
if ~v disp('no files selected'); return; end

for i=1:size(s,2)
    [mergedmetaFileList{i}]=thisdir(s(i)).name;
end

for filen=1:size(mergedmetaFileList,2)
    load(mergedmetaFileList{filen});
    
    for ncell=1:size(cellnumber,2)
        
        [~,p]=max(mergedmeta.meanRepDG.ROIdata{1, cellnumber(ncell)}.summary.meanr);
        n=p+4;
        if n>8
            n=p-4;
        end
        
        pResponse(:,ncell)=mergedmeta.meanRepDB.ROIdata{1,cellnumber(ncell)}.summary.meannororddff(:,p);
        nResponse(:,ncell)=mergedmeta.meanRepDB.ROIdata{1,cellnumber(ncell)}.summary.meannororddff(:,n);
        Time(:,ncell)=mergedmeta.meanRepDB.ROIdata{1,cellnumber(ncell)}.summary.ordfluotime;
        pdir=mergedmeta.meanRepDG.ROIdata{1,cellnumber(ncell)}.summary.pdir;
        DSI=mergedmeta.meanRepDG.ROIdata{1,cellnumber(ncell)}.summary.DSI;
        
        

%%% calculate the onset and offset of specefic ROIs accounting for their location in the FOV    
    pR=strfind(mergedmeta.metadataDB{1,1}.imheader.Footnote.Resolution,'_');
    pixperum=str2num(mergedmeta.metadataDB{1,1}.imheader.Footnote.Resolution(pR(1)+1:pR(2)-1));
    widthum=mergedmeta.metadataDB{1,1}.imheader.General.SizeX/pixperum;
    heightum=mergedmeta.metadataDB{1,1}.imheader.General.SizeY/pixperum;  
    z=750;      % at non-right angles, if the bar start to move within the screen (rather than outside of it, as is the case for right angles).                     
    screenWidthum=3500;  % 3 mm at the plane of the retina
    pixelprojwidth=screenWidthum/600;   % 3 mm / 600 pixels
    radDendTreeum=300;
    corrSpeedFactor=1.2;
    delay=0;
    
    speedterm=pixelprojwidth*corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed;
    t=(screenWidthum-widthum)/2;
    c=radDendTreeum;
    e=widthum;
    f=heightum;
    h=screenWidthum;
    m=sqrt(h^2+h^2)/2;  % half diagonal of the screen 
    
    %%% extract the location of each ROI in the FOV          
    roiLocumX=mergedmeta.metadataDB{1,1}.ROIdata{1,cellnumber(ncell)}.centerX/pixperum;
    roiLocumY=mergedmeta.metadataDB{1,1}.ROIdata{1,cellnumber(ncell)}.centerY/pixperum;
    b=roiLocumX;
    k=roiLocumY;
    g=sqrt((f-k)^2+(b-e/2)^2);
    beta=acosd((b-e/2)/g);
    
    %%% calculate the onset and offset of specefic ROIs accounting for their location in the FOV   

    switch p*45-45
        case 0
            onset1=(delay+t+b-c)/speedterm;
            onset2=(delay+t+b)/speedterm;
        case 45
            x=g*sind(beta-45);
            n=m-x;
            onset1=(delay+n-z-c)/speedterm;
            onset2=(delay+n-z)/speedterm;
        case 180
            onset1=(delay+t+(e-b)-c)/speedterm;
            onset2=(delay+t+(e-b))/speedterm;
        case 225
            x=g*sind(beta-45);
            n=m+x;
            onset1=(delay+n-z-c)/speedterm;
            onset2=(delay+n-z)/speedterm;
        case 90
            onset1=(delay+t+k-c)/speedterm;
            onset2=(delay+t+k)/speedterm;
        case 135
            x=g*cosd(beta-45);
            n=m-x;
            onset1=(delay+n-z-c)/speedterm;
            onset2=(delay+n-z)/speedterm;
        case 270
            onset1=(delay+t+(e-f)+(f-k)-c)/speedterm;
            onset2=(delay+t+(e-f)+(f-k))/speedterm;
        case 315
            x=g*cosd(beta-45);
            n=m+x;
            onset1=(delay+n-z-c)/speedterm;
            onset2=(delay+n-z)/speedterm;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
    
    
        figure
        plot(Time, pResponse,'b');
        hold on
        plot(Time, nResponse,'r');
        hold on
        
    end
    hold off
end

curdir=cd;
inds=find(curdir=='\');
curdir=curdir(inds(2)+1:inds(3)-1);

curfov=cell2mat(mergedmeta.metaFileList{1});
curfov=curfov(1:6);


save([cd,'\identifiedCell_',curdir,'_',curfov,'_cell#',num2str(cellnumber),'.mat'],'Time','pResponse','nResponse','pdir','DSI','onset1','onset2');

end