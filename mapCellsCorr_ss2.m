
%need to be in the parent of the retina folder
%'type' should be 'DS' for DS cells or 'OS' for OS cells
%the function should be called as:  mapCells_ss('DS') OR  mapCells_ss('OS')


function []= mapCellsCorr_ss2()
%function []= mapCellsCorr_ss2(celltype)

Miji(false);

[filedir,go]=bbdirselector('select data folder',cd);
if ~go disp('no folder selected'); return; end

names=getFileList(filedir,'meridian_radius',0,'anywhere');
if ~isempty(names)
    MR=load(num2str(names{1}));
end

try
    names=getFileList(filedir,'injection_site',0,'anywhere');
    if ~isempty(names)
        [inject1,inject2,inject3,inject4]=textread(num2str(names{1}),'%s %s %s %s',1);
    end
    inject=[inject1,inject2,inject3,inject4];
end


try
    namesMap=getFileList(filedir,'map_',0,'anywhere');
    for i=1:size(namesMap,2)
        if strfind(namesMap{i},'.mat')
            prevMap=load(num2str(namesMap{i}));
        end
    end
end

%thedir=cd;
names=getFileList(filedir,'.oib',0,'endonly');
if ~isempty(names)
    %cd(thedir);
    for i = 1:length(names)
        thename=shortfile(names{i});
        a=(['path=[',thename,']']);
        b=([a(1:end-4),'txt',']']);
        c=([a(1:end-4),'tif',']']);
        MIJ.run('Open...', {a});
        MIJ.run('Tiff...',{c});
        MIJ.run('Show Info...');
        pause(1);
        MIJ.run('Text...',{b});
        MIJ.run('Close All');
    end
    MIJ.exit;
end

d=b(7:end-1);
% THIS A TMP FIX FOR JULY 12. IT SHOULD BE COMMENTED.
if strfind(filedir,'July_12_2014')
    header=parseHeaderRegTif('retinaTmp.txt');
else
    header=parseHeaderRegTif(d);
end

e=c(7:end-1);

% sizeX=header.General.SizeX;
% sizeY=header.General.SizeY;
pR=strfind(header.Footnote.Resolution,'_');
pixelsize=1/str2num(header.Footnote.Resolution(pR(1)+1:pR(2)-1));     % this is pixel size in um in the stitched retina image

f1=figure;
retina=imread(e);
imshow(retina,'DisplayRange',[0 4095]);

if exist('prevMap','var')
    
    A=prevMap.map.angCorr;
else
    
    hold on
    title('Select the left cut ','FontWeight','bold')
    [xL,yL]=ginput(1);
    plot(xL(1),yL(1),'Marker','o','MarkerSize',20,'MarkerEdgeColor','m');
    title('');
    
    hold on
    title('Select the right cut ','FontWeight','bold')
    [xR,yR]=ginput(1);    %first point is the center
    plot(xR(1),yR(1),'Marker','o','MarkerSize',20,'MarkerEdgeColor','m');
    title('')
    
    L = createLine([xR  yR], [xL yL]);
    A=180-rad2deg(lineAngle(L));
end

hold on
retinaRot=imrotate(retina,-A,'nearest','crop');
imshow(retinaRot,'DisplayRange',[0 4095]);
%scatter(xL,yL,100,'m');
%scatter(xR,yR,100,'m');
%drawLine(L, 'color', 'm', 'linewidth', 2);

hold on
title('Select the center point ','FontWeight','bold')
[xc,yc]=ginput(1);    %first point is the center
plot(xc(1),yc(1),'Marker','o','MarkerSize',20,'MarkerEdgeColor','y');
title('')
% xotrans=0+xc;
% yotrans=0+yc;
scatter(xc,yc);
Lh = createLine([xc, yc],[1000, yc]);
drawLine(Lh, 'color', 'c', 'linewidth', 2);

hold on
title('Select cut #1 ','FontWeight','bold')
[x1,y1]=ginput(1);
plot(x1(1),y1(1),'Marker','o','MarkerSize',10,'MarkerEdgeColor','m');
title('');

hold on
title('Select cut #2 ','FontWeight','bold')
[x2,y2]=ginput(1);
plot(x2(1),y2(1),'Marker','o','MarkerSize',10,'MarkerEdgeColor','m');
title('')

hold on
title('Select cut #3 ','FontWeight','bold')
[x3,y3]=ginput(1);
plot(x3(1),y3(1),'Marker','o','MarkerSize',10,'MarkerEdgeColor','m');
title('');

hold on
title('Select cut #4 ','FontWeight','bold')
[x4,y4]=ginput(1);
plot(x4(1),y4(1),'Marker','o','MarkerSize',10,'MarkerEdgeColor','m');
title('');

R1 = createRay([xc yc],[x1  y1]);
A1=360-rad2deg(lineAngle(R1));
if A1<0 A1=A1+360; end

R2 = createRay([xc yc],[x2  y2]);
A2=360-rad2deg(lineAngle(R2));
if A2<0 A2=A2+360; end

R3 = createRay([xc yc],[x3  y3]);
A3=360-rad2deg(lineAngle(R3));
if A3<0 A3=A3+360; end

R4 = createRay([xc yc],[x4  y4]);
A4=360-rad2deg(lineAngle(R4));
if A4<0 A4=A4+360; end

if A1>A2 A1=A1-360; end

D1=pixelsize*distancePoints([x1 y1],[xc yc]);
D2=pixelsize*distancePoints([x2 y2],[xc yc]);
D3=pixelsize*distancePoints([x3 y3],[xc yc]);
D4=pixelsize*distancePoints([x4 y4],[xc yc]);

drawRay(R1, 'color', 'm', 'linewidth', 2);
drawRay(R2, 'color', 'm', 'linewidth', 2);
drawRay(R3, 'color', 'm', 'linewidth', 2);
drawRay(R4, 'color', 'm', 'linewidth', 2);

%drawCircle(xc,yc, D);
map.angCorr=A;
map.a1=deg2rad(A1);
map.a2=deg2rad(A2);
map.a3=deg2rad(A3);
map.a4=deg2rad(A4);
map.D1=D1;
map.D2=D2;
map.D3=D3;
map.D4=D4;
map.xc=xc;
map.yc=yc;
map.x1=x1;
map.y1=y1;
map.x2=x2;
map.y2=y2;
map.x3=x3;
map.y3=y3;
map.x4=x4;
map.y4=y4;

[metaFileList]=getFileList(filedir,'metaRep',0,'anywhere');
[mergedmetaFileList]=getFileList(filedir,'mergedmeta',0,'anywhere');

if size(metaFileList,2)>0 && size(mergedmetaFileList,2)==0
    j=1;
    for fov=1:size(metaFileList,2)
        load(metaFileList{fov});
        
        %  map(1,:)={'filename', 'umroiX', 'umroiY',  'thetapixel', 'radiuspixel', 'pdir', 'DSI'};
        
        for ROIn=1:size(metaRep.metadata{1}.ROIdata,2)
            
            %             theta=deg2rad(metaRep.metadata{1,1}.ROIdata{1,ROIn}.theta);
            %             thetaSz=theta-deg2rad(A);  % to be use when combining data from several retinas.
            %             if thetaSz<0
            %                 thetaSz=2*pi+thetaSz;
            %             end
            %             radiuspixel=metaRep.metadata{1,1}.ROIdata{1,ROIn}.radius/pixelsize;
            %             [x, y]=pol2cart(thetaSz, radiuspixel);
            %             xctrans(ROIn)=xc+x;
            %             yctrans(ROIn)=yc-y;
            
            pdir=mergedmeta.meanRepDG.ROIdata{1,ROIn}.summary.pdir; % to be use when combining data from several retinas.
            pdirSz=pdir-A;
            if pdirSz<0
                pdirSz=360+pdirSz;
            end
            
            pori=mergedmeta.meanRepDG.ROIdata{1,ROIn}.summary.pori; % to be use when combining data from several retinas.
            poriSz=pori-A;
            if poriSz<0
                poriSz=360+poriSz;
            end
            if poriSz>180
                poriSz=poriSz-180;
            end
            
            map.cellID(j,:)=([metaRep.metadata{1,1}.filename(:,1:15),'_',mat2str(ROIn)]);
            %             map.theta(j)=rad2deg(theta);
            %             map.thetaSz(j)=rad2deg(thetaSz);
            %             map.radiuspixel(j)=radiuspixel;
            %             map.radiusum(j)=mergedmeta.metadataDG{1,1}.ROIdata{1,ROIn}.radius;
            map.umroiX(j)=metaRep.metadata{1,1}.ROIdata{1,ROIn}.umroiX;
            map.umroiY(j)=metaRep.metadata{1,1}.ROIdata{1,ROIn}.umroiY;
            map.pdir(j)=pdir;
            map.pdirSz(j)=pdirSz;
            map.pori(j)=pori;
            map.poriSz(j)=poriSz;
            map.DSI(j)=metaRep.meanRep.ROIdata{1,ROIn}.meanorddff.DSI;
            map.angCorr=A;
            map.umCutLength=D*pixelsize;
            
            %         map(ROIn,:)={'filename', 'umroiX', 'umroiY',  'thetapixel', 'radiuspixel', 'pdir', 'DSI'};
            %         map(ROIn+1,:)={filename, umroiX, umroiY,  thetapixel, radiuspixel, pdir, DSI};
            j=j+1;
        end
        
        plot(xctrans,yctrans,'o','MarkerSize',2,'MarkerEdgeColor','r');
        
        XYSz=rotateVector([map.umroiX map.umroiY], deg2rad(-A));
        map.umroiXSz=XYSz(:,1);
        map.umroiYSz=XYSz(:,2);
        
    end
        
elseif size(metaFileList,2)==0 && size(mergedmetaFileList,2)>0
    j=1;
    for fov=1:size(mergedmetaFileList,2)
        load(mergedmetaFileList{fov});
   
        %%% initialize parameters for the calculation of the onset of the stimulus
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
        
        for ROIn=1:size(mergedmeta.metadataDG{1}.ROIdata,2)
            
            pdir=mergedmeta.meanRepDG.ROIdata{1,ROIn}.summary.pdir; % to be use when combining data from several retinas.
            pdirSz=pdir+A;
            if pdirSz<0
                pdirSz=360+pdirSz;
            end
            
            pori=mergedmeta.meanRepDG.ROIdata{1,ROIn}.summary.pori; % to be use when combining data from several retinas.
            poriSz=pori+A;
            if poriSz<0
                poriSz=360+poriSz;
            end
            if poriSz>180
                poriSz=poriSz-180;
            end
            
            map.cellID{j,1}=([mergedmeta.metadataDG{1,1}.filename(:,1:15),'_',mat2str(ROIn)]);
            map.umroiX(j,1)=mergedmeta.metadataDG{1,1}.ROIdata{1,ROIn}.umroiX;
            map.umroiY(j,1)=mergedmeta.metadataDG{1,1}.ROIdata{1,ROIn}.umroiY;
            map.DSI(j,1)=mergedmeta.meanRepDG.ROIdata{1,ROIn}.summary.DSI;
            map.isDS(j,1)=mergedmeta.isDS(ROIn);
            map.isOS(j,1)=mergedmeta.isOS(ROIn);
            map.ONorONOFForOFF(j,1)=mergedmeta.ONorONOFForOFF(ROIn);
            try
                map.isCART(j,1)=mergedmeta.isCART(ROIn);
                map.isRBPMS(j,1)=mergedmeta.isRBPMS(ROIn);
                map.isRetro(j,1)=mergedmeta.isRetro(ROIn);
            end
            map.pdir(j,1)=pdir;
            map.pdirSz(j,1)=pdirSz;
            map.pori(j,1)=pori;
            map.poriSz(j,1)=poriSz;
            map.DSI(j,1)=mergedmeta.meanRepDG.ROIdata{1,ROIn}.summary.DSI;
            map.OSI(j,1)=mergedmeta.meanRepDG.ROIdata{1,ROIn}.summary.OSI;
            try
                map.globalOSI(j,1)=mergedmeta.meanRepDG.ROIdata{1,ROIn}.summary.globalOSI;
                map.oriFWHM(j,1)=mergedmeta.meanRepDG.ROIdata{1,ROIn}.summary.oriFWHM;
                map.oriR2(j,1)=mergedmeta.meanRepDG.ROIdata{1,ROIn}.summary.oriR2;
            end
            
            %%% extract the location of each ROI in the FOV
            map.roiLocumX(j,1)=mergedmeta.metadataDB{1,1}.ROIdata{1,ROIn}.centerX/pixperum;
            map.roiLocumY(j,1)=mergedmeta.metadataDB{1,1}.ROIdata{1,ROIn}.centerY/pixperum;
            b=map.roiLocumX(j,1);
            k=map.roiLocumY(j,1);
            g=sqrt((f-k)^2+(b-e/2)^2);
            beta=acosd((b-e/2)/g);
            
            %%% determines the pResponse for all cells %%%%
            [~,p]=max(mergedmeta.meanRepDG.ROIdata{1,ROIn}.summary.meanr);
            map.pdirInd(j,1)=p;
            map.Presponse{j,1}=mergedmeta.meanRepDB.ROIdata{1,ROIn}.summary.meannororddff(:,p);
            
            %%% calculate the onset and offset of specefic ROIs accounting for their location in the FOV
            
            switch p*45-45
                case 0
                    map.onset1(j,1)=(delay+t+b-c)/speedterm;
                    map.onset2(j,1)=(delay+t+b)/speedterm;
                case 45
                    x=g*sind(beta-45);
                    n=m-x;
                    map.onset1(j,1)=(delay+n-z-c)/speedterm;
                    map.onset2(j,1)=(delay+n-z)/speedterm;
                case 180
                    map.onset1(j,1)=(delay+t+(e-b)-c)/speedterm;
                    map.onset2(j,1)=(delay+t+(e-b))/speedterm;
                case 225
                    x=g*sind(beta-45);
                    n=m+x;
                    map.onset1(j,1)=(delay+n-z-c)/speedterm;
                    map.onset2(j,1)=(delay+n-z)/speedterm;
                case 90
                    map.onset1(j,1)=(delay+t+k-c)/speedterm;
                    map.onset2(j,1)=(delay+t+k)/speedterm;
                case 135
                    x=g*cosd(beta-45);
                    n=m-x;
                    map.onset1(j,1)=(delay+n-z-c)/speedterm;
                    map.onset2(j,1)=(delay+n-z)/speedterm;
                case 270
                    map.onset1(j,1)=(delay+t+(e-f)+(f-k)-c)/speedterm;
                    map.onset2(j,1)=(delay+t+(e-f)+(f-k))/speedterm;
                case 315
                    x=g*cosd(beta-45);
                    n=m+x;
                    map.onset1(j,1)=(delay+n-z-c)/speedterm;
                    map.onset2(j,1)=(delay+n-z)/speedterm;
            end
            j=j+1;
        end
    end
    
    XYSz=rotateVector([map.umroiX map.umroiY], deg2rad(-A));
    map.umroiXSz=XYSz(:,1);
    map.umroiYSz=XYSz(:,2);
    map.umroiXSznor=XYSz(:,1)./MR(1,1);
    map.umroiYSznor=XYSz(:,2)./MR(1,1);
    map.PXroiXSz=map.umroiXSz./pixelsize;
    map.PXroiYSz=map.umroiYSz./pixelsize;
    xctrans=xc+map.PXroiXSz;
    yctrans=yc-map.PXroiYSz;
    
    Rc=createRay(repmat([xc yc],size(xctrans,1),1),[xctrans  yctrans]);
    Ac=360-rad2deg(lineAngle(Rc));
    
    sector1=logical((Ac>A1) .* (Ac<A2));
    sector2=logical((Ac>A2) .* (Ac<A3));
    sector3=logical((Ac>A3) .* (Ac<A4));
    sector4=logical(logical(Ac>A4) + logical(Ac<A1));
    
    cellLoc(sector1)=1;
    cellLoc(sector2)=2;
    cellLoc(sector3)=3;
    cellLoc(sector4)=4;
    cellLoc=cellLoc';
    
    names=getFileList(filedir,'retrolabeledCell',0,'anywhere');
    if ~isempty(names)
        retro=importdata(cell2mat(names));
        map=processRetro(map,retro);
    end
    
    DS=map.isDS .* (map.DSI>0.17) .* (map.DSI<1);
    %DS=logical((map.isDS==1) .* (map.ONorONOFForOFF~=0));
    ONDS=logical((map.isDS==1) .* (map.ONorONOFForOFF==1));
    ONOFFDS=logical((map.isDS==1) .* (map.ONorONOFForOFF==2));
    OFFDS=logical((map.isDS==1) .* (map.ONorONOFForOFF==3));
    ONtDS=logical((map.isDS==1) .* (map.ONorONOFForOFF==4));    % ON transient DS cells
    ONOS=logical((map.isOS==1) .* (map.ONorONOFForOFF==1));
    ONOFFOS=logical((map.isOS==1) .* (map.ONorONOFForOFF==2));
    OFFOS=logical((map.isOS==1) .* (map.ONorONOFForOFF==3));
    ONsusOFFOS=logical((map.isOS==1) .* (map.ONorONOFForOFF==5));    % ONsusOFFOS ON sustained OFF OS cells
    
    if isfield(map,'isRetro')
        retroONDS=logical((map.isDS==1) .* (map.ONorONOFForOFF==1).*(map.isRetro==1));
        retroONOFFDS=logical((map.isDS==1) .* (map.ONorONOFForOFF==2).*(map.isRetro==1));
        retroOFFDS=logical((map.isDS==1) .* (map.ONorONOFForOFF==3).*(map.isRetro==1));
        retroOther=logical((map.isDS==0) .* (map.isOS==0) .*(map.isRetro==1));
        retroONtDS=logical((map.isDS==1) .* (map.ONorONOFForOFF==4).*(map.isRetro==1));
    end
    
    percentONDS=100*(sum(ONDS)/length(ONDS));
    percentONOFFDS=100*(sum(ONOFFDS)/length(ONOFFDS));
    percentOFFDS=100*(sum(OFFDS)/length(OFFDS));
    percentONOS=100*(sum(ONOS)/length(ONOS));
    percentONOFFOS=100*(sum(ONOFFOS)/length(ONOFFOS));
    percentOFFOS=100*(sum(OFFOS)/length(OFFOS));
    
    map.ONDS=ONDS;
    map.ONOFFDS=ONOFFDS;
    map.OFFDS=OFFDS;
    map.ONtDS=ONtDS;
    map.ONOS=ONOS;
    map.ONOFFOS=ONOFFOS;
    map.OFFOS=OFFOS;
    map.ONsusOFFOS=ONsusOFFOS;
    
    if isfield(map,'isRetro')
        map.retroONDS=retroONDS;
        map.retroONOFFDS=retroONOFFDS;
        map.retroOFFDS=retroOFFDS;
        map.retroOther=retroOther;
        try
            map.injectionSite=inject;
        end
    end
    
    %         try
    %         map.immuno1=cell2mat(immuno1);
    %         map.immuno2=cell2mat(immuno2);
    %         end
    
    map.percentONDS=percentONDS;
    map.percentONOFFDS=percentONOFFDS;
    map.percentOFFDS=percentOFFDS;
    map.percentONOS=percentONOS;
    map.percentONOFFOS=percentONOFFOS;
    map.percentOFFOS=percentOFFOS;
    
    if isfield(map,'isRetro')
        map.numRetroONDS=sum(retroONDS);
        map.numRetroONOFFDS=sum(retroONOFFDS);
        map.numRetroOFFDS=sum(retroOFFDS);
        map.numRetroOther=sum(retroOther);
    end  
    
    f2=copyobj(gcf,0);  %duplicate figure
    f3=copyobj(gcf,0);  %duplicate figure
    f4=copyobj(gcf,0);  %duplicate figure
    f7=copyobj(gcf,0);  %duplicate figure  
    
    % plot all cells
    set(0, 'CurrentFigure', f2) %make first figure active to plot all cells
    plot(xctrans,yctrans,'o','MarkerSize',2,'MarkerEdgeColor','r');
    
    % plot DS and OS cells and vectors
    
    %make first figure active to plot DS cells
    set(0, 'CurrentFigure', f3)
        
    % plot all DS cells in magenta FOR MANUSCRIPT
    %plot(xctrans(DS==1),yctrans(DS==1),'o','MarkerSize',2,'MarkerEdgeColor','m');
    %text(xctrans(DS),yctrans(DS),[repmat('  ',size(xctrans(DS)),1),num2str(find(DS==1))],'HorizontalAlignment','left','Color','r');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    plot(xctrans(ONDS),yctrans(ONDS),'o','MarkerSize',2,'MarkerEdgeColor','m');
    text(xctrans(ONDS),yctrans(ONDS),[repmat('  ',size(xctrans(ONDS)),1),num2str(find(ONDS==1))],'HorizontalAlignment','left','Color','r');
    %                 plot(xctrans(ONOFFDS),yctrans(ONOFFDS),'o','MarkerSize',2,'MarkerEdgeColor','m');
    %                 %text(xctrans(ONOFFDS),yctrans(ONOFFDS),[repmat('  ',size(xctrans(ONOFFDS)),1),num2str(find(ONOFFDS==1))],'HorizontalAlignment','left','Color','g');
    %                 plot(xctrans(OFFDS),yctrans(OFFDS),'o','MarkerSize',2,'MarkerEdgeColor','m');
    %                 %text(xctrans(OFFDS),yctrans(OFFDS),[repmat('  ',size(xctrans(OFFDS)),1),num2str(find(OFFDS==1))],'HorizontalAlignment','left','Color','b');
    %                 plot(xctrans(ONtDS),yctrans(ONtDS),'o','MarkerSize',2,'MarkerEdgeColor','m');
    %                 %text(xctrans(ONtDS),yctrans(ONtDS),[repmat('  ',size(xctrans(ONtDS)),1),num2str(find(ONtDS==1))],'HorizontalAlignment','left','Color','b');
    
    %standardize direction prefernces to follow matlab conventions
    %from (right=0, top=270, left=180,bottom=90) to (right=0, top=90, left=180,bottom=270)
    for i=1:size(map.pdirSz,1)
        if map.pdirSz(i,1)==0
            pdirSz(i,1)=0;
        else
            pdirSz(i,1)=360-map.pdirSz(i,1);
        end
    end
    
    map.pdirSz=pdirSz;
     
    %%%%%%%%
    %%%%%%%%%
    impdirSz=360-pdirSz;  %correct pdir to be used with quiver
    %%%%%%%%
    %%%%%%%%
    
    imadir=cosd(impdirSz(:,1));
    imbdir=sind(impdirSz(:,1));
    scale=0.2;
    
    % plot all DS cells in magenta FOR MANUSCRIPT
    %quiver(xctrans(DS==1),yctrans(DS==1),imadir(DS==1),imbdir(DS==1),scale,'m','linewidth',1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    quiver(xctrans(ONDS),yctrans(ONDS),imadir(ONDS),imbdir(ONDS),scale,'m','linewidth',1);
    quiver(xctrans(ONOFFDS),yctrans(ONOFFDS),imadir(ONOFFDS),imbdir(ONOFFDS),scale,'m','linewidth',1);
    quiver(xctrans(ONtDS),yctrans(ONtDS),imadir(ONtDS),imbdir(ONtDS),scale,'m','linewidth',1);
    quiver(xctrans(OFFDS),yctrans(OFFDS),imadir(OFFDS),imbdir(OFFDS),scale,'m','linewidth',2);
    %
    adir=cosd(pdirSz(:,1));
    bdir=sind(pdirSz(:,1));
    map.xpdirSz=adir;
    map.ypdirSz=bdir;
    map.SzMapOutputDS=[map.umroiXSznor,map.umroiYSznor,map.xpdirSz,map.ypdirSz];
    
    %make second figure active to plot OS cells
    set(0, 'CurrentFigure', f4)
    plot(xctrans(ONOS),yctrans(ONOS),'o','MarkerSize',3,'MarkerEdgeColor','r');
    %text(xctrans(ONOS),yctrans(ONOS),[repmat('  ',size(xctrans(ONOS)),1),num2str(find(ONOS==1))],'HorizontalAlignment','left','Color','r');
    plot(xctrans(ONOFFOS),yctrans(ONOFFOS),'o','MarkerSize',3,'MarkerEdgeColor','g');
    %text(xctrans(ONOFFOS),yctrans(ONOFFOS),[repmat('  ',size(xctrans(ONOFFOS)),1),num2str(find(ONOFFOS==1))],'HorizontalAlignment','left','Color','g');
    plot(xctrans(OFFOS),yctrans(OFFOS),'o','MarkerSize',3,'MarkerEdgeColor','b');
    %text(xctrans(OFFOS),yctrans(OFFOS),[repmat('  ',size(xctrans(OFFOS)),1),num2str(find(OFFOS==1))],'HorizontalAlignment','left','Color','b');
    plot(xctrans(ONsusOFFOS),yctrans(ONsusOFFOS),'o','MarkerSize',3,'MarkerEdgeColor','b');
    %text(xctrans(ONsusOFFOS),yctrans(ONsusOFFOS),[repmat('  ',size(xctrans(ONsusOFFOS)),1),num2str(find(ONsusOFFOS==1))],'HorizontalAlignment','left','Color','b');
    
    %standardize direction prefernces to follow matlab conventions
    %from (right=0, top=270, left=180,bottom=90) to (right=0, top=90, left=180,bottom=270)
    for i=1:size(map.poriSz,1)
        if map.poriSz(i,1)==0
            poriSz(i,1)=0;
        else
            poriSz(i,1)=180-map.poriSz(i,1);
        end
    end
    
    % converts the orientation in which the gratings move to the orientation of the gratings.
    for i=1:size(poriSz,1)
        if poriSz(i,1)>=90
            poriSz(i,1)=poriSz(i,1)-90;
        elseif poriSz(i,1)<90
            poriSz(i,1)=poriSz(i,1)+90;
        end
    end
    
    map.poriSz=poriSz;
    
    %%%%%%%%%%%
    imporiSz=360-poriSz;  %correct pori to be used with quiver (accounts for the definition of the direction in the drifting gratings and bars)
    %%%%%%%%%%%
    
    imaori=cosd(imporiSz(:,1));
    imbori=sind(imporiSz(:,1));
    scale=0.2;
    quiver(xctrans(ONOS),yctrans(ONOS),imaori(ONOS),imbori(ONOS),scale,'r','linewidth',2,'ShowArrowHead','off');
    quiver(xctrans(ONOFFOS),yctrans(ONOFFOS),imaori(ONOFFOS),imbori(ONOFFOS),scale,'g','linewidth',2,'ShowArrowHead','off');
    quiver(xctrans(OFFOS),yctrans(OFFOS),imaori(OFFOS),imbori(OFFOS),scale,'b','linewidth',2,'ShowArrowHead','off');
    quiver(xctrans(ONsusOFFOS),yctrans(ONsusOFFOS),imaori(ONsusOFFOS),imbori(ONsusOFFOS),scale,'b','linewidth',2,'ShowArrowHead','off');
    
    aori=cosd(poriSz(:,1));
    bori=sind(poriSz(:,1));
    map.xporiSz=aori;
    map.yporiSz=bori;
    map.SzMapOutputOS=[map.umroiXSznor,map.umroiYSznor,map.xporiSz,map.yporiSz];
    
    %make third figure active to plot DS retrolabeled cells
    if isfield(map,'isRetro')
        set(0, 'CurrentFigure', f7)
        hold on
        plot(xctrans(retroONDS),yctrans(retroONDS),'o','MarkerSize',2,'MarkerEdgeColor','r');
        text(xctrans(retroONDS),yctrans(retroONDS),[repmat('  ',size(xctrans(retroONDS)),1),num2str(find(retroONDS==1))],'HorizontalAlignment','left','Color','r');
        plot(xctrans(retroONOFFDS),yctrans(retroONOFFDS),'o','MarkerSize',2,'MarkerEdgeColor','g');
        text(xctrans(retroONOFFDS),yctrans(retroONOFFDS),[repmat('  ',size(xctrans(retroONOFFDS)),1),num2str(find(retroONOFFDS==1))],'HorizontalAlignment','left','Color','g');
        plot(xctrans(retroOFFDS),yctrans(retroOFFDS),'o','MarkerSize',2,'MarkerEdgeColor','b');
        text(xctrans(retroOFFDS),yctrans(retroOFFDS),[repmat('  ',size(xctrans(retroOFFDS)),1),num2str(find(retroOFFDS==1))],'HorizontalAlignment','left','Color','b');
        plot(xctrans(retroONtDS),yctrans(retroONtDS),'o','MarkerSize',2,'MarkerEdgeColor','b');
        text(xctrans(retroONtDS),yctrans(retroONtDS),[repmat('  ',size(xctrans(retroONtDS)),1),num2str(find(retroONtDS==1))],'HorizontalAlignment','left','Color','b');
        plot(xctrans(retroOther),yctrans(retroOther),'o','MarkerSize',2,'MarkerEdgeColor','m');
        text(xctrans(retroOther),yctrans(retroOther),[repmat('  ',size(xctrans(retroOther)),1),num2str(find(retroOther==1))],'HorizontalAlignment','left','Color','m');
        
        quiver(xctrans(retroONDS),yctrans(retroONDS),imadir(retroONDS),imbdir(retroONDS),scale,'r','linewidth',2);
        quiver(xctrans(retroONOFFDS),yctrans(retroONOFFDS),imadir(retroONOFFDS),imbdir(retroONOFFDS),scale,'g','linewidth',2);
        if sum(retroOFFDS); quiver(xctrans(retroOFFDS),yctrans(retroOFFDS),imadir(retroOFFDS),imbdir(retroOFFDS),scale,'b','linewidth',2); end
        if sum(retroOther);quiver(xctrans(retroOther),yctrans(retroOther),imadir(retroOther),imbdir(retroOther),scale,'m','linewidth',2); end
        hold off
    end
end

f5=figure;
set(f5,'position',[450 150 1500 500]);
set(0,'CurrentFigure', f5) %make figure active to plot all cells on polar plots
subplot(1,4,1)
[X,Y]=pol2cart(deg2rad(pdirSz(ONDS)), 1);
compassSS2p(X,Y,1,'r',1);
annotation('textbox', [0.19, 0.85, 0.1, 0.1], 'string','ON DS','LineStyle','none','FontSize',14);
subplot(1,4,2)
[X,Y]=pol2cart(deg2rad(pdirSz(ONOFFDS)), 1);
compassSS2p(X,Y,1,'g',1);
annotation('textbox', [0.38, 0.85, 0.1, 0.1], 'string','ON-OFF DS','LineStyle','none','FontSize',14);
subplot(1,4,3)
[X,Y]=pol2cart(deg2rad(pdirSz(OFFDS)), 1);
compassSS2p(X,Y,1,'b',1);
annotation('textbox', [0.59, 0.85, 0.1, 0.1], 'string','OFF DS','LineStyle','none','FontSize',14);
subplot(1,4,4)
[X,Y]=pol2cart(deg2rad(pdirSz(ONtDS)), 1);
compassSS2p(X,Y,1,'b',1);
annotation('textbox', [0.78, 0.85, 0.1, 0.1], 'string','ON trans DS','LineStyle','none','FontSize',14);
hold off

f6=figure;
set(f6,'position',[450 150 1500 500]);
set(0,'CurrentFigure', f6) %make figure active to plot all cells on polar plots
subplot(1,4,1)
[X,Y]=pol2cart(deg2rad(poriSz(ONOS)), 1);
compassSS2p(X,Y,1,'r',1);
annotation('textbox', [0.19, 0.85, 0.1, 0.1], 'string','ON OS','LineStyle','none','FontSize',14);
subplot(1,4,2)
[X,Y]=pol2cart(deg2rad(poriSz(ONOFFOS)), 1);
compassSS2p(X,Y,1,'g',1);
annotation('textbox', [0.38, 0.85, 0.1, 0.1], 'string','ON-OFF OS','LineStyle','none','FontSize',14);
subplot(1,4,3)
[X,Y]=pol2cart(deg2rad(poriSz(OFFOS)), 1);
compassSS2p(X,Y,1,'b',1);
annotation('textbox', [0.59, 0.85, 0.1, 0.1], 'string','OFF OS','LineStyle','none','FontSize',14);
subplot(1,4,4)
[X,Y]=pol2cart(deg2rad(poriSz(ONsusOFFOS)), 1);
compassSS2p(X,Y,1,'b',1);
annotation('textbox', [0.78, 0.85, 0.1, 0.1], 'string','ONsusOFF OS','LineStyle','none','FontSize',14);
hold off

if isfield(map,'isRetro')
    f8=figure;
    set(f8,'position',[450 150 1500 500]);
    set(0,'CurrentFigure', f8) %make figure active to plot all cells on polar plots
    subplot(1,4,1)
    [X,Y]=pol2cart(deg2rad(pdirSz(retroONDS)), 1);
    compassSS2p(X,Y,1,'r',1);
    annotation('textbox', [0.17, 0.85, 0.1, 0.1], 'string','retro ON DS','LineStyle','none','FontSize',14);
    subplot(1,4,2)
    [X,Y]=pol2cart(deg2rad(pdirSz(retroONOFFDS)), 1);
    compassSS2p(X,Y,1,'g',1);
    annotation('textbox', [0.36, 0.85, 0.1, 0.1], 'string','retro ON-OFF DS','LineStyle','none','FontSize',14);
    subplot(1,4,3)
    [X,Y]=pol2cart(deg2rad(pdirSz(retroOFFDS)), 1);
    compassSS2p(X,Y,1,'b',1);
    annotation('textbox', [0.58, 0.85, 0.1, 0.1], 'string','retro OFF DS','LineStyle','none','FontSize',14);
    subplot(1,4,4)
    [X,Y]=pol2cart(deg2rad(pdirSz(retroOther)), 1);
    compassSS2p(X,Y,1,'m',1);
    annotation('textbox', [0.79, 0.85, 0.1, 0.1], 'string','retro Other','LineStyle','none','FontSize',14);
    hold off
end

% plot subplots of pResponse of DS cells by type
try
    numPresponse=map.Presponse(map.ONDS);
    elements=[];
    for nu=1:length(map.Presponse(map.ONDS))
        elements=[elements,numel(numPresponse{nu})];
    end
    for nu=1:length(map.Presponse(map.ONDS))
        numPresponseT{nu}=numPresponse{nu}(1:min(elements));
    end
    ONDSpResponse=reshape(cell2mat(numPresponseT'),min(elements),[]);
    for p=1:size(ONDSpResponse,2)
        ONDSnorpResponse(:,p)=mat2gray(ONDSpResponse(:,p));
    end
    clear numPresponseT
    
    numPresponse=map.Presponse(map.ONOFFDS);
    elements=[];
    for nu=1:length(map.Presponse(map.ONOFFDS))
        elements=[elements,numel(numPresponse{nu})];
    end
    for nu=1:length(map.Presponse(map.ONOFFDS))
        numPresponseT{nu}=numPresponse{nu}(1:min(elements));
    end
    ONOFFDSpResponse=reshape(cell2mat(numPresponseT'),min(elements),[]);
    for p=1:size(ONOFFDSpResponse,2)
        ONOFFDSnorpResponse(:,p)=mat2gray(ONOFFDSpResponse(:,p));
    end
    clear numPresponseT
    
    numPresponse=map.Presponse(map.ONtDS);
    elements=[];
    for nu=1:length(map.Presponse(map.ONtDS))
        elements=[elements,numel(numPresponse{nu})];
    end
    for nu=1:length(map.Presponse(map.ONtDS))
        numPresponseT{nu}=numPresponse{nu}(1:min(elements));
        %numPresponseT{nu}=padarray(numPresponse{nu},[max(elements)-numel(numPresponse{nu}),0],NaN,'post');
    end
    ONtDSpResponse=reshape(cell2mat(numPresponseT'),min(elements),[]);
    for p=1:size(ONtDSpResponse,2)
        ONtDSnorpResponse(:,p)=mat2gray(ONtDSpResponse(:,p));
    end
    clear numPresponseT
end

f9=figure;
set(f9,'position',[5 150 950 800]);
subplot(2,3,1);
try plot(ONDSnorpResponse,'k'); end
xlim([0 size(ONDSnorpResponse,1)]);
annotation('textbox', [0.20, 0.87, 0.1, 0.1], 'string','ON DS','LineStyle','none','FontSize',14);
subplot(2,3,4);
try plot(mean(ONDSnorpResponse,2)); end
xlim([0 size(ONDSnorpResponse,1)]);
subplot(2,3,2);
try plot(ONOFFDSnorpResponse,'k'); end
xlim([0 size(ONDSnorpResponse,1)]);
annotation('textbox', [0.45, 0.87, 0.1, 0.1], 'string','ON-OFF DS','LineStyle','none','FontSize',14);
subplot(2,3,5);
try plot(mean(ONOFFDSnorpResponse,2)); end
xlim([0 size(ONDSnorpResponse,1)]);
subplot(2,3,3);
try plot(ONtDSnorpResponse,'k'); end
xlim([0 size(ONDSnorpResponse,1)]);
annotation('textbox', [0.7, 0.87, 0.1, 0.1], 'string','transient ON DS','LineStyle','none','FontSize',14);
subplot(2,3,6);
try plot(mean(ONtDSnorpResponse,2)); end
xlim([0 size(ONDSnorpResponse,1)]);

try
    % plot subplots of pResponse of OS cells by type
    numPresponse=map.Presponse(map.ONOS);
    elements=[];
    for nu=1:length(map.Presponse(map.ONOS))
        elements=[elements,numel(numPresponse{nu})];
    end
    for nu=1:length(map.Presponse(map.ONOS))
        numPresponseT{nu}=numPresponse{nu}(1:min(elements));
    end
    ONOSpResponse=reshape(cell2mat(numPresponseT'),min(elements),[]);
    for p=1:size(ONOSpResponse,2)
        ONOSnorpResponse(:,p)=mat2gray(ONOSpResponse(:,p));
    end
    clear numPresponseT
    
    numPresponse=map.Presponse(map.ONOFFOS);
    elements=[];
    for nu=1:length(map.Presponse(map.ONOFFOS))
        elements=[elements,numel(numPresponse{nu})];
    end
    for nu=1:length(map.Presponse(map.ONOFFOS))
        numPresponseT{nu}=numPresponse{nu}(1:min(elements));
    end
    ONOFFOSpResponse=reshape(cell2mat(numPresponseT'),min(elements),[]);
    for p=1:size(ONOFFOSpResponse,2)
        ONOFFOSnorpResponse(:,p)=mat2gray(ONOFFOSpResponse(:,p));
    end
    clear numPresponseT
    
    numPresponse=map.Presponse(map.ONsusOFFOS);
    elements=[];
    for nu=1:length(map.Presponse(map.ONsusOFFOS))
        elements=[elements,numel(numPresponse{nu})];
    end
    for nu=1:length(map.Presponse(map.ONsusOFFOS))
        numPresponseT{nu}=numPresponse{nu}(1:min(elements));
    end
    ONsusOFFOSpResponse=reshape(cell2mat(numPresponseT'),min(elements),[]);
    for p=1:size(ONsusOFFOSpResponse,2)
        ONsusOFFOSnorpResponse(:,p)=mat2gray(ONsusOFFOSpResponse(:,p));
    end
    clear numPresponseT
    
    f10=figure;
    set(f10,'position',[965 150 950 800]);
    subplot(2,3,1);
    plot(ONOSnorpResponse,'k');
    xlim([0 size(ONOSnorpResponse,1)]);
    annotation('textbox', [0.20, 0.87, 0.1, 0.1], 'string','ON OS','LineStyle','none','FontSize',14);
    subplot(2,3,4);
    plot(mean(ONOSnorpResponse,2));
    xlim([0 size(ONOSnorpResponse,1)]);
    subplot(2,3,2);
    plot(ONOFFOSnorpResponse,'k');
    xlim([0 size(ONOSnorpResponse,1)]);
    annotation('textbox', [0.45, 0.87, 0.1, 0.1], 'string','ON-OFF OS','LineStyle','none','FontSize',14);
    subplot(2,3,5);
    plot(mean(ONOFFOSnorpResponse,2));
    xlim([0 size(ONOSnorpResponse,1)]);
    subplot(2,3,3);
    plot(ONsusOFFOSnorpResponse,'k');
    xlim([0 size(ONOSnorpResponse,1)]);
    annotation('textbox', [0.7, 0.87, 0.1, 0.1], 'string','ON sustained OFF OS','LineStyle','none','FontSize',14);
    subplot(2,3,6);
    plot(mean(ONsusOFFOSnorpResponse,2));
    xlim([0 size(ONOSnorpResponse,1)]);
end

% tabledata for DS
tabledata.data=[cellLoc zeros(size(xctrans,1),1)];
tabledata.label={'cell position','NA'};
save('tabledata.mat', 'tabledata');
f=cellTypeSelectionGui2;
waitfor(f);
load tabledata
map.cellPosition=tabledata.data(:,1);
map.SzMapOutputDS(:,5)=map.cellPosition;
delete tabledata.mat

% tabledata for OS
tabledata.data=[cellLoc zeros(size(xctrans,1),1)];
tabledata.label={'cell position','NA'};
save('tabledata.mat', 'tabledata');
f=cellTypeSelectionGui2;
waitfor(f);
load tabledata
map.cellPosition=tabledata.data(:,1);
map.SzMapOutputOS(:,5)=map.cellPosition;
delete tabledata.mat

ind=strfind(filedir,'\');
filename=filedir;
filename(ind(end):length(filename))=[];
filename(1:ind(end-1))=[];
save(['map_',filename,'.mat'],'map');
saveas(f2,['map_allCells_',filename],'fig');
saveas(f3,['map_DS_',filename],'fig');
saveas(f4,['map_OS_',filename],'fig');
try savefig(f7,['map_retroDS_',filename]);end
try savefig(f9,['meanDScells_',filename]);end
try savefig(f10,['meanOScells_',filename]);end
saveas(f5,['polarPlot_allDSCells_',filename],'fig');
saveas(f6,['polarPlot_allOSCells_',filename],'fig');
try saveas(f8,['polarPlot_allRetroCells_',filename],'fig');end
javaaddpath(which('MatlabGarbageCollector.jar'))
clear all; clear java; jheapcl;

end


