
function []= mapHb9Cells()

%% load and process retinal image
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
if strfind(filedir,'Aug_22_2016')
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

%% load xy position and Hb9 position files
[Hb9xyposFileList]=getFileList(filedir,'Hb9position.txt',0,'anywhere');
Hb9xypos=tdfread(Hb9xyposFileList{1});

[xyposFileList]=getFileList(filedir,'xyposition.txt',0,'anywhere');
xypos=tdfread(xyposFileList{1});

[umroiX,umroiY]=roi_locHb9(xypos,Hb9xypos);

%% process the position and preferred direction of cells
if size(Hb9xyposFileList,2)~=0 && size(xyposFileList,2)~=0
        
    map.umroiX=umroiX;
    map.umroiY=umroiY;
    
    for ROIn=1:size(Hb9xypos.angle,1)
        %pdir=Hb9xypos.angle(ROIn);
        pdir=mod(Hb9xypos.angle(ROIn),360);
        pdirSz=pdir+A;
        if pdirSz<0
            pdirSz=360+pdirSz;
        end  
        
        map.cellID{ROIn,1}=([filedir(8:18),'_fov#',mat2str(Hb9xypos.fov(ROIn)),'_cell#',mat2str(Hb9xypos.cell(ROIn))]);
        map.DSI(ROIn,1)=Hb9xypos.distance(ROIn)./Hb9xypos.perimeter(ROIn);
        map.pdir(ROIn,1)=pdir;
        map.pdirSz(ROIn,1)=pdirSz;
        map.isDS(ROIn,1)=1; 
        map.injectionSite{ROIn,1}='';       % for compatibility with downstream functions
        map.isCART(ROIn,1)=0;               % for compatibility with downstream functions
        map.isRBPMS(ROIn,1)=0;              % for compatibility with downstream functions
        map.isRetro(ROIn,1)=0;              % for compatibility with downstream functions
        map.ONorONOFForOFF(ROIn,1)=2;       % for compatibility with downstream functions
        map.Presponse{ROIn,1}=ones(100,1);  % for compatibility with downstream functions
        map.onset1(ROIn,1)=1;               % for compatibility with downstream functions
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
    
end

%% plot figures
f3=copyobj(gcf,0);  %duplicate figure

% plot DS cells and vectors
set(0, 'CurrentFigure', f3) 
plot(xctrans,yctrans,'o','MarkerSize',4,'MarkerEdgeColor','r');
text(xctrans,yctrans,num2str([1:size(xctrans,1)]'),'HorizontalAlignment','left','Color','r','FontSize',10);

%standardize direction prefernces to follow matlab conventions
%from (right=0, top=270, left=180,bottom=90) to (right=0, top=90, left=180,bottom=270)
% for i=1:size(map.pdirSz,1)
%     if map.pdirSz(i,1)==0
%         pdirSz(i,1)=0;
%     else
%         pdirSz(i,1)=360-map.pdirSz(i,1);
%     end
% end

pdirSz=map.pdirSz;


%%%%%%%%
%%%%%%%%%
impdirSz=360-pdirSz;  %correct pdir to be used with quiver
%%%%%%%%
%%%%%%%%


imadir=cosd(impdirSz(:,1));
imbdir=sind(impdirSz(:,1));
scale=0.2;
quiver(xctrans,yctrans,imadir,imbdir,scale,'r','linewidth',1);

adir=cosd(pdirSz(:,1));
bdir=sind(pdirSz(:,1));
map.xpdirSz=adir;
map.ypdirSz=bdir;
map.SzMapOutputDS=[map.umroiXSznor,map.umroiYSznor,map.xpdirSz,map.ypdirSz];

% polar plot
f5=figure;
%set(f5,'position',[450 150 1500 500]);
set(0,'CurrentFigure', f5) %make figure active to plot all cells on polar plots 
[X,Y]=pol2cart(deg2rad(pdirSz), 1);
compassSS2p(X,Y,1,'r',1);
annotation('textbox', [0.05, 0.85, 0.1, 0.1], 'string','Hb9 ON-OFF DS','LineStyle','none','FontSize',14);

%% tabledata for DS
tabledata.data=[cellLoc zeros(size(xctrans,1),1)];
tabledata.label={'cell position','NA'};
save('tabledata.mat', 'tabledata');
f=cellTypeSelectionGui2;
waitfor(f);
load tabledata
map.cellPosition=tabledata.data(:,1);
map.SzMapOutputDS(:,5)=map.cellPosition;
delete tabledata.mat

%% save file and figures
ind=strfind(filedir,'\');
filename=filedir;
filename(ind(end):length(filename))=[];
filename(1:ind(end-1))=[];
save(['map_',filename,'.mat'],'map');
saveas(f3,['map_DS_',filename],'fig');
saveas(f5,['polarPlot_allDSCells_',filename],'fig');
javaaddpath(which('MatlabGarbageCollector.jar'))
clear all; clear java; jheapcl;

end


