

function []= poolSTHXYZ(discPoints,cellType,fieldType)


% The current directory should be 'All maps'
% should be called as:    poolSTHXYZ(cellType,fieldType)

[filedir,go]=bbdirselector('select data folder',cd);
if ~go disp('no folder selected'); return; end


%pooling and plotting vector fields in generalized (X,Y,Z) coordinates
namesMap=getFileList(filedir,['XYZvecComp_',num2str(discPoints),'_',cellType,'_',fieldType],0,'anywhere');
for i=1:size(namesMap,2)   
    if strfind(namesMap{i},'.mat')
   XYZvecCompAll(i)=load(num2str(namesMap{i}));
    end
end

X=[];
Y=[];
Z=[];
VecX=[];
VecY=[];
VecZ=[];
angleXYZ=[];
for i=1:size(XYZvecCompAll,2)
X=[X;XYZvecCompAll(1,i).X'];
Y=[Y;XYZvecCompAll(1,i).Y'];
Z=[Z;XYZvecCompAll(1,i).Z'];
VecX=[VecX;XYZvecCompAll(1,i).VecX'];
VecY=[VecY;XYZvecCompAll(1,i).VecY'];
VecZ=[VecZ;XYZvecCompAll(1,i).VecZ'];
angleXYZ=[angleXYZ;XYZvecCompAll(1,i).angleXYZ'];
end

pooledmap.XYZvecComp.X=X;
pooledmap.XYZvecComp.Y=Y;
pooledmap.XYZvecComp.Z=Z;
pooledmap.XYZvecComp.VecX=VecX;
pooledmap.XYZvecComp.VecY=VecY;
pooledmap.XYZvecComp.VecZ=VecZ;
pooledmap.XYZvecComp.angleXYZ=angleXYZ;

figure;
scale=0.2;
quiver3(X,Y,Z,VecX,VecY,VecZ,scale);

%pooling and plotting vector fields in generalized STH coordinates
namesMap=getFileList(filedir,['STHgenvecComp_',num2str(discPoints),'_',cellType,'_',fieldType],0,'anywhere');
for i=1:size(namesMap,2)   
    if strfind(namesMap{i},'.mat')
   STHgenvecCompAll(i)=load(num2str(namesMap{i}));
    end
end

Xprime=[];
Yprime=[];
VecXprime=[];
VecYprime=[];
angleSTHGen=[];
for i=1:size(STHgenvecCompAll,2)
Xprime=[Xprime;STHgenvecCompAll(1,i).Xprime'];
Yprime=[Yprime;STHgenvecCompAll(1,i).Yprime'];
VecXprime=[VecXprime;STHgenvecCompAll(1,i).VecXprime'];
VecYprime=[VecYprime;STHgenvecCompAll(1,i).VecYprime'];
angleSTHGen=[angleSTHGen;STHgenvecCompAll(1,i).angleSTHGen'];
end

pooledmap.STHgenvecComp.Xprime=Xprime;
pooledmap.STHgenvecComp.Yprime=Yprime;
pooledmap.STHgenvecComp.VecXprime=VecXprime;
pooledmap.STHgenvecComp.VecYprime=VecYprime;
pooledmap.STHgenvecComp.angleSTHGen=angleSTHGen;

figure;
scale=0.2;
quiver(Xprime,Yprime,VecXprime,VecYprime,scale);


%pooling and plotting vector fields in STH coordinates
namesMap=getFileList(filedir,['STHvecComp_',num2str(discPoints),'_',cellType,'_',fieldType],0,'anywhere');
for i=1:size(namesMap,2)   
    if strfind(namesMap{i},'.mat')
   STHvecCompAll(i)=load(num2str(namesMap{i}));
    end
end

S=[];
TH=[];
VecS=[];
VecTH=[];
angleSTH=[];
for i=1:size(STHvecCompAll,2)
S=[S;STHvecCompAll(1,i).S'];
TH=[TH;STHvecCompAll(1,i).TH'];
VecS=[VecS;STHvecCompAll(1,i).VecS'];
VecTH=[VecTH;STHvecCompAll(1,i).VecTH'];
angleSTH=[angleSTH;STHvecCompAll(1,i).angleSTH'];
end

pooledmap.STHvecComp.S=S;
pooledmap.STHvecComp.TH=TH;
pooledmap.STHvecComp.VecS=VecS;
pooledmap.STHvecComp.VecTH=VecTH;
pooledmap.STHvecComp.angleSTH=angleSTH;

figure;
scale=0.2;
quiver(S,TH/(2*pi),VecS,VecTH,scale);


%pooling and plotting vector fields in UV coordinates
namesMap=getFileList(filedir,['UVvecComp_',num2str(discPoints),'_',cellType,'_',fieldType],0,'anywhere');
for i=1:size(namesMap,2)   
    if strfind(namesMap{i},'.mat')
   UVvecCompAll(i)=load(num2str(namesMap{i}));
    end
end

U=[];
V=[];
VecU=[];
VecV=[];
angleUV=[];
for i=1:size(UVvecCompAll,2)
U=[U;UVvecCompAll(1,i).U'];
V=[V;UVvecCompAll(1,i).V'];
VecU=[VecU;UVvecCompAll(1,i).VecU'];
VecV=[VecV;UVvecCompAll(1,i).VecV'];
angleUV=[angleUV;UVvecCompAll(1,i).angleUV'];
end

pooledmap.UVvecComp.U=U;
pooledmap.UVvecComp.V=V;
pooledmap.UVvecComp.VecU=VecU;
pooledmap.UVvecComp.VecV=VecV;
pooledmap.UVvecComp.angleUV=angleUV;

figure;
scale=0.2;
quiver(U,V,VecU,VecV,scale);


%pooling and plotting vector fields in UVstandard coordinates
namesMap=getFileList(filedir,['UVStvecComp_',num2str(discPoints),'_',cellType,'_',fieldType],0,'anywhere');
for i=1:size(namesMap,2)   
    if strfind(namesMap{i},'.mat')
   UVStvecCompAll(i)=load(num2str(namesMap{i}));
    end
end

USt=[];
VSt=[];
VecUSt=[];
VecVSt=[];
angleUVSt=[];
for i=1:size(UVStvecCompAll,2)
USt=[USt;UVStvecCompAll(1,i).USt'];
VSt=[VSt;UVStvecCompAll(1,i).VSt'];
VecUSt=[VecUSt;UVStvecCompAll(1,i).VecUSt'];
VecVSt=[VecVSt;UVStvecCompAll(1,i).VecVSt'];
angleUVSt=[angleUVSt;UVStvecCompAll(1,i).angleUVSt'];
end

pooledmap.UVStvecComp.USt=USt;
pooledmap.UVStvecComp.VSt=VSt;
pooledmap.UVStvecComp.VecUSt=VecUSt;
pooledmap.UVStvecComp.VecVSt=VecVSt;
pooledmap.UVStvecComp.angleUVSt=angleUVSt;

%%%%%% plot standard retina surface %%%%%%%%
f1=figure;
StSec1Data=dlmread(['StandardSector1_',num2str(discPoints)]);
StSec2Data=dlmread(['StandardSector2_',num2str(discPoints)]);
StSec3Data=dlmread(['StandardSector3_',num2str(discPoints)]);
StSec4Data=dlmread(['StandardSector4_',num2str(discPoints)]);

n=size(StSec1Data,1);

StRHO1=StSec1Data(1:n,1:n);
StF1=StSec1Data(1:n,n+1:2*n);
StU1=StRHO1.*cos(StF1);
StV1=StRHO1.*sin(StF1);

StRHO2=StSec2Data(1:n,1:n);
StF2=StSec2Data(1:n,n+1:2*n);
StU2=StRHO2.*cos(StF2);
StV2=StRHO2.*sin(StF2);

StRHO3=StSec3Data(1:n,1:n);
StF3=StSec3Data(1:n,n+1:2*n);
StU3=StRHO3.*cos(StF3);
StV3=StRHO3.*sin(StF3);

StRHO4=StSec4Data(1:n,1:n);
StF4=StSec4Data(1:n,n+1:2*n);
StU4=StRHO4.*cos(StF4);
StV4=StRHO4.*sin(StF4);

hold on;
surf(StU1,StV1,zeros(size(StU1,1),size(StU1,1)), 'EdgeColor', 'none');
surf(StU2,StV2,zeros(size(StU2,1),size(StU2,1)), 'EdgeColor', 'none');
surf(StU3,StV3,zeros(size(StU3,1),size(StU3,1)), 'EdgeColor', 'none');
surf(StU4,StV4,zeros(size(StU4,1),size(StU4,1)), 'EdgeColor', 'none');
hold on;

pooledmap.standardRetina.StU1=StU1;
pooledmap.standardRetina.StU2=StU2;
pooledmap.standardRetina.StU3=StU3;
pooledmap.standardRetina.StU4=StU4;

f2=copyobj(gcf,0);  %duplicate figure f1 for polar plots overlay f1
set(f2,'position',[450 150 950 800]);
%f3=figure;  %for subplots


set(0,'CurrentFigure', f1) %make the UVstandard figure active to select a ROI 
scale=0.2;
quiver(USt,VSt,VecUSt,VecVSt,scale);
hold on;

% selecting a ROI

%if get(handles.addClusterbutton,'value')==1
nu=1;
while nu<5  
        set(0,'CurrentFigure', f1) %make the UVstandard figure active to select a ROI 
        set(f1,'position',[-1450 150 950 800]);
        h = impoly();
        nodes = getPosition(h);
        in=inpolygon(U,V,nodes(:,1),nodes(:,2));
        UStP=USt(in);
        VStP=VSt(in);
        k=convhull(UStP,VStP);
        
if isfield(pooledmap,'clusterdata')
    curClustern=length(pooledmap.clusterdata)+1;
else
    curClustern=1;
end
%set(handles.ROInTXT,'string',sprintf('current cluster: %d',curClustern), 'userdata',curClustern);

        set(0,'CurrentFigure', f2) %make the polar plots figure active
        [GEOM,~,~] =polygeom(UStP(k),VStP(k));
        %plot(UStP(k),VStP(k),'r-',UStP,VStP,'bo');
        %area(UStP(k),VStP(k));
        plot(UStP(k),VStP(k),'r-');
        scale=0.1;
        quiver(repmat(GEOM(1,2),size(VecUSt(in),1),1),repmat(GEOM(1,3),size(VecUSt(in),1),1),VecUSt(in),VecVSt(in),scale);
        text(UStP(VStP==min(VStP)),min(VStP)-0.1,num2str(curClustern));
%         set(0,'CurrentFigure', f3) %make the polar plots figure active
%         subplot(1,5,nu)
%         compassSS2p(VecUSt(in),VecVSt(in),1,'r',1);
        
end
%end

saveas(f1,['pooledMap_',num2str(discPoints),'_',cellType,'_',fieldType],'fig');
saveas(f2,['pooledMapClusters_',num2str(discPoints),'_',cellType,'_',fieldType],'fig');
%saveas(f3,['pooledMapPolar_',num2str(discPoints),'_',cellType,'_',fieldType],'fig');



end


