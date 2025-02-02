
function pooledmap = addCluster(pooledmap,f1,f2,f3)

set(0,'CurrentFigure', f1) %make the UVstandard figure active to select a ROI
%set(f1,'position',[-1450 150 950 800]);
set(f1,'position',[450 150 950 800]);
h = impoly();
nodes = getPosition(h);
in=inpolygon(pooledmap.UVStvecComp.USt,pooledmap.UVStvecComp.VSt,nodes(:,1),nodes(:,2));
plot([nodes(:,1);nodes(1,1)],[nodes(:,2);nodes(1,2)],'ro-','LineWidth',2);
UStP=pooledmap.UVStvecComp.USt(in);
VStP=pooledmap.UVStvecComp.VSt(in);
k=convhull(UStP,VStP);

if isfield(pooledmap,'clusterdata')
    curClustern=size(pooledmap.clusterdata,2)+1;
else
    curClustern=1;
end
    

hText=text(UStP(find(VStP==min(VStP),1)),min(VStP)-0.1,num2str(curClustern));

set(0,'CurrentFigure', f2) %make the polar plots figure active
t1=text(UStP(find(VStP==min(VStP),1)),min(VStP)-0.1,num2str(curClustern));
uistack(t1,'top');
[GEOM,~,~] =polygeom(UStP(k),VStP(k));
%plot(UStP(k),VStP(k),'r-',UStP,VStP,'bo');
%area(UStP(k),VStP(k));
plot(UStP(k),VStP(k),'r-','LineWidth',2);
scale=0.1;
quiver(repmat(GEOM(1,2),size(pooledmap.UVStvecComp.VecUSt(in),1),1),repmat(GEOM(1,3),size(pooledmap.UVStvecComp.VecUSt(in),1),1),pooledmap.UVStvecComp.VecUSt(in),pooledmap.UVStvecComp.VecVSt(in),scale,'k');
%quiver(repmat(GEOM(1,2),size(pooledmap.UVStvecComp.VecUSt(in),1),1),repmat(GEOM(1,3),size(pooledmap.UVStvecComp.VecUSt(in),1),1),pooledmap.UVStvecComp.VecUStcompASC(in),pooledmap.UVStvecComp.VecVStcompASC(in),scale,'b');
%quiver(repmat(GEOM(1,2),size(pooledmap.UVStvecComp.VecUSt(in),1),1),repmat(GEOM(1,3),size(pooledmap.UVStvecComp.VecUSt(in),1),1),pooledmap.UVStvecComp.VecUStcompPSC(in),pooledmap.UVStvecComp.VecVStcompPSC(in),scale,'g');
%quiver(repmat(GEOM(1,2),size(pooledmap.UVStvecComp.VecUSt(in),1),1),repmat(GEOM(1,3),size(pooledmap.UVStvecComp.VecUSt(in),1),1),pooledmap.UVStvecComp.VecUStcompLSC(in),pooledmap.UVStvecComp.VecVStcompLSC(in),scale,'r');
hold on



set(0,'CurrentFigure', f3) %make the polar plots figure active
[THETA,RHO]=cart2pol(pooledmap.UVStvecComp.VecUSt(in),pooledmap.UVStvecComp.VecVSt(in))
subplot(3,5,curClustern)
rose(THETA);

pooledmap.clusterdata{curClustern}.originX=repmat(GEOM(1,2),size(pooledmap.UVStvecComp.VecUSt(in),1),1);
pooledmap.clusterdata{curClustern}.originY=repmat(GEOM(1,3),size(pooledmap.UVStvecComp.VecUSt(in),1),1);
pooledmap.clusterdata{curClustern}.cellID=pooledmap.cellID(in);
pooledmap.clusterdata{curClustern}.typeClust=pooledmap.typeClust(in);
try pooledmap.clusterdata{curClustern}.Itrans=pooledmap.Itrans(in); end
pooledmap.clusterdata{curClustern}.UVStvecComp.USt=pooledmap.UVStvecComp.USt(in);
pooledmap.clusterdata{curClustern}.UVStvecComp.VSt=pooledmap.UVStvecComp.VSt(in);
pooledmap.clusterdata{curClustern}.UVStvecComp.VecUSt=pooledmap.UVStvecComp.VecUSt(in);
pooledmap.clusterdata{curClustern}.UVStvecComp.VecVSt=pooledmap.UVStvecComp.VecVSt(in);
pooledmap.clusterdata{curClustern}.UVStvecComp.VecUStcompASC=pooledmap.UVStvecComp.VecUStcompASC(in);
pooledmap.clusterdata{curClustern}.UVStvecComp.VecVStcompASC=pooledmap.UVStvecComp.VecVStcompASC(in);
pooledmap.clusterdata{curClustern}.UVStvecComp.angleUVStASC=pooledmap.UVStvecComp.angleUVStASC(in);
pooledmap.clusterdata{curClustern}.UVStvecComp.VecUStcompPSC=pooledmap.UVStvecComp.VecUStcompPSC(in);
pooledmap.clusterdata{curClustern}.UVStvecComp.VecVStcompPSC=pooledmap.UVStvecComp.VecVStcompPSC(in);
pooledmap.clusterdata{curClustern}.UVStvecComp.angleUVStPSC=pooledmap.UVStvecComp.angleUVStPSC(in);
pooledmap.clusterdata{curClustern}.UVStvecComp.VecUStcompLSC=pooledmap.UVStvecComp.VecUStcompLSC(in);
pooledmap.clusterdata{curClustern}.UVStvecComp.VecVStcompLSC=pooledmap.UVStvecComp.VecVStcompLSC(in);
pooledmap.clusterdata{curClustern}.UVStvecComp.angleUVStLSC=pooledmap.UVStvecComp.angleUVStLSC(in);

pooledmap.clusterdata{curClustern}.XYZvecComp.X=pooledmap.XYZvecComp.X(in);
pooledmap.clusterdata{curClustern}.XYZvecComp.Y=pooledmap.XYZvecComp.Y(in);
pooledmap.clusterdata{curClustern}.XYZvecComp.Z=pooledmap.XYZvecComp.Z(in);
pooledmap.clusterdata{curClustern}.XYZvecComp.VecX=pooledmap.XYZvecComp.VecX(in);
pooledmap.clusterdata{curClustern}.XYZvecComp.VecY=pooledmap.XYZvecComp.VecY(in);
pooledmap.clusterdata{curClustern}.XYZvecComp.VecZ=pooledmap.XYZvecComp.VecZ(in);
pooledmap.clusterdata{curClustern}.XYZvecComp.VecXcompASC=pooledmap.XYZvecComp.VecXcompASC(in);
pooledmap.clusterdata{curClustern}.XYZvecComp.VecYcompASC=pooledmap.XYZvecComp.VecYcompASC(in);
pooledmap.clusterdata{curClustern}.XYZvecComp.VecZcompASC=pooledmap.XYZvecComp.VecZcompASC(in);
pooledmap.clusterdata{curClustern}.XYZvecComp.angleXYZASC=pooledmap.XYZvecComp.angleXYZASC(in);
pooledmap.clusterdata{curClustern}.XYZvecComp.VecXcompPSC=pooledmap.XYZvecComp.VecXcompPSC(in);
pooledmap.clusterdata{curClustern}.XYZvecComp.VecYcompPSC=pooledmap.XYZvecComp.VecYcompPSC(in);
pooledmap.clusterdata{curClustern}.XYZvecComp.VecZcompPSC=pooledmap.XYZvecComp.VecZcompPSC(in);
pooledmap.clusterdata{curClustern}.XYZvecComp.angleXYZPSC=pooledmap.XYZvecComp.angleXYZPSC(in);
pooledmap.clusterdata{curClustern}.XYZvecComp.VecXcompLSC=pooledmap.XYZvecComp.VecXcompLSC(in);
pooledmap.clusterdata{curClustern}.XYZvecComp.VecYcompLSC=pooledmap.XYZvecComp.VecYcompLSC(in);
pooledmap.clusterdata{curClustern}.XYZvecComp.VecZcompLSC=pooledmap.XYZvecComp.VecZcompLSC(in);
pooledmap.clusterdata{curClustern}.XYZvecComp.angleXYZLSC=pooledmap.XYZvecComp.angleXYZLSC(in);



end


