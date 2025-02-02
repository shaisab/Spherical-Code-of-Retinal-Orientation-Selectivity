
function pooledmap = addClusterGrid(pooledmap,COLORstate,f1,f2,f3)

set(0,'CurrentFigure', f1) %make the UVstandard figure active to select a ROI
%set(f1,'position',[-1450 150 950 800]);
set(f1,'position',[450 150 950 800]);

[X,Y] = meshgrid(-1:.2:1, -1:.2:1);
% X=reshape(X,[numel(X),1]);
% Y=reshape(Y,[numel(Y),1]);

for g=1:size(X,1)-1
    for n=1:size(X,2)-1
nodes=[X(g,n),Y(g,n);X(g,n+1),Y(g,n+1);X(g+1,n),Y(g+1,n);X(g+1,n+1),Y(g+1,n+1)];

set(0,'CurrentFigure', f2) %make the polar plots figure active
plot([nodes(1,1),nodes(2,1),nodes(4,1),nodes(3,1)],[nodes(1,2),nodes(2,2),nodes(4,2),nodes(3,2)],'r-');


in=inpolygon(pooledmap.UVvecComp.U,pooledmap.UVvecComp.V,nodes(:,1),nodes(:,2));




UStP=pooledmap.UVStvecComp.USt(in);
VStP=pooledmap.UVStvecComp.VSt(in);

if size(VStP,1)>2 && numel(unique(VStP))>2 && numel(unique(UStP))>2
k=convhull(UStP,VStP);

if isfield(pooledmap,'clusterdata')
    curClustern=size(pooledmap.clusterdata,2)+1;
else
    curClustern=1;
end

%hText=text(UStP(find(VStP==min(VStP),1)),min(VStP)-0.05,num2str(curClustern));


%t1=text(UStP(find(VStP==min(VStP),1)),min(VStP)-0,num2str(curClustern));

[GEOM,~,~] =polygeom(UStP(k),VStP(k));

t1=text(nodes(2,1)-0.04,nodes(2,2)+0.03,num2str(curClustern));
%t1=text(GEOM(1,2),GEOM(1,3),num2str(curClustern));
uistack(t1,'top');
%plot(UStP(k),VStP(k),'r-',UStP,VStP,'bo');
%area(UStP(k),VStP(k));
%plot(UStP(k),VStP(k),'r-');
scale=0.1;

switch COLORstate
    case 0
ind=find(in);
for c=1:size(ind,1)
    if pooledmap.Itrans(ind(c))==1
        quiver(repmat(GEOM(1,2),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),repmat(GEOM(1,3),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),pooledmap.UVStvecComp.VecUSt(ind(c)),pooledmap.UVStvecComp.VecVSt(ind(c)),scale,'b');
    elseif pooledmap.Itrans(ind(c))==2
        quiver(repmat(GEOM(1,2),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),repmat(GEOM(1,3),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),pooledmap.UVStvecComp.VecUSt(ind(c)),pooledmap.UVStvecComp.VecVSt(ind(c)),scale,'g');
    elseif pooledmap.Itrans(ind(c))==3
        quiver(repmat(GEOM(1,2),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),repmat(GEOM(1,3),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),pooledmap.UVStvecComp.VecUSt(ind(c)),pooledmap.UVStvecComp.VecVSt(ind(c)),scale,'r');
    elseif pooledmap.Itrans(ind(c))==4
        quiver(repmat(GEOM(1,2),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),repmat(GEOM(1,3),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),pooledmap.UVStvecComp.VecUSt(ind(c)),pooledmap.UVStvecComp.VecVSt(ind(c)),scale,'k');
    end
end   

    case 3
     quiver(repmat(GEOM(1,2),size(pooledmap.UVStvecComp.VecUSt(in),1),1),repmat(GEOM(1,3),size(pooledmap.UVStvecComp.VecUSt(in),1),1),pooledmap.UVStvecComp.VecUSt(in),pooledmap.UVStvecComp.VecVSt(in),scale,'k');   
end
hold on



set(0,'CurrentFigure', f3) %make the polar plots figure active
[THETA,RHO]=cart2pol(pooledmap.UVStvecComp.VecUSt(in),pooledmap.UVStvecComp.VecVSt(in))
subplot(floor(size(X,1)),floor(size(X,2)),curClustern)
rose(THETA);

pooledmap.clusterdata{curClustern}.originX=repmat(GEOM(1,2),size(pooledmap.UVStvecComp.VecUSt(in),1),1);
pooledmap.clusterdata{curClustern}.originY=repmat(GEOM(1,3),size(pooledmap.UVStvecComp.VecUSt(in),1),1);
pooledmap.clusterdata{curClustern}.cellID=pooledmap.cellID(in);
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
    end
end
end


