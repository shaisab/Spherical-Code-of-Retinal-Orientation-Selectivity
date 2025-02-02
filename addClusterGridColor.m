
function pooledmap = addClusterGridColor(pooledmap,COLORstate,f1,f2,f3,f4)

% % create spherical retina figure
%         f4=figure
%         f4=RetinaPlot(f4);
        

set(0,'CurrentFigure', f1) %make the UVstandard figure active to select a ROI
set(f1,'position',[450 150 950 800]);

[X,Y] = meshgrid(-1:.2:1, -1:.2:1);

for g=1:size(X,1)-1
    for n=1:size(X,2)-1
        nodes=[X(g,n),Y(g,n);X(g,n+1),Y(g,n+1);X(g+1,n+1),Y(g+1,n+1);X(g+1,n),Y(g+1,n)];
        
        
        
        set(0,'CurrentFigure', f2) %make the polar plots figure active
        %plot([nodes(1,1),nodes(2,1),nodes(3,1),nodes(4,1)],[nodes(1,2),nodes(2,2),nodes(3,2),nodes(4,2)],'m--','LineWidth',0.5);
        
        
        [in on]=inpolygon(pooledmap.UVStvecComp.USt,pooledmap.UVStvecComp.VSt,nodes(:,1),nodes(:,2));
%         [~,ind_in,~]=intersect(in,on);
%         in(ind_in)=0;

        [GEOM,~,~] =polygeom(nodes(:,1),nodes(:,2));

        [Uout(g,n),Vout(g,n)] = StandardRemove(GEOM(2),GEOM(3));
        [x(g,n),y(g,n),z(g,n)]=Standard2Retina(Uout(g,n),Vout(g,n));     % these are the x,y,z coordiantes of the centroid of the bins
        

        UStP=pooledmap.UVStvecComp.USt(in);
        VStP=pooledmap.UVStvecComp.VSt(in);
        
        % make vectors tangent to spherical retina
%         phys=dlmread('phys');
%         R=phys(1);
%         ind=find(in);
%         V=[pooledmap.XYZvecComp.VecX(ind),pooledmap.XYZvecComp.VecY(ind),pooledmap.XYZvecComp.VecZ(ind)];
%         for Vn=1:size(V,1)
%         [VT(Vn,:)]=TangentConversion(x(g,n),y(g,n),z(g,n),R,V(Vn,:));
%         end
        
        
        %if size(VStP,1)>0
        
        if isfield(pooledmap,'clusterdata')
            curClustern=size(pooledmap.clusterdata,2)+1;
        else
            curClustern=1;
        end
        
%         if curClustern==100
%         t1=text(nodes(2,1)-0.065,nodes(2,2)+0.03,num2str(curClustern));
%         else
%         t1=text(nodes(2,1)-0.045,nodes(2,2)+0.03,num2str(curClustern));    
%         end


        scale=0.08;
        %scale=0.07;
        %scale=0.09;
        %scale=0.08;
        
        switch COLORstate
            case 0
                c1=0;
                c2=0;
                c3=0;
                c4=0;
                ind=find(in);
                for c=1:size(ind,1)
                    if pooledmap.Itrans(ind(c))==1
                        c1=c1+1;
                        % standard retina
                        set(0,'CurrentFigure', f2) %make the polar plots figure active
                        quiver(repmat(GEOM(2),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),repmat(GEOM(3),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),pooledmap.UVStvecComp.VecUSt(ind(c)),pooledmap.UVStvecComp.VecVSt(ind(c)),scale,'b','ShowArrowHead','off');
                        % sphrtical retina
                        set(0,'CurrentFigure', f4) %make the spherical retina figure active
                        quiver3(repmat(x(g,n),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),repmat(y(g,n),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),repmat(z(g,n),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),...
                        pooledmap.XYZvecComp.VecX(ind(c)),pooledmap.XYZvecComp.VecY(ind(c)),pooledmap.XYZvecComp.VecZ(ind(c))+0.3,scale,'k',ShowArrowHead,'off');
                    elseif pooledmap.Itrans(ind(c))==2
                        c2=c2+1;
                        % standard retina
                        set(0,'CurrentFigure', f2) %make the polar plots figure active
                        quiver(repmat(GEOM(2),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),repmat(GEOM(3),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),pooledmap.UVStvecComp.VecUSt(ind(c)),pooledmap.UVStvecComp.VecVSt(ind(c)),scale,'g','ShowArrowHead','off');
                        % sphrtical retina
                        set(0,'CurrentFigure', f4) %make the spherical retina figure active
                        quiver3(repmat(x(g,n),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),repmat(y(g,n),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),repmat(z(g,n),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),...
                        pooledmap.XYZvecComp.VecX(ind(c)),pooledmap.XYZvecComp.VecY(ind(c)),pooledmap.XYZvecComp.VecZ(ind(c))+0.3,scale,'k','ShowArrowHead','off');
                    elseif pooledmap.Itrans(ind(c))==3
                        c3=c3+1;
                        % standard retina
                        set(0,'CurrentFigure', f2) %make the polar plots figure active
                        quiver(repmat(GEOM(2),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),repmat(GEOM(3),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),pooledmap.UVStvecComp.VecUSt(ind(c)),pooledmap.UVStvecComp.VecVSt(ind(c)),scale,'r','ShowArrowHead','off');
                        % sphrtical retina
                        set(0,'CurrentFigure', f4) %make the spherical retina figure active
                        quiver3(repmat(x(g,n),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),repmat(y(g,n),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),repmat(z(g,n),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),...
                        pooledmap.XYZvecComp.VecX(ind(c)),pooledmap.XYZvecComp.VecY(ind(c)),pooledmap.XYZvecComp.VecZ(ind(c))+0.3,scale,'k','ShowArrowHead','off');
                    elseif pooledmap.Itrans(ind(c))==4
                        c4=c4+1;
                        % standard retina
                        set(0,'CurrentFigure', f2) %make the polar plots figure active
                        quiver(repmat(GEOM(2),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),repmat(GEOM(3),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),pooledmap.UVStvecComp.VecUSt(ind(c)),pooledmap.UVStvecComp.VecVSt(ind(c)),scale,'k','ShowArrowHead','off');
                        % sphrtical retina
                        set(0,'CurrentFigure', f4) %make the spherical retina figure active
                        quiver3(repmat(x(g,n),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),repmat(y(g,n),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),repmat(z(g,n),size(pooledmap.UVStvecComp.VecUSt(ind(c)),1),1),...
                        pooledmap.XYZvecComp.VecX(ind(c)),pooledmap.XYZvecComp.VecY(ind(c)),pooledmap.XYZvecComp.VecZ(ind(c))+0.3,scale,'r','ShowArrowHead','off');
                    end
                end
                try pooledmap.cellCount{g,n}=[c1,c2,c3,c4]; end
            case 3
                ind=find(in);
                % standard retina       GOLD [1,0.6314,0]
                quiver(repmat(GEOM(2),size(pooledmap.UVStvecComp.VecUSt(ind),1),1),repmat(GEOM(3),size(pooledmap.UVStvecComp.VecUSt(ind),1),1),pooledmap.UVStvecComp.VecUSt(ind),pooledmap.UVStvecComp.VecVSt(ind),scale,'Color',[0,0,0],'ShowArrowHead','off');
                % sphrtical retina
                set(0,'CurrentFigure', f4) %make the spherical retina figure active
                quiver3(repmat(x(g,n),size(pooledmap.UVStvecComp.VecUSt(ind),1),1),repmat(y(g,n),size(pooledmap.UVStvecComp.VecUSt(ind),1),1),repmat(z(g,n),size(pooledmap.UVStvecComp.VecUSt(ind),1),1),...
                pooledmap.XYZvecComp.VecX(ind),pooledmap.XYZvecComp.VecY(ind),pooledmap.XYZvecComp.VecZ(ind)+0.3,scale,'k','ShowArrowHead','off');
%                 if exist('VT','var')
%                     if ~isnan(VT(1)) && ~isnan(x(g,n)) && size(pooledmap.UVStvecComp.VecUSt(ind),1)
%                 quiver3(repmat(x(g,n),size(pooledmap.UVStvecComp.VecUSt(ind),1),1),repmat(y(g,n),size(pooledmap.UVStvecComp.VecUSt(ind),1),1),repmat(z(g,n),size(pooledmap.UVStvecComp.VecUSt(ind),1),1),...
%                                 VT(:,1),VT(:,2),VT(:,3),scale,'k'); 
%                     end
%                 end
        end
        %axis equal
        xlim([-1 1]);ylim([-1 1]);
        hold on
        
        
        
        set(0,'CurrentFigure', f3) %make the polar plots figure active
        [THETA,RHO]=cart2pol(pooledmap.UVStvecComp.VecUSt(in),pooledmap.UVStvecComp.VecVSt(in))
        subplot(floor(size(X,1)),floor(size(X,2)),curClustern)
        rose(THETA);
        
        % pooledmap.clusterdata{curClustern}.originX=repmat(GEOM(1,2),size(pooledmap.UVStvecComp.VecUSt(in),1),1);
        % pooledmap.clusterdata{curClustern}.originY=repmat(GEOM(1,3),size(pooledmap.UVStvecComp.VecUSt(in),1),1);
        pooledmap.clusterdata{curClustern}.cellID=pooledmap.cellID(in);
        pooledmap.clusterdata{curClustern}.poriSz=pooledmap.poriSz(in);
        Orientation_alignment=abs(asind(sind(pooledmap.clusterdata{curClustern}.poriSz)));
        pooledmap.clusterdata{curClustern}.Orientation_alignment_index=sum(Orientation_alignment>60|Orientation_alignment<30)/length(Orientation_alignment);
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
        
        
        %end
    end
end
hold off
end


