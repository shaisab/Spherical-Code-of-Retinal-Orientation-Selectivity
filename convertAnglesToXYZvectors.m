
function ang = convertAnglesToXYZvectors(data)

phys=dlmread('phys');       % get R from the phys file
R=phys(1);

numerator=-data.XYZvecComp.Y.*data.XYZvecComp.VecX+data.XYZvecComp.X.*data.XYZvecComp.VecY;
denominator=sqrt(2.*data.XYZvecComp.Z.*R-data.XYZvecComp.Z.^2);

%convert the x,y,z components to angles in degrees between 0 and 360
ang=real(acosd(numerator./denominator));
figure
circ_plot(deg2rad(ang));

cond=data.XYZvecComp.VecX.*((data.XYZvecComp.X.*(R-data.XYZvecComp.Z))./sqrt(data.XYZvecComp.X.^2+data.XYZvecComp.Y.^2))+...
    data.XYZvecComp.VecY.*((data.XYZvecComp.Y.*(R-data.XYZvecComp.Z))./sqrt(data.XYZvecComp.X.^2+data.XYZvecComp.Y.^2))+...
        data.XYZvecComp.VecZ.*sqrt((data.XYZvecComp.X.^2+data.XYZvecComp.Y.^2)./R.^2);

for j=1:size(cond,1)
    if cond(j)>=0
        ang(j)=ang(j);
    elseif cond(j)<0
        ang(j)=360-ang(j);
    end
end

%% plot quiver polar plot for cluster 1,VecUSt and VecVSt % change clusern
clustern=35;
figure; quiver(repmat(0,size(pooledmap.clusterdata{1, clustern}.UVStvecComp.VecUSt,1),1),repmat(0,size(pooledmap.clusterdata{1, clustern}.UVStvecComp.VecUSt,1),1),...
    pooledmap.clusterdata{1, clustern}.UVStvecComp.VecUSt,pooledmap.clusterdata{1, clustern}.UVStvecComp.VecVSt ,1,'k','LineWidth',1.5);
set(findobj(gcf, 'type','axes'), 'Visible','off')
set(gcf, 'color', [1 1 1]);
%axis[-1 1 -1 1]
axis equal;
hold off





end