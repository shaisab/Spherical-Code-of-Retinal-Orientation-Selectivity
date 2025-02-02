% open the appropriate pooledmap files with Itrans

function [] = clusterBinnedCellsModel()

[Name,~] = uigetfile('*.mat','Select a pooledmap file with color-codded vectors');
            load(Name);

for i=1:size(pooledmap.clusterdata,2)
    if pooledmap.clusterdata{1,i}.UVStvecComp.VecUSt
    [theta,rho]=cart2pol(pooledmap.clusterdata{1,i}.UVStvecComp.VecUSt,pooledmap.clusterdata{1,i}.UVStvecComp.VecVSt);
    ang=mod(theta,2*pi);
%     figure;
%     circ_plot(ang);
    meanAng(i,1:4)=[mod(circ_mean(ang(pooledmap.clusterdata{1,i}.Itrans==1)),2*pi),mod(circ_mean(ang(pooledmap.clusterdata{1,i}.Itrans==2)),2*pi),...
        mod(circ_mean(ang(pooledmap.clusterdata{1,i}.Itrans==3)),2*pi),mod(circ_mean(ang(pooledmap.clusterdata{1,i}.Itrans==4)),2*pi)]
%     hold on
%     circ_plot(meanAng(i,1:4),'r'); 
    meanAng=rad2deg(meanAng);
    angDiff(i)=abs(meanAng(i,1)-meanAng(i,2))
    if angDiff(i)>180
        angDiff(i)=360-angDiff(i);
    end
    else
    meanAng(i,1:4)=[nan,nan,nan,nan];
    angDiff(i)=nan;    
    end
   
end


angDiff=reshape(angDiff,10,[])';
angDiff=abs(flipud(angDiff))
figure
imagesc(angDiff)

% interpolate
%resolution=5;
xold=linspace(-1,1,10);
yold=linspace(-1,1,10);
[Aold,Bold]=meshgrid(xold,yold);

newRes=0.01;
xnew=-1:newRes:1;
ynew=-1:newRes:1;
[Anew,Bnew]=meshgrid(xnew,ynew);
data=griddata(Aold,Bold,angDiff,Anew,Bnew);
    

figure
contourf(xnew,ynew,flipud(data))
figure
imagesc(xnew,ynew,data)  
axis equal
axis off

c=colorbar('northoutside');
c.Label.String='angle difference (\circ)';
caxis([0 180]);
c.XTick=[0 45 90 135 180];
c.FontSize=14;
    
% [cid,ang,mu]=circ_clust(ang,4)
% hang=histogram(rad2deg(ang),360);
% valueHist=[hang.Values,0];
% periodogram(valueHist)
end


