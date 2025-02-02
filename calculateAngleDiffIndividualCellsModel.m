% open the appropriate pooledmap files with Itrans

function [] = clusterBinnedCellsModel()

%% load file and identify the four vectors originating from each retinal location
[Name,~] = uigetfile('*.mat','Select a pooledmap file with color-codded vectors');
            load(Name);

j=1;
for i=1:size(pooledmap.UVStvecComp.USt,1)
    ind=intersect(find(pooledmap.UVStvecComp.USt==pooledmap.UVStvecComp.USt(i)),find(pooledmap.UVStvecComp.VSt==pooledmap.UVStvecComp.VSt(i)));
    if ~isempty(ind)
        for k=1:size(ind,1)
            if pooledmap.Itrans(ind(k))==1
                PU(j,1)=pooledmap.UVStvecComp.USt(ind(k));
                PV(j,1)=pooledmap.UVStvecComp.VSt(ind(k));
                PvecU(j,1)=pooledmap.UVStvecComp.VecUSt(ind(k));
                PvecV(j,1)=pooledmap.UVStvecComp.VecVSt(ind(k));
            elseif pooledmap.Itrans(ind(k))==2
                PU(j,2)=pooledmap.UVStvecComp.USt(ind(k));
                PV(j,2)=pooledmap.UVStvecComp.VSt(ind(k));
                PvecU(j,2)=pooledmap.UVStvecComp.VecUSt(ind(k));
                PvecV(j,2)=pooledmap.UVStvecComp.VecVSt(ind(k));
            elseif pooledmap.Itrans(ind(k))==3
                PU(j,3)=pooledmap.UVStvecComp.USt(ind(k));
                PV(j,3)=pooledmap.UVStvecComp.VSt(ind(k));
                PvecU(j,3)=pooledmap.UVStvecComp.VecUSt(ind(k));
                PvecV(j,3)=pooledmap.UVStvecComp.VecVSt(ind(k));
            elseif pooledmap.Itrans(ind(k))==4
                PU(j,4)=pooledmap.UVStvecComp.USt(ind(k));
                PV(j,4)=pooledmap.UVStvecComp.VSt(ind(k));
                PvecU(j,4)=pooledmap.UVStvecComp.VecUSt(ind(k));
                PvecV(j,4)=pooledmap.UVStvecComp.VecVSt(ind(k));  
            end 
        end
        j=j+1;
    end
end            
            
            
%% calculate angle difference for at each retina location                     

[theta,rho]=cart2pol(PvecU,PvecV);
ang=mod(theta,2*pi);
ang=rad2deg(ang);

for i=1:size(PU,1)
%     hold on
%     circ_plot(ang(i,1:4),'r');    
    angDiff(i,1)=abs(ang(i,1)-ang(i,2));
    if angDiff(i,1)>180
        angDiff(i,1)=360-angDiff(i,1);
    end    
end

% interpolat to range between -1 and 1
outputRes=0.01;
alpha=-1:outputRes:1;
beta=-1:outputRes:1;
[A,B]=meshgrid(alpha,beta);
CR=griddata(PU(:,1),PV(:,2),angDiff(:,1),A,B);

CR=abs(flipud(CR))
figure
imagesc(CR)
axis equal
axis off


figure
contourf(alpha,beta,flipud(CR))
axis equal
axis off

c=colorbar('northoutside');
c.Label.String='angle difference (\circ)';
caxis([0 180]);
c.XTick=[0 45 90 135 180];
c.FontSize=14;
    
end


