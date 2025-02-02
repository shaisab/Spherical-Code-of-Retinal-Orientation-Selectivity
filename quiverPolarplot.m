
%% plot a polar plot for all Hb9 physiology cells (from allCells file)

[THETA,RHO]=cart2pol(allCells.UVStvecComp.VecUSt,allCells.UVStvecComp.VecVSt)
[X,Y]=pol2cart(THETA,RHO-0.02);
figure;
compassSS2p(X,Y,1,'r',1);

%% plot a polar plot for all CART negative cells (from pooledmap file)

[THETA,RHO]=cart2pol(pooledmap.UVStvecComp.VecUSt,pooledmap.UVStvecComp.VecVSt)
[X,Y]=pol2cart(THETA,RHO-0.02);
figure;
compassSS2p(X,Y,1,'r',1);

%% plot a polar plot for all Hb9 morphology cells (from allCells file) with the length of vectors proportional to DSI
figure;
[X,Y]=pol2cart([deg2rad(map.pdirSz);deg2rad(map.pdirSz(1,1))],[map.DSI;1]);
%[X,Y]=pol2cart(deg2rad(map.pdirSz),map.DSI);
figure;
compassSS2p(X,Y,1,'b',1);

%% plot a polar plot for all Hb9 morphology cells (from map file) with the length of vectors proportional to DSI

[THETA,RHO]=cart2pol(allCells.UVStvecComp.VecUSt,allCells.UVStvecComp.VecVSt)
[X,Y]=pol2cart(THETA,1.55*allCells.DSI);
figure;
compassSS2p(X,Y,1,'b',1);

%% plot a quiver polar plot for all cells in pooledmap
figure; 
h=quiver(repmat(0,size(pooledmap.UVStvecComp.VecUSt,1),1),repmat(0,size(pooledmap.UVStvecComp.VecUSt,1),1),pooledmap.UVStvecComp.VecUSt,pooledmap.UVStvecComp.VecVSt,1,'k');
set(findobj(gcf, 'type','axes'), 'Visible','off')
set(gcf, 'color', [1 1 1]);
axis equal;
h.ShowArrowHead='off';

figure;
[angle,~]=cart2pol(pooledmap.UVStvecComp.VecUSt,pooledmap.UVStvecComp.VecVSt);
rose(angle,40);
%% plot quiver polar plot for cluster 1,VecUSt and VecVSt % change clusern
clustern=56;
figure; 
h=quiver(repmat(0,size(pooledmap.clusterdata{1, clustern}.UVStvecComp.VecUSt,1),1),repmat(0,size(pooledmap.clusterdata{1, clustern}.UVStvecComp.VecUSt,1),1),pooledmap.clusterdata{1, clustern}.UVStvecComp.VecUSt,pooledmap.clusterdata{1, clustern}.UVStvecComp.VecVSt ,1,'k','LineWidth',1.5);
set(findobj(gcf, 'type','axes'), 'Visible','off')
set(gcf, 'color', [1 1 1]);
%axis[-1 1 -1 1]
axis equal;
h.ShowArrowHead='off';

figure;
[angle,~]=cart2pol(pooledmap.clusterdata{1, clustern}.UVStvecComp.VecUSt,pooledmap.clusterdata{1, clustern}.UVStvecComp.VecVSt);
rose(angle,40);

%% plot quiver polar plot for of 4 bins 1,VecUSt and VecVSt % change clusern
clustern=[45,46,55,56];     % central cluster
%clustern=[33,34,43,44];     % ventral-temporal cluster
%clustern=[24,25,34,35];

vecX=[pooledmap.clusterdata{1, clustern(1,1)}.UVStvecComp.VecUSt;pooledmap.clusterdata{1, clustern(1,2)}.UVStvecComp.VecUSt;...
    pooledmap.clusterdata{1, clustern(1,3)}.UVStvecComp.VecUSt;pooledmap.clusterdata{1, clustern(1,4)}.UVStvecComp.VecUSt];
vecY=[pooledmap.clusterdata{1, clustern(1,1)}.UVStvecComp.VecVSt;pooledmap.clusterdata{1, clustern(1,2)}.UVStvecComp.VecVSt;...
    pooledmap.clusterdata{1, clustern(1,3)}.UVStvecComp.VecVSt;pooledmap.clusterdata{1, clustern(1,4)}.UVStvecComp.VecVSt];

figure; 
h=quiver(zeros(size(vecX,1),1),zeros(size(vecX,1),1),vecX,vecY,1,'k','LineWidth',1.5);
set(findobj(gcf, 'type','axes'), 'Visible','off')
set(gcf, 'color', [1 1 1]);
%axis[-1 1 -1 1]
axis equal;
h.ShowArrowHead='off';

figure;
[angle,~]=cart2pol(vecX,vecY);
rose(angle,40);

%% plot quiver polar plot for cluster 1,ASC predictions % change clusern
clustern=5;
figure; quiver(repmat(0,size(pooledmap.clusterdata{1, clustern}.VecUStcompASC,1),1),repmat(0,size(pooledmap.clusterdata{1, clustern}.VecUStcompASC,1),1),pooledmap.clusterdata{1, clustern}.VecUStcompASC,pooledmap.clusterdata{1, clustern}.VecVStcompASC ,1,'b');
hold on
quiver(repmat(0,size(pooledmap.clusterdata{1, clustern}.VecUStcompASC,1),1),repmat(0,size(pooledmap.clusterdata{1, clustern}.VecUStcompASC,1),1),-pooledmap.clusterdata{1, clustern}.VecUStcompASC,-pooledmap.clusterdata{1, clustern}.VecVStcompASC ,1,'b');
quiver(repmat(0,size(pooledmap.clusterdata{1, clustern}.VecUStcompASC,1),1),repmat(0,size(pooledmap.clusterdata{1, clustern}.VecUStcompASC,1),1),pooledmap.clusterdata{1, clustern}.VecUStcompPSC,pooledmap.clusterdata{1, clustern}.VecVStcompPSC ,1,'m');
quiver(repmat(0,size(pooledmap.clusterdata{1, clustern}.VecUStcompASC,1),1),repmat(0,size(pooledmap.clusterdata{1, clustern}.VecUStcompASC,1),1),-pooledmap.clusterdata{1, clustern}.VecUStcompPSC,-pooledmap.clusterdata{1, clustern}.VecVStcompPSC ,1,'m');
quiver(repmat(0,size(pooledmap.clusterdata{1, clustern}.VecUStcompASC,1),1),repmat(0,size(pooledmap.clusterdata{1, clustern}.VecUStcompASC,1),1),pooledmap.clusterdata{1, clustern}.VecUStcompLSC,pooledmap.clusterdata{1, clustern}.VecVStcompLSC ,1,'r');
quiver(repmat(0,size(pooledmap.clusterdata{1, clustern}.VecUStcompASC,1),1),repmat(0,size(pooledmap.clusterdata{1, clustern}.VecUStcompASC,1),1),-pooledmap.clusterdata{1, clustern}.VecUStcompLSC,-pooledmap.clusterdata{1, clustern}.VecVStcompLSC ,1,'r');
set(findobj(gcf, 'type','axes'), 'Visible','off')
set(gcf, 'color', [1 1 1]);
axis equal;
hold off


