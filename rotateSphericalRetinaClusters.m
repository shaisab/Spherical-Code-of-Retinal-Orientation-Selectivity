% rotate spherical retina with clusters

%% set viewing angle to stright down
az = 0;
el = 90;
view(az, el);
set(gca, 'visible', 'off');
set(gca,'color',[1 1 1]);

%% set viewing angle to side and down
az = 0;
el = 45;
view(az, el);
set(gca, 'visible', 'off');
set(gca,'color',[1 1 1]);

%% set viewing angle to side and down
az = 180;
%el = 50;
el = 45;
view(az, el);
set(gca, 'visible', 'off');
set(gca,'color',[1 1 1]);

%% rotate camera orbit
camorbit(3,20)
