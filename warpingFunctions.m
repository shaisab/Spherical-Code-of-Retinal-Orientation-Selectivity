

load mandrill
colormap(map)
[x,y,z] = cylinder;
Xhalf = [ones(480,375)*max(max(X))/2, X, ...
ones(480,125)*max(max(X))/2];
surface(x,y,z, 'FaceColor','texturemap',...
'EdgeColor','none','Cdata',flipud(Xhalf))
view(3)

[x,y,z] = sphere;
h = surf(x,y,z);
set(h, 'FaceColor','texturemap','EdgeColor','none','Cdata',flipud(Xhalf))

f = imread('ONOFFDS superimposed on equal contribution four translation channels based on ONOFF locations.tif');
warp(f);
rotate3d