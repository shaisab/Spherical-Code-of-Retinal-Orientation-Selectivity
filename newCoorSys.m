

P=rad2deg(atan2(STHdata(:,4),STHdata(:,3)));
P1=[P(P>0); 360+P(P<0)];
pdir=rad2deg(STHdata(:,2))-P1;
pdir1=[pdir(pdir>0); 360+pdir(pdir<0)];
[X,Y]=pol2cart(deg2rad(pdir1), 1);
compassSS2p(X,Y,1,'r',3);