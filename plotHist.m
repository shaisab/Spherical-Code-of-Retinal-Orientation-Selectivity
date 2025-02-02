
%% plot hitogram of canals from polledmap files
% need to load pooledmap file and change the canal identity below.

angleXYZ=rad2deg(pooledmap.XYZvecComp.angleXYZASC)';

%angleXYZ90=[angleXYZ(find(angleXYZ<=90)),180-angleXYZ(find(angleXYZ>90))];

figure;
histogram(angleXYZ,15);
ax=gca; xlim([0 180]); ax.XTick=[0:45:180];
set(gca,'FontSize',24)
%xlabel('angle difference empirical - prediced (degrees)','FontSize',20);ylabel('number of cells','FontSize',20);
xlabel(' \theta (degrees)','FontSize',30);ylabel('number of cells','FontSize',30);
set(ax,'box','off','color',[1 1 1]);

drawnow expose





%% plot hitogram of canals
% need to load vecComp file for the requested canal
angleXYZ=rad2deg(XYZvecComp.angleXYZ);
angleXYZ90=[angleXYZ(find(angleXYZ<=90)),180-angleXYZ(find(angleXYZ>90))];

figure;
subplot(2,1,1); histogram(angleXYZ,15); 
ax=gca; xlim([0 180]); ax.XTick=[0:45:180];
subplot(2,1,2); histogram(angleXYZ90,15);
ax=gca; xlim([0 90]); ax.XTick=[0:45:90];
drawnow expose

%% plot hitogram of canals
% need to load vecComp file for the requested canal
angleXYZ=rad2deg(XYZvecComp.angleXYZ);
angleXYZ90=[angleXYZ(find(angleXYZ<=90)),180-angleXYZ(find(angleXYZ>90))];

figure;
histogram(angleXYZ90,15);
ax=gca; xlim([0 90]); ax.XTick=[0:45:90];
drawnow expose
%% plot spherical contour with one color
% NEED to first run ContourPlotAlphaBeta until the spherical retina
% section


figure;
hold on;
surf(Xshift,Yshift,Zshift);
shading interp;
axis equal
alpha(.8)
set(findobj(gcf, 'type','axes'), 'Visible','off')

thetafine=linspace(0,2*pi,40);
sfine=linspace(pi/2,pi,10);

[THfine,Sfine]=meshgrid(thetafine,sfine);

Xfine=(R+0.01)*cos(THfine).*sin(Sfine);
Yfine=(R+0.01)*sin(THfine).*sin(Sfine);
Zfine=-(R+0.01)*cos(Sfine);

surf(Xfine,Yfine,Zfine,'FaceColor','none');



%%