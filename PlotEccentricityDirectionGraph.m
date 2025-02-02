hold on;
set(gcf, 'color', [1 1 1]);
%plot([0 360],[90 90],'--m','LineWidth',3);
grid on
axis equal
ax = gca;
ax.XTick = 0:45:360;
ax.YTick = 0:45:180;
ax.XTickLabel = {'P','315','I','225','A','135','S','45','P'};
ax.FontSize=14;
set(ax,'Xdir','reverse');
%ax.YTickLabel = {'','','','',''};
%ax.YTickLabel = {'anti-OD','135', 'margin','45','OD'};
ax.YTickLabel = {'OA','45','margin','135','anti-OA'};

%xlabel('Direction (deg.)','FontSize',20);     % for single maps
xlabel('Direction (deg.)','FontSize',16);
%ax.YTickLabel = {'','','','',''};

%ylabh=ylabel('Eccentricity (deg.)','FontSize',20);     % for single maps
ylabh=ylabel('Eccentricity (deg.)','FontSize',16);
set(ylabh, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);

c=colorbar('northoutside');
%c=colorbar('eastoutside');
%c.AxisLocation='out';
c.Label.String='Goodness of fit (%)';
%c.Label.String='Spiking output';
%c.FontSize=18;     % for single maps
c.FontSize=14;
c.Limits=[0 40];


