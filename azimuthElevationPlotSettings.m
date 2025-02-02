function []= azimuthElevationPlotSettings()

hold on
ax=gca;
set(gcf, 'color', [1 1 1]);
plot([-180 180],[0 0],'m','LineWidth',1);
plot([0 0],[-90 90],'m','LineWidth',1);
xlim([-180 180]);ylim([-90 90]);
axis equal tight
grid on
ax.XTick = -180:45:180;
ax.YTick = -90:45:90;
ax.XTickLabel = {'behind','-135','left','-45','ahead','45','right','135','behind'};
ax.FontSize=12;
%xlabel('Azimuth (deg.)','FontSize',16);
%ax.YTickLabel = {'','','','',''};

%ylabh=ylabel('Elevation (deg.)','FontSize',16);
%set(ylabh, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
ax.YTickLabel = {'zenith','45', 'horizon','-45','nadir'};       % for confusion maps
%ax.YTickLabel = {'nadir','-45', 'horizon','45','zenith'};      % for
%global maps

% %c=colorbar('northoutside');
% c=colorbar('eastoutside');
% c.AxisLocation='out';
% c.Label.String='Spiking output';
% %c.FontSize=18;     % for single maps
% c.FontSize=14;
% c.Limits=[0 0.9];


end