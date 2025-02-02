

alphaON=[5,100,180,285];
%alphaON=360-alphaON;    %equivalent to fliplr
betaON=[140,105,55,70];
alphaONOFF=[10,95,185,275];
%alphaONOFF=360-alphaONOFF;    %equivalent to fliplr
betaONOFF=[120,105,60,70];
coeffON=[0.167274554714880,0.350309766818000,0.233410692512288,0.249004985954832];
coeffONOFF=[0.414804628989209,0.159656043199595,0.175192315576808,0.250347012234388];



figure;
s1=scatter(alphaON,betaON,600.*coeffON,'b','filled');
hold on;
s2=scatter(alphaONOFF,betaONOFF,600.*coeffONOFF,'r','filled');
xlim([0 360]);ylim([0 180]);
hold on;
set(gcf, 'color', [1 1 1]);
plot([45 405],[90 90],'--m','LineWidth',2);
%grid on
ax = gca;
ax.XTick = 0:45:360;
ax.YTick = 0:45:180;
ax.XTickLabel = {'P','315','I','225','A','135','S','45'};
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
legend([s1,s2],{'ON','ON-OFF',},'Location','northeast','FontSize',14);
legend('boxoff');





