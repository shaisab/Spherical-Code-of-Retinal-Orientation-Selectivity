

figure
set(gcf, 'color', [1 1 1]);
%bar([xn180(1),xn180(3);xn180(2),xn180(4)],0.8,'grouped');
bar([xn(1),xn(3);xn(2),xn(4)],0.8,'grouped');
%bar([xn(2),xn(4);xn(3),xn(1)],0.8,'grouped');       % use for ON rotation fit
%ax.XTickLabel={'nasal','temporal';'dorsal','ventral'};
xlim([0.65 2.4]);

set(gca,'FontSize',16)
xlabel('Channel','FontSize',20);ylabel('Relative contribution','FontSize',20);
box off
title(['\itR\rm^2 = ',num2str(R2,2)]);
%title(['\itR\rm^2 = ',num2str(R2_180,2)]);