function plotPolarAndQuiver(map)

for i=1:length(map.poriSz)
    if isnan(map.poriSz(i))
    elseif map.oriR2(i)<0.6
        rho=[map.meanr(i,:),map.meanr(i,1)];

 
% create a polar plot
f1=figure;
set(f1, 'color', [1 1 1]);
set(f1,'position',[50 50 500 500]);

P=polarSS2p(deg2rad(map.theta(i,:)),rho./max(rho),1,[0.3010 0.7450 0.9330],6);
set(findall(gcf, 'String', '45', '-or','String','135', '-or','String','225', '-or','String','315') ,'String', '  ');
set(findall(gcf, 'String', '0'),'String', ' N');
set(findall(gcf, 'String', '90'),'String', ' D');
set(findall(gcf, 'String', '180'),'String', ' T');
set(findall(gcf, 'String', '270'),'String', ' V');
ax = ancestor(P(1),'axes');
th = findall(ax,'Type','Text');
set(th, 'FontSize', 26);

hold on;

% create a quiver plot
[X,Y]=pol2cart(deg2rad([map.poriSz(i),map.poriSz(i)+180]'), 1);
h=quiver(zeros(size(X,1),1),zeros(size(X,1),1),[X],[Y],1,'r','LineWidth',8);
%set(findobj(gcf, 'type','axes'), 'Visible','off')
%set(gcf, 'color', [1 1 1]);
axis equal;
h.ShowArrowHead='off';

hold off
    end
end
end






