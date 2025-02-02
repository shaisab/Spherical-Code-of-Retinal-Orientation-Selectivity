function [] = plotChannelDistribution(distdata,datatype,singleChannels,dataname)


if datatype==0  % positive percentages only
lowerc=0;
upperc=100;
elseif datatype==1  % positive and negative percentages
lowerc=-40;
upperc=40;
end    
        
xbin=1:10;
ybin=1:10;

switch singleChannels
    case 1
f1=figure;
set(f1,'position',[450 150 950 800]);
subplot(2,2,1)
plotStandardRetina
hold on
subplot(2,2,1)
contourf(xbin,ybin,distdata.c1,5)
colormap(jet)
caxis([lowerc,upperc]);
colorbar
title('channel 1');
hold off

subplot(2,2,2)
plotStandardRetina
hold on
subplot(2,2,2)
contourf(xbin,ybin,distdata.c2,5)
colormap(jet)
caxis([lowerc,upperc]);

colorbar
title('channel 2');
hold off

subplot(2,2,3)
plotStandardRetina
hold on
subplot(2,2,3)
contourf(xbin,ybin,distdata.c3,5)
colormap(jet)
caxis([lowerc,upperc]);

colorbar
title('channel 3');
hold off

subplot(2,2,4)
plotStandardRetina
hold on
subplot(2,2,4)
contourf(xbin,ybin,distdata.c4,5)
colormap(jet)
caxis([lowerc,upperc]);
colorbar
title('channel 4');
hold off

subtitle(dataname);

    case 0
f1=figure;
set(f1,'position',[450 150 950 800]);        
plotStandardRetina
hold on
contourf(xbin,ybin,distdata,5)
colormap(jet)
caxis([lowerc,upperc]);
colorbar
title('diff count between ON and ONOFF DS cells');
hold off     
end

end