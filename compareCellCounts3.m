
function [] = compareCellCounts3()



% plot contout plots of diff percentage between ONOFF and ON cells, for each of four translation channels
[FileName,PathName] = uigetfile('*.mat','Select a cellCountONandONOFF file');
load(FileName);
% pooledmap_ONOFF.cellCount=flipud(pooledmap_ONOFF.cellCount);
% pooledmap_ON.cellCount=flipud(pooledmap_ON.cellCount);

for i=1:size(pooledmap_ON.cellCount,1)
    for j=1:size(pooledmap_ON.cellCount,2)
        percentONOFF{i,j}=100*(pooledmap_ONOFF.cellCount{i,j}./sum(pooledmap_ONOFF.cellCount{i,j}));
        percentON{i,j}=100*(pooledmap_ON.cellCount{i,j}./sum(pooledmap_ON.cellCount{i,j}));
        diff{i,j}=percentONOFF{i,j}-percentON{i,j};         
    end
end


for c=1:size(diff{1,1},2)
    fieldname=['c',num2str(c)];
    for i=1:size(diff,1)
        for j=1:size(diff,2)
            singleDiff.(fieldname)(i,j)=diff{i,j}(c);
        end
    end
end

xbin=1:10;
ybin=1:10;

[StU1,StV1,StU2,StV2,StU3,StV3,StU4,StV4]=plotStandardRetina3;
layer=ones(10,10);

for i=1:10,
    for j=1:10,
%         min([
%     end
% end
size(xbin)
size(singleDiff.c1)
f1=figure;
set(f1,'position',[450 150 950 800]);
subplot(2,2,1)
plotStandardRetina;
hold on
subplot(2,2,1)
contourf(xbin,ybin,singleDiff.c1)
colormap(jet)
caxis([-40,40]);
colorbar
title('channel 1');
hold off

subplot(2,2,2)
plotStandardRetina;
hold on
subplot(2,2,2)
contourf(xbin,ybin,singleDiff.c2)
colormap(jet)
caxis([-40,40]);
colorbar
title('channel 2');
hold off

subplot(2,2,3)
plotStandardRetina;
hold on
subplot(2,2,3)
contourf(xbin,ybin,singleDiff.c3)
colormap(jet)
caxis([-40,40]);
colorbar
title('channel 3');
hold off

subplot(2,2,4)
plotStandardRetina;
hold on
subplot(2,2,4)
contourf(xbin,ybin,singleDiff.c4)
colormap(jet)
caxis([-40,40]);
colorbar
title('channel 4');
hold off

    end
end
