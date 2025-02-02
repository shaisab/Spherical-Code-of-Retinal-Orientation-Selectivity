
function [] = compareCellCounts()

f1=figure;
set(f1,'position',[450 150 950 800]);

[FileName,PathName] = uigetfile('*.mat','Select a cellCountONandONOFF file');
load(FileName);
pooledmap_ONOFF.cellCount=flipud(pooledmap_ONOFF.cellCount);
pooledmap_ON.cellCount=flipud(pooledmap_ON.cellCount);

for i=1:size(pooledmap_ON.cellCount,1)
    for j=1:size(pooledmap_ON.cellCount,2)
        percentONOFF{i,j}=100*(pooledmap_ONOFF.cellCount{i,j}./sum(pooledmap_ONOFF.cellCount{i,j}));
        percentON{i,j}=100*(pooledmap_ON.cellCount{i,j}./sum(pooledmap_ON.cellCount{i,j}));
        diff{i,j}=percentONOFF{i,j}-percentON{i,j}; 
        
        
        
        
        
        
    end
end



figure
y=diff{3,5};                  %the data.
fHand = figure;
aHand = axes('parent', fHand);
hold(aHand, 'on')
%colors = hsv(numel(y));
colors = [0 0 1; 0 1 0; 1 0 0; 0 0 0];
for i = 1:numel(y)
    bar(i, y(i), 'parent', aHand, 'facecolor', colors(i,:));
end
set(gca, 'XTick', 1:numel(y), 'XTickLabel', {'c1', 'c2', 'c3', 'c4'})
ylim([-40 40]);



end
