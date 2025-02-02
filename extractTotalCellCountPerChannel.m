
%extractTotalCellCountPerChannel

% need to open the appropriate pooledmap file (with Itrans)

N=[];
D=[];
T=[];
V=[];
t=reshape(pooledmap.cellCount,[],1);
for i=1:size(t,1)
N=[N;pooledmap.cellCount{i}(1)];
D=[D;pooledmap.cellCount{i}(2)];
T=[T;pooledmap.cellCount{i}(3)];
V=[V;pooledmap.cellCount{i}(4)];
end

sumN=sum(N);
sumD=sum(D);
sumT=sum(T);
sumV=sum(V);

save('totalCellCountPerChannel_translation_ON_DSI0.3.mat');


%% pooling all data from all DSI conditions
onoff05=[sumN,sumD,sumT,sumV];
save('summary_all DSI conditions across all channels.mat');

groups={'nasal','dorsal','temporal','ventral'};
per=[12.09	43.20	55.26	52.50];     %This is the percentage of cell remanining DS when moving from DSI threshold of 0.3 to 0.5 
per6bin=[5.13	40.00	56.92	47.54];     %This is the percentage of cell remanining DS when moving from DSI threshold of 0.3 to 0.5 based on 6 central bins only 

figure;
hold on
bar(1,per(1),'b');
bar(2,per(2),'g');
bar(3,per(3),'r');
bar(4,per(4),'m');
set(gca,'FontSize',16)
xlabel('Channel','FontSize',20);ylabel('Cells identified as DS (%)','FontSize',20);
set(ax,'box','off');
box off
set(gcf, 'color', [1 1 1]);
ax = gca;
ax.XTick = 1:1:4;
ax.XTickLabel = {'nasal','dorsal','temporal','ventral'};



figure;
hold on
bar(1,per6bin(1),'b');
bar(2,per6bin(2),'g');
bar(3,per6bin(3),'r');
bar(4,per6bin(4),'m');
set(gca,'FontSize',16)
xlabel('Channel','FontSize',20);ylabel('Cells identified as DS (%)','FontSize',20);
set(ax,'box','off');
box off
set(gcf, 'color', [1 1 1]);
ax = gca;
ax.XTick = 1:1:4;
ax.XTickLabel = {'nasal','dorsal','temporal','ventral'};






