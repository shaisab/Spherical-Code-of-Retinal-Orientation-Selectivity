
for i=1:size(metaRep.meanRep.ROIdata,2)
DSI(i,1)=metaRep.meanRep.ROIdata{1, i}.summary.DSI;
end

for i=1:size(metaRep.meanRep.ROIdata,2)
isDS(i,1)=metaRep.meanRep.ROIdata{1, i}.summary.isDS;
end

mean(DSI(isDS==1))


for i=1:size(metaRep.meanRep.ROIdata,2)
OSI(i,1)=metaRep.meanRep.ROIdata{1, i}.summary.OSI;
end

for i=1:size(metaRep.meanRep.ROIdata,2)
isOS(i,1)=metaRep.meanRep.ROIdata{1, i}.summary.isOS;
end

mean(OSI(isOS==1))


%%%%%%%%%%
for i=1:size(metaRep.meanRep.ROIdata,2)
ONorONOFForOFF(i,1)=metaRep.meanRep.ROIdata{1, i}.summary.ONorONOFForOFF;
end

nON=sum(ONorONOFForOFF==1);
nONOFF=sum(ONorONOFForOFF==2);
nOFF=sum(ONorONOFForOFF==3);
