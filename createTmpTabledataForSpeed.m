

[mergedmetaFileList]=getFileList(cd,'mergedmeta_fov',0,'anywhere');
for k=1:size(mergedmetaFileList,2)
load(mergedmetaFileList{k})

a1=mergedmeta.isDS;
a2=mergedmeta.isOS;
a3=mergedmeta.ONorONOFForOFF;
for i=1:size(mergedmeta.meanRepDG.ROIdata,2)
a4(i,1)=mergedmeta.meanRepDG.ROIdata{1,i}.summary.DSI; 
a5(i,1)=mergedmeta.meanRepDG.ROIdata{1,i}.summary.OSI;  
a6(i,1)=mergedmeta.meanRepDG.ROIdata{1,i}.summary.globalOSI;  
end
try a=[a1,a2,a3,a4,a5,a6]; end

tmp=num2str(mergedmetaFileList{k});
filename=tmp(end-9:end-4);
dirname=tmp(1:end-21);
save(([dirname,'TMP_','tabledataForSpeed_',filename,'.mat']), 'a');
clear a a1 a2 a3 a4 a5 a6
end