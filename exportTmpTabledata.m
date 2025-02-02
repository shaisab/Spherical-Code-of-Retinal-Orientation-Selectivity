


[mergedmetaFileList]=getFileList(cd,'mergedmeta_fov',0,'anywhere');
[TMP_tabledataFileList]=getFileList(cd,'TMP_tabledata_fov',0,'anywhere');
for k=1:size(mergedmetaFileList,2)
load(mergedmetaFileList{k})

for j=1:size(TMP_tabledataFileList,2) 
if strcmp(mergedmetaFileList{k}(end-9:end-4),TMP_tabledataFileList{j}(end-9:end-4))
load(TMP_tabledataFileList{j})

if sum(mergedmeta.isCART)+sum(mergedmeta.isRBPMS)+sum(mergedmeta.isRetro)==0
try mergedmeta.isCART=a(:,4); end
try mergedmeta.isRBPMS=a(:,5); end
try mergedmeta.isRetro=a(:,6); end

save(([mergedmetaFileList{k}]), 'mergedmeta');
end
end
end
end
