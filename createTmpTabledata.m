

[mergedmetaFileList]=getFileList(cd,'mergedmeta_fov',0,'anywhere');
for k=1:size(mergedmetaFileList,2)
load(mergedmetaFileList{k})

a1=mergedmeta.isDS;
a2=mergedmeta.isOS;
a3=mergedmeta.ONorONOFForOFF;
try a4=mergedmeta.isCART; end
try a5=mergedmeta.isRBPMS; end
try a6=mergedmeta.isRetro; end
try a=[a1,a2,a3,a4,a5,a6]; end

tmp=num2str(mergedmetaFileList{k});
filename=tmp(end-9:end-4);
dirname=tmp(1:end-21);
save(([dirname,'TMP_','tabledata_',filename,'.mat']), 'a');
end