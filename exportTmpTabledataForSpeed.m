


[TMP_tabledataFileList]=getFileList(cd,'TMP_tabledataForSpeed_fov',0,'anywhere');
s=[];
for j=1:size(TMP_tabledataFileList,2)
    fieldname=num2str([TMP_tabledataFileList{j}(end-31:end-25), TMP_tabledataFileList{j}(end-9:end-7),TMP_tabledataFileList{j}(end-5:end-4)]);
    %s=[s;load(TMP_tabledataFileList{j})];
    %s.(fieldname)=load(TMP_tabledataFileList{j});
    s=[s;cell2mat(struct2cell(load(TMP_tabledataFileList{j})))];
end
