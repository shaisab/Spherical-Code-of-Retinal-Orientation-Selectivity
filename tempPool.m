function tempPool

[filedir,go]=bbdirselector('select data folder',cd);
if ~go disp('no folder selected'); return; end

namesASC=getFileList(filedir,['vecComp_','50_DS_ASC'],0,'anywhere');
for k=1:size(namesASC,2)
    ind=strfind(namesASC{k},'_');
    expDate{k}=namesASC{k}(ind(end-2)+1:ind(end)+4);
    
namesMap=getFileList(filedir,['vecComp_','50','_','DS'],0,'anywhere');
    for h=1:size(namesMap,2)  
    if sum(strfind(namesMap{h},['ASC_',expDate{k}])) && sum(strfind(namesMap{h},'ASC'))
        vecCompASC{k}=load(num2str(namesMap{h}));
    elseif sum(strfind(namesMap{h},['PSC_',expDate{k}])) && sum(strfind(namesMap{h},'PSC'))
        vecCompPSC{k}=load(num2str(namesMap{h}));
    elseif sum(strfind(namesMap{h},['LSC_',expDate{k}])) && sum(strfind(namesMap{h},'LSC'))
        vecCompLSC{k}=load(num2str(namesMap{h}));
    end
  end
end         %TMP end
end