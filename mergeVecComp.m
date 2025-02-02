function mergeVecComp(discPoints)

[filedir,go]=bbdirselector('select data folder',cd);
if ~go disp('no folder selected'); return; end

namesMap=getFileList(filedir,['vecComp_',num2str(discPoints),'_','DS'],0,'anywhere');
for k=1:size(namesMap,2)
    if strfind(namesMap{k},'ASC')
        vecCompASC=load(num2str(namesMap{k}));
    elseif strfind(namesMap{k},'PSC')
        vecCompPSC=load(num2str(namesMap{k}));
    elseif strfind(namesMap{k},'LSC')
        vecCompLSC=load(num2str(namesMap{k}));
    end
end
    
    
      
    
end