
function []= mergeMetaRepFiles()

olddir=cd;
thisdir = dir;
str = {thisdir.name};
[s,v] = listdlg('PromptString','Select repetition folders :',...
    'SelectionMode','multiple',...
    'ListString',str);
if ~v disp('no folder selected'); return; end

for i=1:size(s,2)
    filedir=thisdir(s(i)).name;
    [metaFileList{i}]=getFileList(filedir,'metaRep',0,'anywhere');
end

if size(metaFileList,2)~=size(s,2)
    h=warndlg('The number of meta files is not correct');
    uiwait(h);
    return
end
mergedmeta.metaFileList=metaFileList;

for k=1:size(metaFileList,2)
    load([olddir,'\',num2str(cell2mat(metaFileList{k}))]);
    if strcmp(metaRep.metadata{1}.log.expType,'DG')
        mergedmeta.meanRepDG=metaRep.meanRep;
        mergedmeta.metadataDG=metaRep.metadata;
        for ROIn=1:size(metaRep.metadata{1,1}.ROIdata,2)
            mergedmeta.isDS(ROIn,1)=metaRep.meanRep.ROIdata{ROIn}.summary.isDS;
            mergedmeta.isOS(ROIn,1)=metaRep.meanRep.ROIdata{ROIn}.summary.isOS;
            mergedmeta.isVertical(ROIn,1)=metaRep.meanRep.ROIdata{ROIn}.summary.isVertical;
            mergedmeta.isHorizontal(ROIn,1)=metaRep.meanRep.ROIdata{ROIn}.summary.isHorizontal;
            mergedmeta.oriFWHM(ROIn,1)=metaRep.meanRep.ROIdata{ROIn}.summary.oriFWHM;
            mergedmeta.pori2(ROIn,1)=metaRep.meanRep.ROIdata{ROIn}.summary.pori2;
            mergedmeta.symmetry_ratio(ROIn,1)=metaRep.meanRep.ROIdata{ROIn}.summary.symmetry_ratio;
            mergedmeta.sumCorrRep2(ROIn,1)=metaRep.meanRep.ROIdata{ROIn}.summary.sumCorrRep2;
        end
    elseif strcmp(metaRep.metadata{1}.log.expType,'DB')
        mergedmeta.meanRepDB=metaRep.meanRep;
        mergedmeta.metadataDB=metaRep.metadata;
        for ROIn=1:size(metaRep.metadata{1,1}.ROIdata,2)
            mergedmeta.ONorONOFForOFF(ROIn,1)=metaRep.meanRep.ROIdata{ROIn}.summary.ONorONOFForOFF;
        end
    end
end

save((['mergedmeta_',filedir(1:6),'.mat']), 'mergedmeta');

end
