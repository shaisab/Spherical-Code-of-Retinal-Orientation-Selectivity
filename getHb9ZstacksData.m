


filedir=cd;
[Hb9xyposFileList]=getFileList(filedir,'Hb9position.txt',0,'anywhere');
Hb9xypos=tdfread(Hb9xyposFileList{1});  

[xyposFileList]=getFileList(filedir,'xyposition.txt',0,'anywhere');
xypos=tdfread(xyposFileList{1});
    
metadata.ScanXYZ.relXPosition=xypos.relXPosition(xyind);
metadata.ScanXYZ.relYPosition=xypos.relYPosition(xyind);