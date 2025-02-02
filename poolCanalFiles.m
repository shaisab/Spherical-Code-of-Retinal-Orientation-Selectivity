
function []= poolCanalFiles(varargin);
clc;

%cd('helperfiles'); a=cd; addpath(a); cd ...

[thedir,go]=bbdirselector(cd,'select data folder');
if ~go return; end

convert_oib_to_tif(thedir);

names=getFileList(thedir,'.lsm',0,'endonly');
if length(names)>0
    
    olddir=cd;
    cd(thedir);
    for i = 1:length(names)
        thename=shortfile(names{i})
        thename=thename(1:end-4);
        if exist(thename)~=7
            mkdir(thename)
        end
        [STATUS,MESSAGE,MESSAGEID] = movefile(names{i},thename,'f');
    end
    cd(olddir);
else  
    names=getFileList(thedir,'.tif',0,'endonly');
    olddir=cd;
    cd(thedir);
    for i = 1:length(names)
        thename=shortfile(names{i});
        thename=thename(1:end-4);
        if exist(names{i})~=7
            mkdir(thename)
        end
        [STATUS,MESSAGE,MESSAGEID] = movefile(names{i},thename,'f');
    end
    
    names=getFileList(thedir,'.oib',0,'endonly');
    for i = 1:length(names)
        thename=shortfile(names{i});
        thename=thename(1:end-4);
        [STATUS,MESSAGE,MESSAGEID] = movefile(names{i},thename,'f');
    end
    
    matnames=getFileList(thedir,'.mat',0,'endonly');
     for i = 1:length(matnames)
        log=load(matnames{i});
        [~, ~, ~, logH(i),logMN(i),logS(i)]=datevec(log.startTime);
        thename=shortfile(matnames{i})
        thename=thename(1:22);
        [STATUS,MESSAGE,MESSAGEID] = movefile(matnames{i},thename,'f');
     end
     
     names=getFileList(thedir,'.txt',0,'endonly');
     for i = 1:length(names)
        thename=shortfile(names{i});
        thename=thename(1:end-4);
        if exist(thename)==7
        [STATUS,MESSAGE,MESSAGEID] = movefile(names{i},thename,'f');
        end
     end
     
     names=getFileList(thedir,'.abf',0,'endonly');
     for i = 1:length(names)
        thename=shortfile(names{i});
        thename=thename(1:end-4);
        [abfData, si, h] = abfload(names{i},'channels','a');
        [time_string,abfH(i),abfMN(i),abfS(i)]=secs2hms(h.lFileStartTime);
        for j=1:size(logH,1)
            if abfH(i)==logH(j)
                if abs((abfS(i)+(abfMN(i)*60))-(logS(j)+(logMN(j)*60)))<15  % allows maximum of 15 seconds difference between the start time of the log file and the abf file
                  [STATUS,MESSAGE,MESSAGEID] = movefile(names{i},shortfile(matnames{j}(1:end-12)),'f');  
                end
            end
        end
     end
  
%ss_organizeScimData(thedir);     

cd(olddir);
javaaddpath(which('MatlabGarbageCollector.jar'))
clear all; clear java; jheapcl;

end