
function []= extractBinValues()

filedir=cd;
physDBFileList=getFileList(filedir,'physdata_DB',0,'anywhere');
physDGFileList=getFileList(filedir,'physdata_DG',0,'anywhere');
        
i=1;
if ~isempty(physDBFileList)
    for ncell=1:size(physDBFileList,2)
        curname=physDBFileList{ncell};
        inds=find(curname=='_');
             DS.physDB(i)=load(num2str(physDBFileList{ncell}));   
             DS.physDG(i)=load(num2str(physDGFileList{ncell}));   
             DS.spiketimes{1,i}=DS.physDB(i).physdata.all.spiketimesAll;
             DS.spikeArray{1,i}=DS.physDB(i).physdata.all.spikeArrayAll;
             DS.pdir(1,i)=DS.physDG(i).physdata.all.pdir;
             i=i+1;
    end
end


%% extract the bin values for each DS cell and each direction (repetitions pooled)
ftmp=figure;
tmp={};
for c=1:size(DS.spiketimes,2)                 % cell#
    for n=1:size(DS.spiketimes{1,1},2)        % direction#
        for m=1:size(DS.spiketimes{1,1},1)    % repitition#          
            tmp=[tmp,DS.spiketimes{1,c}{m,n}];       
        end
        spiketimesT=(cell2mat(tmp).*DS.physDB(1).physdata.h(1).si)./1000000;          % spike times in seconds
        htmp=histogram(spiketimesT);                  % plot PSTH
        htmp.BinEdges=0:1:floor(DS.physDB(1).physdata.rep{1,1}.phystime(end,1)-DS.physDB(1).physdata.rep{1,1}.phystime(1,1));
        htmp.NumBins=floor(DS.physDB(1).physdata.rep{1,1}.phystime(end,1)-DS.physDB(1).physdata.rep{1,1}.phystime(1,1))*3;
        
        binsPerSec=size(htmp.Values,2)/DS.physDB(1).physdata.all.phystime(end,1);
        timeHist=(0:htmp.NumBins-1)/binsPerSec;
        timeHist=[0,timeHist];
        valueHist=[htmp.Values,0];
        plot(timeHist,valueHist);
        DS.valueHistAll{c,n}=valueHist;        % rows represents cells and columns represent direction
        tmp={};
    end   
end
close (ftmp);

%% extract the skipetimes at the pdir of DS cells and plot the mean and std line histogram
for c=1:size(DS.spiketimes,2)                 % cell#
    p=round(DS.pdir(c)/45)+1;
    DS.pResponse(c,:)=cell2mat(DS.valueHistAll(c,p));
end

for k=1:size(DS.pResponse,1)
DS.pResponseNor(k,:)=mat2gray(DS.pResponse(k,:));
end


DS.meanHist=mean(DS.pResponseNor,1);     % calculate mean lineHist across cells
DS.stdHist=std(DS.pResponseNor,1);       % calculate standard deviation of lineHist across cells
f1=figure;
%set(f1,'position',[435 50 1550 700]);
shadedErrorBar(timeHist,DS.meanHist,DS.stdHist,'k');
yl2(1,1:2)=ylim;
%ylim([0 yl2(2)]);
ylim([0 1]);
xlabel('time (seconds)','FontSize',18,'FontWeight','bold','Color','k');
ylabel('normalized frequency  (A.U.)','FontSize',18,'FontWeight','bold','Color','k');
% calculate and plot the onset and offset of stimulus        
load(physDBFileList{1,1});            
plotOnsetOffset(physdata,1);        % coresponds to direction of 0 degrees. This is used as approximation becase the mean preferred responses are from different 
                                    %directions and there is no accurate way to calculate the onset and offset of the stimulus to such a composite group. 

%% save files
save('physdataSummary.mat','DS','-v7.3');


end