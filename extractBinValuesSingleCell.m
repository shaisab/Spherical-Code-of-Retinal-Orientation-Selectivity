
function DS = extractBinValuesSingleCell(DB,DG)
        
i=1;
if ~isempty(DB)
    for ncell=1:size(DB,2)
             DS.physDB(i)=DB;   
             DS.physDG(i)=DG;   
             DS.spiketimes{1,i}=DS.physDB(i).physdata.all.spiketimesAll;
             DS.spikeArray{1,i}=DS.physDB(i).physdata.all.spikeArrayAll;
             DS.pdir(1,i)=DS.physDG(i).physdata.all.pdir;
             i=i+1;
    end
end

%% extract bin values for each DS cell and each direction (repetitions pooled)
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
        htmp.NumBins=floor(DS.physDB(1).physdata.rep{1,1}.phystime(end,1)-DS.physDB(1).physdata.rep{1,1}.phystime(1,1))*2;
        
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
    if p==9; p=1; end
    DS.pResponse(c,:)=cell2mat(DS.valueHistAll(c,p));
end

end