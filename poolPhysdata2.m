
function []= poolPhysdata2()

filedir=cd;
physFileList=getFileList(filedir,'physdata_DB',0,'anywhere');
        
i=1;
j=1;
if ~isempty(physFileList)
    for ncell=1:size(physFileList,2)
        curname=physFileList{ncell};
        inds=find(curname=='_');
        if strcmp(curname(inds(end)+1:length(curname)-4),'ONDS') 
             ONDS.phys(i)=load(num2str(physFileList{ncell}));   
             ONDS.spiketimes{1,i}=ONDS.phys(i).physdata.all.spiketimesAll;
             ONDS.spikeArray{1,i}=ONDS.phys(i).physdata.all.spikeArrayAll;
             ONDS.pdir(1,i)=ONDS.phys(i).physdata.pdir;
             i=i+1;
        elseif strcmp(curname(inds(end)+1:length(curname)-4),'ONOFFDS')
            ONOFFDS.phys(j)=load(num2str(physFileList{ncell}));
            ONOFFDS.spiketimes{1,j}=ONOFFDS.phys(j).physdata.all.spiketimesAll;
            ONOFFDS.spikeArray{1,j}=ONOFFDS.phys(j).physdata.all.spikeArrayAll;
            ONOFFDS.pdir(1,j)=ONOFFDS.phys(j).physdata.pdir;
            j=j+1;
        end
    end
end


%% extract the bin values for each ONDS cell and each direction (repetitions pooled)
ftmp=figure;
tmp={};
for c=1:size(ONDS.spiketimes,2)                 % cell#
    for n=1:size(ONDS.spiketimes{1,1},2)        % direction#
        for m=1:size(ONDS.spiketimes{1,1},1)    % repitition#          
            tmp=[tmp,ONDS.spiketimes{1,c}{m,n}];       
        end
        spiketimesT=(cell2mat(tmp).*ONDS.phys(1).physdata.h(1).si)./1000000;          % spike times in seconds
        htmp=histogram(spiketimesT);                  % plot PSTH
        htmp.BinEdges=0:1:floor(ONDS.phys(1).physdata.rep{1,1}.phystime(end,1)-ONDS.phys(1).physdata.rep{1,1}.phystime(1,1));
        htmp.NumBins=floor(ONDS.phys(1).physdata.rep{1,1}.phystime(end,1)-ONDS.phys(1).physdata.rep{1,1}.phystime(1,1))*2;
        
        binsPerSec=size(htmp.Values,2)/ONDS.phys(1).physdata.all.phystime(end,1);
        timeHist=(0:htmp.NumBins-1)/binsPerSec;
        timeHist=[0,timeHist];
        valueHist=[htmp.Values,0];
        plot(timeHist,valueHist);
        ONDS.valueHistAll{c,n}=valueHist;        % rows represents cells and columns represent direction
        tmp={};
    end   
end
close (ftmp);

%% extract the skipetimes at the pdir of ONDS cells and plot the mean and std line histogram
for c=1:size(ONDS.spiketimes,2)                 % cell#
    p=round(ONDS.pdir(c)/45)+1;
    ONDS.pResponse(c,:)=cell2mat(ONDS.valueHistAll(c,p));
end

for k=1:size(ONDS.pResponse,1)
ONDS.pResponseNor(k,:)=mat2gray(ONDS.pResponse(k,:));
end


ONDS.meanHist=mean(ONDS.pResponseNor,1);     % calculate mean lineHist across cells
ONDS.stdHist=std(ONDS.pResponseNor,1);       % calculate standard deviation of lineHist across cells
f1=figure;
set(f1,'position',[435 50 1550 700]);
subplot(1,2,1);
shadedErrorBar(timeHist,ONDS.meanHist,ONDS.stdHist,'k');
yl2(1,1:2)=ylim;
%ylim([0 yl2(2)]);
ylim([0 1]);
xlabel('time (seconds)','FontSize',18,'FontWeight','bold','Color','k');
ylabel('normalized frequency  (A.U.)','FontSize',18,'FontWeight','bold','Color','k');
% calculate and plot the onset and offset of stimulus        
load(physFileList{1,1});            
plotOnsetOffset(physdata,1);        % coresponds to direction of 0 degrees. This is used as approximation becase the mean preferred responses are from different 
                                    %directions and there is not accurate way to calculate the onset and offset of the stimulus to such a composite group. 
%% extract the bin values for each ONOFFDS cell and each direction (repetitions pooled)
ftmp=figure;
tmp={};
for c=1:size(ONOFFDS.spiketimes,2)                 % cell#
    for n=1:size(ONOFFDS.spiketimes{1,1},2)        % direction#
        for m=1:size(ONOFFDS.spiketimes{1,1},1)    % repitition#          
            tmp=[tmp,ONOFFDS.spiketimes{1,c}{m,n}];       
        end
        spiketimesT=(cell2mat(tmp).*ONOFFDS.phys(1).physdata.h(1).si)./1000000;          % spike times in seconds
        htmp=histogram(spiketimesT);                  % plot PSTH
        htmp.BinEdges=0:1:floor(ONOFFDS.phys(1).physdata.rep{1,1}.phystime(end,1)-ONOFFDS.phys(1).physdata.rep{1,1}.phystime(1,1));
        htmp.NumBins=floor(ONOFFDS.phys(1).physdata.rep{1,1}.phystime(end,1)-ONOFFDS.phys(1).physdata.rep{1,1}.phystime(1,1))*2;
        
        binsPerSec=size(htmp.Values,2)/ONOFFDS.phys(1).physdata.all.phystime(end,1);
        timeHist=(0:htmp.NumBins-1)/binsPerSec;
        timeHist=[0,timeHist];
        valueHist=[htmp.Values,0];
        plot(timeHist,valueHist);
        ONOFFDS.valueHistAll{c,n}=valueHist;        % rows represents cells and columns represent direction
        tmp={};
    end   
end
close (ftmp);

%% extract the skipetimes at the pdir of ONOFFDS cells and plot the mean and std line histogram
for c=1:size(ONOFFDS.spiketimes,2)                 % cell#
    p=round(ONOFFDS.pdir(c)/45)+1;
    ONOFFDS.pResponse(c,:)=cell2mat(ONOFFDS.valueHistAll(c,p));
end

for k=1:size(ONOFFDS.pResponse,1)
ONOFFDS.pResponseNor(k,:)=mat2gray(ONOFFDS.pResponse(k,:));
end


ONOFFDS.meanHist=mean(ONOFFDS.pResponseNor,1);     % calculate mean lineHist across cells
ONOFFDS.stdHist=std(ONOFFDS.pResponseNor,1);       % calculate standard deviation of lineHist across cells
subplot(1,2,2);
if size(ONOFFDS.pResponse,1)==1
  plot(timeHist,ONOFFDS.meanHist,'k');  
else
shadedErrorBar(timeHist,ONOFFDS.meanHist,ONOFFDS.stdHist,'k');
end
yl3(1,1:2)=ylim;
%ylim([0 max([yl2(2),yl3(2)])]);
ylim([0 1]);
xlabel('time (seconds)','FontSize',18,'FontWeight','bold','Color','k')
%ylabel('normalized frequency  (A.U.)','FontSize',12,'FontWeight','bold','Color','k')
% calculate and plot the onset and offset of stimulus        
load(physFileList{1,1});            
plotOnsetOffset(physdata,1);        % coresponds to direction of 0 degrees. This is used as approximation becase the mean preferred responses are from different 
                                    %directions and there is not accurate way to calculate the onset and offset of the stimulus to such a composite group. 

%% save files
save('physdataSummary.mat','ONDS','ONOFFDS');


end