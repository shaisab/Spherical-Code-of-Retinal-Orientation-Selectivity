
function []= poolPhysdata3()
%% extract data
filedir=cd;
physFileList=getFileList(filedir,'physdata_DB',0,'anywhere');     
i=1;
j=1;
if ~isempty(physFileList)
    for ncell=1:size(physFileList,2)
             DB.phys(i)=load(num2str(physFileList{ncell}));   
             DB.spiketimes{1,i}=DB.phys(i).physdata.all.spiketimesAll;
             DB.spikeArray{1,i}=DB.phys(i).physdata.all.spikeArrayAll;
             i=i+1;
    end
end

physFileList=getFileList(filedir,'physdata_DG',0,'anywhere');      
i=1;
j=1;
if ~isempty(physFileList)
    for ncell=1:size(physFileList,2)
             DG.phys(i)=load(num2str(physFileList{ncell}));   
             DG.pdir(1,i)=DG.phys(i).physdata.all.pdir;
             i=i+1;
    end
end


%% extract the bin values for each DS cell and each direction (repetitions pooled)
ftmp=figure;
tmp={};
for c=1:size(DB.spiketimes,2)                 % cell#
    for n=1:size(DB.spiketimes{1,1},2)        % direction#
        for m=1:size(DB.spiketimes{1,1},1)    % repitition#          
            tmp=[tmp,DB.spiketimes{1,c}{m,n}];       
        end
        spiketimesT=(cell2mat(tmp).*DB.phys(1).physdata.h(1).si)./1000000;          % spike times in seconds
        htmp=histogram(spiketimesT,'BinLimits',[0 floor(DB.phys(1).physdata.rep{1,1}.phystime(end)-DB.phys(1).physdata.rep{1,1}.phystime(1))],...
            'BinWidth',0.1);
        xlim([0 (DB.phys(1).physdata.rep{1,1}.phystime(end)-DB.phys(1).physdata.rep{1,1}.phystime(1))]);
        
        %htmp=histogram(spiketimesT,floor(DB.phys(1).physdata.rep{1,1}.phystime(end)-DB.phys(1).physdata.rep{1,1}.phystime(1))*5);                  % plot PSTH
        %htmp.NumBins=floor(DB.phys(1).physdata.rep{1,1}.phystime(end,1)-DB.phys(1).physdata.rep{1,1}.phystime(1,1))*2;
        %binsPerSec=10*(size(htmp.Values,2)/DB.phys(1).physdata.all.phystime(end,1));
        binsPerSec=1/htmp.BinWidth;
        timeHist=(0:1:htmp.NumBins-1)/binsPerSec;
        timeHist=[0,timeHist];
        valueHist=[htmp.Values,0];
        plot(timeHist,valueHist);
        DB.valueHistAll{c,n}=valueHist;        % rows represents cells and columns represent direction
        tmp={};
    end   
end
close (ftmp);

%% extract the skipetimes at the pdir of ONDS cells and plot the mean and std line histogram
for c=1:size(DB.spiketimes,2)                 % cell#
    p=round(DG.pdir(c)/45)+1;
    DB.pResponse(c,:)=cell2mat(DB.valueHistAll(c,p));
end

for k=1:size(DB.pResponse,1)
DB.pResponseNor(k,:)=mat2gray(DB.pResponse(k,:));
end

% tmp fix %%%%%%%%%%%%%%%%%%%%%%
DB.pResponseNor=DB.pResponse;           % plots the non-normalized responses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DB.meanHist=mean(DB.pResponseNor,1);     % calculate mean lineHist across cells
DB.stdHist=std(DB.pResponseNor,1);       % calculate standard deviation of lineHist across cells
DB.seHist=std(DB.pResponseNor,1)./sqrt(size(DB.pResponseNor,1));       % calculate standard deviation of lineHist across cells
f1=figure;
set(f1, 'color', [1 1 1]);
hold on;
% set(f1,'position',[435 50 1550 700]);
% subplot(1,2,1);
shadedErrorBar(timeHist,DB.meanHist,DB.seHist,'k');
plot(timeHist,DB.meanHist,'k','LineWidth',1);
yl2(1,1:2)=ylim;
set(gca,'FontSize',14);
%ylim([0 yl2(2)]);
ylim([0 5]);
xlim([4 16]);
xlabel('Time (sec)','FontSize',18,'Color','k');
%ylabel('normalized frequency  (A.U.)','FontSize',18,'FontWeight','bold','Color','k');
ylabel('Firing rate (spikes 100 msec^-^1)','FontSize',18,'Color','k');
% calculate and plot the onset and offset of stimulus        
load(physFileList{1,1});            
plotOnsetOffset(physdata,1);        % coresponds to direction of 0 degrees. This is used as approximation because the mean preferred responses are from different 
                                    %directions and there is not accurate way to calculate the onset and offset of the stimulus to such a composite group. 
box off
%% save files
save('physdataSummary.mat','DB','DG');


end