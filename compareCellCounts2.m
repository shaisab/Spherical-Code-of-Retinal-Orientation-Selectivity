
function [] = compareCellCounts2()


% plot the contour of diff percentage between ONOFF and ON cells, for each of four translation channels
[FileName,PathName] = uigetfile('*.mat','Select a cellCountONandONOFF file');
load(FileName);
% pooledmap_ONOFF.cellCount=flipud(pooledmap_ONOFF.cellCount);
% pooledmap_ON.cellCount=flipud(pooledmap_ON.cellCount);

%% generate a single matrix for the countONOFF (all channels pooled) 
for a=1:10
    for b=1:10
countONOFF(a,b)=100*(sum(pooledmap_ONOFF.cellCount{a,b})./sum(sum([pooledmap_ONOFF.cellCount{:}])));       
    end
end
countONOFF(find(countONOFF==0))=nan;
%% generate a single matrix for the countONOFF (all channels pooled) 
for a=1:10
    for b=1:10
countON(a,b)=100*(sum(pooledmap_ON.cellCount{a,b})./sum(sum([pooledmap_ON.cellCount{:}]))); 
    end
end
countON(find(countON==0))=nan;
%% generate a single matrix for the diff count of ON and ONOFF cells
diffCount=10*(countONOFF-countON);      % multiply by 10 to roughly range between -40 and 40 

%% calculate percentages for each channel, for either ON or ONOFF
for i=1:size(pooledmap_ON.cellCount,1)
    for j=1:size(pooledmap_ON.cellCount,2)
        percentONOFF{i,j}=100*(pooledmap_ONOFF.cellCount{i,j}./sum(pooledmap_ONOFF.cellCount{i,j}));
        percentON{i,j}=100*(pooledmap_ON.cellCount{i,j}./sum(pooledmap_ON.cellCount{i,j}));
        diff{i,j}=percentONOFF{i,j}-percentON{i,j};         
    end
end

%% generate a single matrix for the percentONOFF for each channel 
for c=1:size(percentONOFF{1,1},2)
    fieldname=['c',num2str(c)];
    for i=1:size(percentONOFF,1)
        for j=1:size(percentONOFF,2)
            singlePercentONOFF.(fieldname)(i,j)=percentONOFF{i,j}(c);
        end
    end
end

%% generate a single matrix for the percentON for each channel 
for c=1:size(percentON{1,1},2)
    fieldname=['c',num2str(c)];
    for i=1:size(percentON,1)
        for j=1:size(percentON,2)
            singlePercentON.(fieldname)(i,j)=percentON{i,j}(c);
        end
    end
end
%% generate a single matrix for the diff percentage for each channel
for c=1:size(diff{1,1},2)
    fieldname=['c',num2str(c)];
    for i=1:size(diff,1)
        for j=1:size(diff,2)
            singleDiff.(fieldname)(i,j)=diff{i,j}(c);
        end
    end
end

%% plot channel distribution files

plotChannelDistribution(singleDiff,1,1,'channel distribution for diff percentages between ON and ONOFF DS cells')
plotChannelDistribution(singlePercentONOFF,0,1,'channel distribution for ONOFF DS cells')
plotChannelDistribution(singlePercentON,0,1,'channel distribution for ON DS cells')
% plotChannelDistribution(countON,0,0,'cell count for ON DS cells')
% plotChannelDistribution(countONOFF,0,0,'cell count for ON DS cells')
plotChannelDistribution(diffCount,1,0,'diff count between ON and ONOFF DS cells')

end
