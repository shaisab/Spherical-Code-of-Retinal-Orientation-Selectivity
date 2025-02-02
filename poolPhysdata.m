
function []= poolPhysdata()

filedir=cd;
physFileList=getFileList(filedir,['physdata_DB'],0,'anywhere');
        
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
             i=i+1;
        elseif strcmp(curname(inds(end)+1:length(curname)-4),'ONOFFDS')
            ONOFFDS.phys(j)=load(num2str(physFileList{ncell}));
            ONOFFDS.spiketimes{1,i}=ONOFFDS.phys(i).physdata.all.spiketimesAll;
            ONOFFDS.spikeArray{1,i}=ONOFFDS.phys(i).physdata.all.spikeArrayAll;
            j=j+1;
        end
    end
end

%% for each direction, pool all spiketimes across cells to be used to plot the PSTH
tmp={};
for n=1:size(ONDS.spiketimes{1,1},2)        % direction#
    for m=1:size(ONDS.spiketimes{1,1},1)    % repitition#
        for c=1:size(ONDS.spiketimes,2)     % cell#           
            tmp=[tmp,ONDS.spiketimes{1,c}{m,n}];  
        end
    end
    ONDS.pooled.spiketimes{n}=cell2mat(tmp);
    tmp={};
end

%% plot PSTH for ONDS across directions
for k=1:size(ONDS.pooled.spiketimes,2)                                   % loop on each of the 8 directions
    spiketimesDir=cell2mat(ONDS.pooled.spiketimes(:,k)');                % pool all spiketimes for all runs and reps
    spiketimesT=(spiketimesDir.*ONDS.phys(1).physdata.h(1).si)./1000000;          % spike times in seconds
    h(k)=subplot(3,3,k);
    hi=histogram(spiketimesT);                  % plot PSTH
    hi.BinEdges=0:1:floor(ONDS.phys(1).physdata.rep{1,1}.phystime(end,1)-ONDS.phys(1).physdata.rep{1,1}.phystime(1,1));
    hi.NumBins=floor(ONDS.phys(1).physdata.rep{1,1}.phystime(end,1)-ONDS.phys(1).physdata.rep{1,1}.phystime(1,1))*5;
    hi.FaceColor='k';
    hi.EdgeColor='k';
    title(['direction ',num2str(k*45-45)]);
    xlabel('time (seconds)');
    ylabel('number of spikes');
    hold on
    plot(ONDS.phys(1).physdata.all.phystime(:,1),3*ONDS.phys(1).physdata.rep{1,1}.signal{1,1}(:,1)-4,'k');
    yl(k,1:2)=ylim;
    hold on
    plotOnsetOffset(ONDS.phys(1).physdata,k);
end
%standardize the ylim of all subplots
for j=1:size(ONDS.pooled.spiketimes,2)
    ylim(h(j),[-5 max(yl(:,2))]);
end


%% plot line histogram

% binsPerSec=size(hi.Values,2)/ONDS.phys(1).physdata.all.phystime(end,1);
% timeHist=(0:hi.NumBins-1)/binsPerSec;
% timeHist=[0,timeHist];
% valueHist=[hi.Values,0];
% plot(timeHist,valueHist);


%% extract the bin values for each cell and each direction (repetitions pooled)
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
        htmp.NumBins=floor(ONDS.phys(1).physdata.rep{1,1}.phystime(end,1)-ONDS.phys(1).physdata.rep{1,1}.phystime(1,1))*1;
        
        binsPerSec=size(htmp.Values,2)/ONDS.phys(1).physdata.all.phystime(end,1);
        timeHist=(0:htmp.NumBins-1)/binsPerSec;
        timeHist=[0,timeHist];
        valueHist=[htmp.Values,0];
        plot(timeHist,valueHist);
        valueHistAll{c,n}=valueHist;        % rows represents cells and columns represent direction
        tmp={};
    end   
end
close (ftmp);




for n=1:size(ONDS.spiketimes{1,1},2)        % direction#
meanHist(n,:)=mean(cell2mat(valueHistAll(:,n)),1);     % calculate mean lineHist across cells
stdHist(n,:)=std(cell2mat(valueHistAll(:,n)),1);       % calculate standard deviation of lineHist across cells
end


f12=figure;
set(f12,'position',[5 150 950 800]);
for k=1:size(ONDS.pooled.spiketimes,2)                                   % loop on each of the 8 directions
    h(k)=subplot(3,3,k);
     plot(timeHist,meanHist);
    
end


% plots mean responses of all ONDS cells

subplot(2,2,1);

shadedErrorBar(1:size(DSnorpResponseT,1),mean(DSnorpResponseT(:,idxT==1),2),std(DSnorpResponseT(:,idxT==1),[],2),'k');
xlim([0 size(DSnorpResponseT,1)]); ylim([0 1]);
title(['n = ',num2str(sum(idxT==1))]);
subplot(2,2,4);
shadedErrorBar(1:size(DSnorpResponseT,1),mean(DSnorpResponseT(:,idxT==2),2),std(DSnorpResponseT(:,idxT==2),[],2),'k');
xlim([0 size(DSnorpResponseT,1)]); ylim([0 1]);
title(['n = ',num2str(sum(idxT==2))]);
subtitle('all DS cells');









TimeONDS(:,~any(TimeONDS,1))=[];  %columns
pResponseONDS(:,~any(pResponseONDS,1))=[];  %columns
nResponseONDS(:,~any(nResponseONDS,1))=[];  %columns
TimeONOFFDS(:,~any(TimeONOFFDS,1))=[];  %columns
pResponseONOFFDS(:,~any(pResponseONOFFDS,1))=[];  %columns
nResponseONOFFDS(:,~any(nResponseONOFFDS,1))=[];  %columns

for i=1:size(pResponseONDS,2)
    pResponseONDSS(:,i)=smooth(pResponseONDS(:,i),15,'moving');
    %pResponseONDSnorS(:,i)=mat2gray(pResponseONDSS(:,i));
    pResponseONDSnorS(:,i)=(pResponseONDSS(:,i)-mean(pResponseONDSS(10:40,i),1))./(max(pResponseONDSS(:,i))-mean(pResponseONDSS(10:40,i),1));
    pResponseONDSnor(:,i)=mat2gray(pResponseONDS(:,i));
    
    nResponseONDSnor(:,i)=mat2gray(nResponseONDS(:,i));
    
    [peakLocONmax{i}, peakMagONmax{i}]=peakfinder(pResponseONDSnorS(:,i),0.1,0.1,1,false);
    %[peakLocONmin{i}, peakMagONmin{i}]=peakfinder(pResponseONDSnorS(91:end,i),0.05,0.9,-1,false);
    peakMinMagONDS(i)=pResponseONDSnorS(130,i);     %data point number 130 is where the minima occur for both ON and ONOFF cells
end
for i=1:size(pResponseONOFFDS,2)
    pResponseONOFFDSS(:,i)=smooth(pResponseONOFFDS(:,i),15,'moving');
    %pResponseONOFFDSnorS(:,i)=mat2gray(pResponseONOFFDSS(:,i));
    pResponseONOFFDSnorS(:,i)=(pResponseONOFFDSS(:,i)-mean(pResponseONOFFDSS(10:40,i),1))./(max(pResponseONOFFDSS(:,i))-mean(pResponseONOFFDSS(10:40,i),1));
    pResponseONOFFDSnor(:,i)=mat2gray(pResponseONOFFDS(:,i));
    
    nResponseONOFFDSnor(:,i)=mat2gray(nResponseONOFFDS(:,i));
    
    [peakLocONOFFmax{i}, peakMagONOFFmax{i}]=peakfinder(pResponseONOFFDSnorS(:,i),0.1,0.1,1,false);
    %[peakLocONOFFmin{i}, peakMagONOFFmin{i}]=peakfinder(pResponseONOFFDSnorS(91:end,i),0.05,0.5,-1,false);
    peakMinMagONOFFDS(i)=pResponseONOFFDSnorS(130,i);
end


elementToExtract = 1;
% extract the chosen element of all the cells in peakLocON.
peakLocONmaxf=cellfun(@(x) x(elementToExtract),peakLocONmax);
peakLocONOFFmaxf=cellfun(@(x) x(elementToExtract),peakLocONOFFmax);
peakMagONmaxf=cellfun(@(x) x(elementToExtract),peakMagONmax);
peakMagONOFFmaxf=cellfun(@(x) x(elementToExtract),peakMagONOFFmax);

peakLocONmaxfSec=TimeONDS(peakLocONmaxf);
peakLocONOFFmaxfSec=TimeONOFFDS(peakLocONOFFmaxf);


yDataONDS=pResponseONDSnorS(peakLocONmaxf:peakLocONmaxf+14,:);
yDataONOFFDS=pResponseONOFFDSnorS(peakLocONOFFmaxf:peakLocONOFFmaxf+14,:);
% Set up fittype and options.
ft = fittype( 'poly1' );
opts = fitoptions( ft );
% Fit model to data.
for i=1:size(pResponseONDS,2)
[fitresult, gof] = fit((1:size(yDataONDS,1))', yDataONDS(:,i), ft, opts );
slopeONDS(i)=fitresult.p1;
end

for i=1:size(pResponseONOFFDS,2)
[fitresult, gof] = fit((1:size(yDataONOFFDS,1))', yDataONOFFDS(:,i), ft, opts );
slopeONOFFDS(i)=fitresult.p1;
end

% i=1
% [fitresult, gof] = fit((1:size(yDataONDS,1))', yDataONDS(:,i), ft, opts );
% % Plot fit with data.
% figure( 'Name', 'untitled fit ON' );
% h = plot( fitresult, (1:size(yDataONDS,1))', yDataONDS(:,i) );
% legend( h, 'tmpOO', 'untitled fit ON', 'Location', 'NorthEast' );
% % Label axes
% ylabel( 'tmpOO' );
% grid on



% elementToExtract = 1;
% % extract the chosen element of all the cells in peakLocON.
% peakLocONminf=cellfun(@(x) x(elementToExtract),peakLocONmin);
% peakLocONOFFminf=cellfun(@(x) x(elementToExtract),peakLocONOFFmin);
% peakMagONminf=cellfun(@(x) x(elementToExtract),peakMagONmin);
% peakMagONOFFminf=cellfun(@(x) x(elementToExtract),peakMagONOFFmin);

maxminRatioONDS=peakMagONmaxf./peakMinMagONDS;
maxminRatioONOFFDS=peakMagONOFFmaxf./peakMinMagONOFFDS;

sumONDS=sum(pResponseONDSnorS(60:130,:),1);
sumONOFFDS=sum(pResponseONOFFDSnorS(60:130,:),1);


%%%%% Confidence intervals %%%%%%%
statistic = @(x)(mean(x)); 
ci_maxminRatioONDS=bootci(10000,{statistic,maxminRatioONDS},'alpha',0.01,'type','bca'); 
ci_maxminRatioONOFFDS=bootci(10000,{statistic,maxminRatioONOFFDS},'alpha',0.01,'type','bca');

ci_peakLocONmaxfSec=bootci(10000,{statistic,peakLocONmaxfSec},'alpha',0.01,'type','bca'); 
ci_peakLocONOFFmaxfSec=bootci(10000,{statistic,peakLocONOFFmaxfSec},'alpha',0.01,'type','bca'); 

ci_slopeONDS=bootci(10000,{statistic,slopeONDS},'alpha',0.01,'type','bca');
ci_slopeONOFFDS=bootci(10000,{statistic,slopeONOFFDS},'alpha',0.01,'type','bca');

ci_sumONDS=bootci(10000,{statistic,sumONDS},'alpha',0.01,'type','bca');
ci_sumONOFFDS=bootci(10000,{statistic,sumONOFFDS},'alpha',0.01,'type','bca');

% ci_peakMinMagONDS=bootci(10000,{statistic,peakMinMagONDS},'alpha',0.01,'type','bca');
% ci_peakMinMagONOFFDS=bootci(10000,{statistic,peakMinMagONOFFDS},'alpha',0.01,'type','bca');

figure
scatter(slopeONDS,maxminRatioONDS,'r');
hold on
scatter(slopeONOFFDS,maxminRatioONOFFDS,'g');

figure
scatter(peakLocONmaxfSec,sumONDS,'r');
hold on
scatter(peakLocONOFFmaxfSec,sumONOFFDS,'g');

% plot responses at the preffered direction
f1=figure;
set(f1,'position',[435 50 1550 700]);
subplot(1,2,1)
shadedErrorBar(TimeONDS(:,1),mean(pResponseONDSnor,2),std(pResponseONDSnor,0,2),'k');
hold on
plot(TimeONDS,mean(pResponseONDSnor,2),'r','linewidth',2);
ylim([0 1]);
xlabel('Time (seconds)','FontSize',12,'FontWeight','bold','Color','k')
ylabel('\it\DeltaF/F  (A.U.)','FontSize',12,'FontWeight','bold','Color','k')
subplot(1,2,2)
try 
shadedErrorBar(TimeONOFFDS(:,1),mean(pResponseONOFFDSnor,2),std(pResponseONOFFDSnor,0,2),'k');
hold on
plot(TimeONOFFDS,mean(pResponseONOFFDSnor,2),'r','linewidth',2);
ylim([0 1]);
hold off
xlabel('Time (seconds)','FontSize',12,'FontWeight','bold','Color','k')
ylabel('\it\DeltaF/F  (A.U.)','FontSize',12,'FontWeight','bold','Color','k')
end

try
f2=figure;
set(f2,'position',[435 50 1050 900]);
plot(TimeONDS,mean(pResponseONDSnor,2),'r','linewidth',2);
hold on
plot(TimeONOFFDS,mean(pResponseONOFFDSnor,2),'g','linewidth',2);
ylim([0 1]);
hold off
end

% plot responses at the null direction
f3=figure;
set(f3,'position',[435 50 1550 700]);
subplot(1,2,1)
shadedErrorBar(TimeONDS(:,1),mean(nResponseONDSnor,2),std(nResponseONDSnor,0,2),'k');
hold on
plot(TimeONDS,mean(nResponseONDSnor,2),'r','linewidth',2);
ylim([0 1]);
xlabel('Time (seconds)','FontSize',12,'FontWeight','bold','Color','k')
ylabel('\it\DeltaF/F  (A.U.)','FontSize',12,'FontWeight','bold','Color','k')

subplot(1,2,2)
try
shadedErrorBar(TimeONOFFDS(:,1),mean(nResponseONOFFDSnor,2),std(nResponseONOFFDSnor,0,2),'k');
hold on
plot(TimeONOFFDS,mean(nResponseONOFFDSnor,2),'r','linewidth',2);
ylim([0 1]);
hold off
xlabel('Time (seconds)','FontSize',12,'FontWeight','bold','Color','k')
ylabel('\it\DeltaF/F  (A.U.)','FontSize',12,'FontWeight','bold','Color','k')
end

try
f4=figure;
set(f4,'position',[435 50 1050 900]);
plot(TimeONDS,mean(nResponseONDSnor,2),'r','linewidth',2);
hold on
plot(TimeONOFFDS,mean(nResponseONOFFDSnor,2),'g','linewidth',2);
ylim([0 1]);
hold off
end

summary.fileID=physFileList;
summary.TimeONDS=TimeONDS;
summary.pResponseONDS=pResponseONDS;
summary.pResponseONDSnor=pResponseONDSnor;
summary.pResponseONDSnorS=pResponseONDSnorS;
summary.pResponseONOFFDS=pResponseONOFFDS;
summary.pResponseONOFFDSnor=pResponseONOFFDSnor;
summary.pResponseONOFFDSnorS=pResponseONOFFDSnorS;
summary.TimeONOFFDS=TimeONOFFDS;
summary.nResponseONDS=nResponseONDS;
summary.nResponseONDSnor=nResponseONDSnor;
summary.nResponseONOFFDS=nResponseONOFFDS;
summary.nResponseONOFFDSnor=nResponseONOFFDSnor;
summary.peakLocONmaxfSec=peakLocONmaxfSec;
summary.peakLocONOFFmaxfSec=peakLocONOFFmaxfSec;
summary.peakMagONmaxf=peakMagONmaxf;
summary.peakMagONOFFmaxf=peakMagONOFFmaxf;
summary.slopeONDS=slopeONDS;
summary.slopeONOFFDS=slopeONOFFDS;
summary.maxminRatioONDS=maxminRatioONDS;
summary.maxminRatioONOFFDS=maxminRatioONOFFDS;
summary.sumONDS=sumONDS;
summary.sumONOFFDS=sumONOFFDS;
summary.ci_maxminRatioONDS=ci_maxminRatioONDS;
summary.ci_maxminRatioONOFFDS=ci_maxminRatioONOFFDS;
summary.ci_peakLocONmaxfSec=ci_peakLocONmaxfSec;
summary.ci_peakLocONOFFmaxfSec=ci_peakLocONOFFmaxfSec;
summary.ci_slopeONDS=ci_slopeONDS;
summary.ci_slopeONOFFDS=ci_slopeONOFFDS;
summary.ci_sumONDS=ci_sumONDS;
summary.ci_sumONOFFDS=ci_sumONOFFDS;


ind=strfind(filedir,'\');
filename=filedir;
filename=filename(ind(end)+1:length(filename));
save(['summary_',filename,'.mat'],'summary');
saveas(f1,['pResponse_ONDS_ONOFFDS_',filename],'fig');
saveas(f3,['nResponse_ONDS_ONOFFDS_',filename],'fig');
try saveas(f2,['pResponse_ONDS_ONOFFDS_overlay_',filename],'fig');end
try saveas(f4,['nResponse_ONDS_ONOFFDS_overlay_',filename],'fig');end


end