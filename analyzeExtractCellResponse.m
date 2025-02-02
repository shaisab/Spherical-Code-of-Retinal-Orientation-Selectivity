
function []= analyzeExtractCellResponse()


filedir=cd;
physFileList=getFileList(filedir,['physdata_DB'],0,'anywhere');         % used for plotting the onset and offset of stimulus

[filedir,go]=bbdirselector('select data folder',cd);
if ~go disp('no folder selected'); return; end

names=getFileList(filedir,'identifiedCell',0,'anywhere');

TimeONDS=zeros(300,300);
pResponseONDS=zeros(300,300);
nResponseONDS=zeros(300,300);
TimeONOFFDS=zeros(300,300);
pResponseONOFFDS=zeros(300,300);
nResponseONOFFDS=zeros(300,300);

if ~isempty(names)
    for ncell=1:size(names,2)
        curname=cell2mat(names(ncell));
        inds=find(curname=='_');
        if strcmp(curname(inds(end)+1:length(curname)-4),'ONDS')
            ONDSres=load(num2str(names{ncell}));
            
            dif=size(TimeONDS,1)-size(ONDSres.Time,1);
            if dif>=0
                TimeONDS(size(ONDSres.Time,1)+1:end,:)=[];
                pResponseONDS(size(ONDSres.pResponse,1)+1:end,:)=[];
                nResponseONDS(size(ONDSres.nResponse,1)+1:end,:)=[];
            elseif dif<0
                ONDSres.Time(size(TimeONDS,1)+1:end,:)=[];
                ONDSres.pResponse(size(pResponseONDS,1)+1:end,:)=[];
                ONDSres.nResponse(size(nResponseONDS,1)+1:end,:)=[];
            end
            
            TimeONDS(:,ncell)=ONDSres.Time;
            pResponseONDS(:,ncell)=ONDSres.pResponse;
            nResponseONDS(:,ncell)=ONDSres.nResponse;
            
        elseif strcmp(curname(inds(end)+1:length(curname)-4),'ONOFFDS')
            ONOFFDSres=load(num2str(names{ncell}));
            
            dif=size(TimeONOFFDS,1)-size(ONOFFDSres.Time,1);
            if dif>=0
                TimeONOFFDS(size(ONOFFDSres.Time,1)+1:end,:)=[];
                pResponseONOFFDS(size(ONOFFDSres.pResponse,1)+1:end,:)=[];
                nResponseONOFFDS(size(ONOFFDSres.nResponse,1)+1:end,:)=[];
            elseif dif<0
                ONOFFDSres.Time(size(TimeONOFFDS,1)+1:end,:)=[];
                ONOFFDSres.pResponse(size(pResponseONOFFDS,1)+1:end,:)=[];
                ONOFFDSres.nResponse(size(nResponseONOFFDS,1)+1:end,:)=[];
            end
            
            TimeONOFFDS(:,ncell)=ONOFFDSres.Time;
            pResponseONOFFDS(:,ncell)=ONOFFDSres.pResponse;
            nResponseONOFFDS(:,ncell)=ONOFFDSres.nResponse;
        end
        
    end
end

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
shadedErrorBar(TimeONDS(:,1)+1,mean(pResponseONDSnor,2),std(pResponseONDSnor,0,2),'k');
hold on
plot(TimeONDS+1,mean(pResponseONDSnor,2),'r','linewidth',2);
ylim([0 1]);
xlim([0 20]);
xlabel('time (seconds)','FontSize',18,'FontWeight','bold','Color','k')
ylabel('\it\DeltaF/F  (A.U.)','FontSize',18,'FontWeight','bold','Color','k')
% calculate and plot the onset and offset of stimulus        
load(physFileList{1,1});            
plotOnsetOffset(physdata,1);        % coresponds to direction of 0 degrees. This is used as approximation becase the mean preferred responses are from different 
                                    %directions and there is not accurate way to calculate the onset and offset of the stimulus to such a composite group. 
subplot(1,2,2)
try 
shadedErrorBar(TimeONOFFDS(:,1)+1,mean(pResponseONOFFDSnor,2),std(pResponseONOFFDSnor,0,2),'k');
hold on
plot(TimeONOFFDS+1,mean(pResponseONOFFDSnor,2),'r','linewidth',2);
ylim([0 1]);
xlim([0 20]);
hold off
xlabel('time (seconds)','FontSize',18,'FontWeight','bold','Color','k')
%ylabel('\it\DeltaF/F  (A.U.)','FontSize',12,'FontWeight','bold','Color','k')
% calculate and plot the onset and offset of stimulus        
load(physFileList{1,1});            
plotOnsetOffset(physdata,1);        % coresponds to direction of 0 degrees. This is used as approximation becase the mean preferred responses are from different 
                                    %directions and there is not accurate way to calculate the onset and offset of the stimulus to such a composite group. 
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

summary.fileID=names;
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