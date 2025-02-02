% Active directory should be the Tstacks or map folder, where all the mergedmeta files are located.


function []= plotMergedmetaCells2() 

thisdir = dir;
str = {thisdir.name};
[s,v] = listdlg('PromptString','Select mergedmeta file :',...
    'SelectionMode','multiple',...
    'ListString',str);
if ~v disp('no files selected'); return; end

for i=1:size(s,2)
    [mergedmetaFileList{i}]=thisdir(s(i)).name;
end

for filen=1:size(mergedmetaFileList,2)
    load(mergedmetaFileList{filen});
      
    cellType={'DS','OS'};
    for j=1:size(cellType,2)
        type=cellType{j};
        switch type
            case 'DS'
                property=mergedmeta.isDS;
            case 'OS'
                property=mergedmeta.isOS;
        end
        
        
        %     if strcmp('DS',DSorOS)
        %         property=mergedmeta.isDS;
        %     elseif strcmp('OS',DSorOS)
        %         property=mergedmeta.isOS;
        %     end
        
        %rep=1;
        %     filedir='E:\Calcium imaging data\Selected identified cells';
        %     SNSnames=getFileList(filedir,'summary_Selected identified cells',0,'anywhere');
        filedir='D:\OS data\All SNS cells';
        SNSnames=getFileList(filedir,'summary_All SNS cells',0,'anywhere');
        if ~isempty(SNSnames)
            SNS=load(num2str(SNSnames{1}));
        end
        
        for ROIn=1:size(mergedmeta.meanRepDG.ROIdata,2)
            
            if property(ROIn)
                % Plot DG
                r=mergedmeta.meanRepDG.ROIdata{ROIn}.summary.meannororddff;
                repDG=size(mergedmeta.metadataDG,2);
                a1=mergedmeta.meanRepDG.ROIdata{ROIn}.summary.dir;
                d=a1;
                a=mergedmeta.meanRepDG.ROIdata{ROIn}.summary.meanr;
                pdir=mergedmeta.meanRepDG.ROIdata{ROIn}.summary.pdir;
                sd=mergedmeta.meanRepDG.ROIdata{ROIn}.summary.SDpdir;
                DSI=mergedmeta.meanRepDG.ROIdata{ROIn}.summary.DSI;
                DI=mergedmeta.meanRepDG.ROIdata{ROIn}.summary.DI;
                pori=mergedmeta.meanRepDG.ROIdata{ROIn}.summary.pori;
                OSI=mergedmeta.meanRepDG.ROIdata{ROIn}.summary.OSI;
                globalOSI=mergedmeta.meanRepDG.ROIdata{ROIn}.summary.globalOSI;
                oriR2=mergedmeta.meanRepDG.ROIdata{ROIn}.summary.oriR2;
                SDpOri=mergedmeta.meanRepDG.ROIdata{ROIn}.summary.SDmaxori;
                sumCorrRep2=mergedmeta.meanRepDG.ROIdata{ROIn}.summary.sumCorrRep2;
                symmetry_ratio=mergedmeta.meanRepDG.ROIdata{ROIn}.summary.symmetry_ratio;
                try original_rmse=mergedmeta.meanRepDG.ROIdata{ROIn}.summary.original_rmse; end
                try shuffle_rmse=mergedmeta.meanRepDG.ROIdata{ROIn}.summary.shuffle_rmse; end
                
                dispfigh=figure(1000+ROIn);clf;
                set(dispfigh,'position',[5 150 950 800]);
                for p=1:size(r,2)
                    subplot(ceil((size(r,2)+1)/3),3,p);
                    for nrep=1:repDG
                        for nruns=1:mergedmeta.metadataDG{nrep}.log.numberofruns
                            plot(mergedmeta.metadataDG{1,nrep}.ROIdata{1,ROIn}.ordfluotime{1,nruns},mergedmeta.metadataDG{1,nrep}.ROIdata{1,ROIn}.nororddff{1,nruns}(:,p),'Color',[0.7,0.7,0.7],'linewidth',1);
                            hold on
                        end
                    end
                    plot(mergedmeta.meanRepDG.ROIdata{1,1}.summary.ordfluotime,r(:,p),'r','linewidth',2);
                    title(strcat(num2str(a1(p),'%6.0f'), ' degrees'));
                    maxy=max(max(r(:,1:size(r,2))));miny=min(min(r(:,1:size(r,2))));
                    for nrep=1:repDG
                        for nruns=1:mergedmeta.metadataDG{nrep}.log.numberofruns
                            maxym(nrep*nruns)=max(max(mergedmeta.metadataDG{1,nrep}.ROIdata{1,ROIn}.nororddff{1,nruns}(:,1:size(r,2))));
                            minym(nrep*nruns)=min(min(mergedmeta.metadataDG{1,nrep}.ROIdata{1,ROIn}.nororddff{1,nruns}(:,1:size(r,2))));
                        end
                    end
                    maxym=max(maxym);
                    minym=min(minym);
                    minF=min([miny minym]);
                    maxF=max([maxy maxym]);
                    ylim([minF maxF]);
                    xlim([0 max(mergedmeta.meanRepDG.ROIdata{1,1}.summary.ordfluotime)]);
                end
                
                subplot(ceil((size(r,2)+1)/3),3,p+1);
                dpm=[d,d(1)]; apm=[a,a(1)];
                maxap=0;
                for nrep=1:repDG
                    for nruns=1:mergedmeta.metadataDG{nrep}.log.numberofruns
                        ap=[mergedmeta.metadataDG{1,nrep}.ROIdata{1,ROIn}.meanr{1,nruns},mergedmeta.metadataDG{1,nrep}.ROIdata{1,ROIn}.meanr{1,nruns}(1)];
                        maxap=max([maxap apm ap]);
                    end
                end
                
                for nrep=1:repDG
                    for nruns=1:mergedmeta.metadataDG{nrep}.log.numberofruns
                        ap=[mergedmeta.metadataDG{1,nrep}.ROIdata{1,ROIn}.meanr{1,nruns},mergedmeta.metadataDG{1,nrep}.ROIdata{1,ROIn}.meanr{1,nruns}(1)];
                        P=polarSS2p(deg2rad(dpm), ap/max(ap),1,'b',1);
                        hold on
                    end
                end
                P=polarSS2p(deg2rad(dpm), apm/max(apm),1,'b',2);
                hold on;
                [X,Y]=pol2cart(deg2rad(pdir), 1);
                compassSS2p(X,Y,1,'r',2);
                [X,Y]=pol2cart(deg2rad(pori), 1);
                compassSS2p(X,Y,1,'g',2);
                
                
                for nrep=1:repDG
                    for nruns=1:mergedmeta.metadataDG{nrep}.log.numberofruns
                        hold on
                        [X,Y]=pol2cart(deg2rad(mergedmeta.metadataDG{1,nrep}.ROIdata{1,ROIn}.pdir{1,nruns}), 1);
                        compassSS2p(X,Y,1,'r',1);
                    end
                end
                
                %title(strcat('n = ',num2str(nr,'%6.0f'),',pDir = ',num2str(pdir,'%6.0f'),',SDpDir = ',num2str(sd,'%6.0f'), ', DSI = ',num2str(DSI,'%6.3f'), ', DI = ',num2str(DI,'%6.3f')));
                annotation('textbox', [0.55, 0.27, 0.1, 0.1], 'string', ['n = ',num2str(mergedmeta.meanRepDG.ROIdata{ROIn}.summary.n,'%6.0f'),', pDir = ',num2str(pdir,'%6.0f'),', SDpDir = ',num2str(sd,'%6.0f'),', pOri = ',num2str(pori,'%6.0f'),', oriR2 = ',num2str(oriR2,'%6.3f'),', SDpOri = ',num2str(SDpOri,'%6.0f')],'LineStyle','none')
                annotation('textbox', [0.6, 0, 0.1, 0.1], 'string', ['DSI = ',num2str(DSI,'%6.3f'), ', DI = ',num2str(DI,'%6.3f'), ', OSI = ',num2str(OSI,'%6.3f'), ', globalOSI = ',num2str(globalOSI,'%6.3f')],'LineStyle','none')
                annotation('textbox', [0.2, 0, 0.7, 0.05], 'string', ['isDS = ',num2str(mergedmeta.meanRepDG.ROIdata{ROIn}.summary.isDS,'%d'), ', isOS = ',num2str(mergedmeta.meanRepDG.ROIdata{ROIn}.summary.isOS,'%d'), ', ONorONOFForOFF = ',num2str(mergedmeta.meanRepDG.ROIdata{ROIn}.summary.ONorONOFForOFF,'%d')],'LineStyle','none','FontSize',12,'Color','r');
                
                
                
                %             annotation('textbox', [0.6, 0, 0.1, 0.1], 'string', ['DSI = ',num2str(DSI,'%6.3f'), ', DI = ',num2str(DI,'%6.3f'), ', OSI = ',num2str(OSI,'%6.3f')],'LineStyle','none')
                
                hold off;
                drawnow nocallbacks
              
                
                
                
                
                % Plot DB
                r=mergedmeta.meanRepDB.ROIdata{ROIn}.summary.meannororddff;
                repDB=size(mergedmeta.metadataDB,2);
                a1=mergedmeta.meanRepDB.ROIdata{ROIn}.summary.dir;
                d=a1;
                a=mergedmeta.meanRepDB.ROIdata{ROIn}.summary.meanr;
                pdir=mergedmeta.meanRepDB.ROIdata{ROIn}.summary.pdir;
                sd=mergedmeta.meanRepDB.ROIdata{ROIn}.summary.SDpdir;
                DSI=mergedmeta.meanRepDB.ROIdata{ROIn}.summary.DSI;
                DI=mergedmeta.meanRepDB.ROIdata{ROIn}.summary.DI;
                pori=mergedmeta.meanRepDB.ROIdata{ROIn}.summary.pori;
                OSI=mergedmeta.meanRepDB.ROIdata{ROIn}.summary.OSI;
                globalOSI=mergedmeta.meanRepDB.ROIdata{ROIn}.summary.globalOSI;
                oriR2=mergedmeta.meanRepDB.ROIdata{ROIn}.summary.oriR2;
                SDpOri=mergedmeta.meanRepDB.ROIdata{ROIn}.summary.SDmaxori;
                sumCorrRep2=mergedmeta.meanRepDB.ROIdata{ROIn}.summary.sumCorrRep2;
                symmetry_ratio=mergedmeta.meanRepDB.ROIdata{ROIn}.summary.symmetry_ratio;
                try original_rmse=mergedmeta.meanRepDB.ROIdata{ROIn}.summary.original_rmse; end
                try shuffle_rmse=mergedmeta.meanRepDB.ROIdata{ROIn}.summary.shuffle_rmse; end
                
                pR=strfind(mergedmeta.metadataDB{1,1}.imheader.Footnote.Resolution,'_');
                pixperum=str2num(mergedmeta.metadataDB{1,1}.imheader.Footnote.Resolution(pR(1)+1:pR(2)-1));
                widthum=mergedmeta.metadataDB{1,1}.imheader.General.SizeX/pixperum;
                heightum=mergedmeta.metadataDB{1,1}.imheader.General.SizeY/pixperum;
                roiLocumX=mergedmeta.metadataDB{repDB}.ROIdata{ROIn}.centerX/pixperum;
                roiLocumY=mergedmeta.metadataDB{repDB}.ROIdata{ROIn}.centerY/pixperum;
                
                
                
                
                
                
                
                %conditions for cell type identification (ONDS vs ONOFFDS based CI of SNS cells)
                %%%%%%%%%%%%%%%%%%%%%
                if  ~isfield(mergedmeta.meanRepDB.ROIdata{ROIn}.summary,'isONDS')
                    
                    if  mergedmeta.meanRepDG.ROIdata{ROIn}.summary.isDS==1
                     %if  mergedmeta.isDS==1     % for LAP4 treatments only
                        [~,p]=max(mergedmeta.meanRepDG.ROIdata{1,ROIn}.summary.meanr);
                        n=p+4;
                        if n>8
                            n=p-4;
                        end
                        pResponse=mergedmeta.meanRepDB.ROIdata{1,ROIn}.summary.meannororddff(:,p);
                    elseif mergedmeta.meanRepDG.ROIdata{ROIn}.summary.isOS==1
                        [~,p]=max(mergedmeta.meanRepDG.ROIdata{1,ROIn}.summary.meanr);
                        p2=p+4;
                        if p2>8
                            p2=p-4;
                        end
                        n=p+2;
                        if n>8
                            n=p-2;
                        end
                        pResponse=mean([mergedmeta.meanRepDB.ROIdata{1,ROIn}.summary.meannororddff(:,p),mergedmeta.meanRepDB.ROIdata{1,ROIn}.summary.meannororddff(:,p2)],2);
                    end
                    
                    
                    Time=mergedmeta.meanRepDB.ROIdata{1,ROIn}.summary.ordfluotime;
                    %figure; plot(Time, pResponse,'b');
                    
                    %pResponseNor=mat2gray(pResponse);
                    %pResponseNorS=smooth(pResponseNor,9,'moving');
                    
                    pResponseS=smooth(pResponse,15,'moving');
                    %pResponseS2=smooth(pResponse,5,'moving');
                    %pResponseNorS=mat2gray(pResponseS);
                    pResponseNorS=(pResponseS-mean(pResponseS(10:40,1),1))./(max(pResponseS)-mean(pResponseS(10:40,1),1));
                    %pResponseNorS2=(pResponseS2-mean(pResponseS2(10:40,1),1))./(max(pResponseS2)-mean(pResponseS2(10:40,1),1));
                    
                    peakMinMag=pResponseNorS(130,1);     %data point number 130 is where the minima occur for both ON and ONOFF cells
                    
                    dispfigh=figure(3000+ROIn);clf;
                    screendims=get(0,'MonitorPositions');
                    if numel(screendims)==8
                        set(dispfigh,'position',[5 150 950 800]);
                    elseif numel(screendims)==4
                        set(dispfigh,'position',[5 150 950 800]);
                    end
                        hold on;
                    peakfinder(pResponseNorS,0.1,0.1,1,false);
                    drawnow nocallbacks
                   
                    %peakfinder(pResponseNorS2,0.1,0.1,1,false);
                    [peakLocMax,peakMagMax]=peakfinder(pResponseNorS,0.1,0.1,1,false);
                    if sum((peakLocMax>145).*(peakLocMax<195),1)~=0
                        offPeak=1;
                    else
                        offPeak=0;
                    end
                    
                    %%% calculate slope of OFF peak
                    try peakLocOff=peakLocMax(end); end
                    
                    try
                        offind=find(peakLocMax>145);
                        peakLocOff=peakLocMax(offind(1),1);
                    end
                    try
                    yOffData=pResponseNorS(peakLocOff:peakLocOff+15,:);
                    % Set up fittype and options.
                    ft = fittype('poly1');
                    opts = fitoptions(ft);
                    % Fit model to data.
                    for i=1:size(pResponse,2)
                        [fitresultOff, gofOff] = fit((1:size(yOffData,1))', yOffData(:,i), ft, opts );
                        slopeOff(i)=fitresultOff.p1;
                    end
                    catch
                       slopeOff(i)=0; 
                    end
                    
                    
                    %%% calcultes slope of ON peak
                    peakMagMax(peakLocMax<60)=[];
                    peakLocMax(peakLocMax<60)=[];
                    try
                    peakLocMaxf=peakLocMax(1);
                    peakMagMaxf=peakMagMax(1);
                    peakLocMaxfSec=Time(peakLocMaxf);
                    yData=pResponseNorS(peakLocMaxf:peakLocMaxf+14,:);
                    
                    % Set up fittype and options.
                    ft = fittype('poly1');
                    opts = fitoptions(ft);
                    % Fit model to data.
                    for i=1:size(pResponse,2)
                        [fitresult, gof] = fit((1:size(yData,1))', yData(:,i), ft, opts );
                        slope(i)=fitresult.p1;
                    end
                    catch
                        slope(i)=0;
                    end
                    
                    
                    maxminRatio=peakMagMaxf/peakMinMag;
                    sumResponse=sum(pResponseNorS(60:130,:),1);
                    
                    
                    %if maxminRatio>SNS.summary.ci_maxminRatioONDS(1) && maxminRatio<SNS.summary.ci_maxminRatioONDS(2)
                    if maxminRatio<SNS.summary.ci_maxminRatioONDS(2) || peakLocMaxfSec>8
                        isONDS(1,1)=1;
                        isONOFFDS(1,1)=0;
                        %elseif maxminRatio>SNS.summary.ci_maxminRatioONOFFDS(1) && maxminRatio<SNS.summary.ci_maxminRatioONOFFDS(2)
                    elseif maxminRatio>SNS.summary.ci_maxminRatioONOFFDS(1)
                        isONDS(1,1)=0;
                        isONOFFDS(1,1)=1;
                    else
                        isONDS(1,1)=0;
                        isONOFFDS(1,1)=0;
                    end
                    
                    %if peakLocMaxfSec>SNS.summary.ci_peakLocONmaxfSec(1) && peakLocMaxfSec<SNS.summary.ci_peakLocONmaxfSec(2)
                    if peakLocMaxfSec>SNS.summary.ci_peakLocONmaxfSec(1)
                        isONDS(2,1)=1;
                        isONOFFDS(2,1)=0;
                        %elseif peakLocMaxfSec>SNS.summary.ci_peakLocONOFFmaxfSec(1) && peakLocMaxfSec<SNS.summary.ci_peakLocONOFFmaxfSec(2)
                    elseif peakLocMaxfSec<SNS.summary.ci_peakLocONOFFmaxfSec(2)
                        isONDS(2,1)=0;
                        isONOFFDS(2,1)=1;
                    else
                        isONDS(2,1)=0;
                        isONOFFDS(2,1)=0;
                    end
                    
                    %if peakMinMag>SNS.summary.ci_peakMinMagONDS(1) && peakMinMag<SNS.summary.ci_peakMinMagONDS(2)
                    if slope>SNS.summary.ci_slopeONDS(1)
                        isONDS(3,1)=1;
                        isONOFFDS(3,1)=0;
                        %elseif peakMinMag>SNS.summary.ci_peakMinMagONOFFDS(1) && peakMinMag<SNS.summary.ci_peakMinMagONOFFDS(2)
                    elseif slope<SNS.summary.ci_slopeONOFFDS(2)
                        isONDS(3,1)=0;
                        isONOFFDS(3,1)=1;
                    else
                        isONDS(3,1)=0;
                        isONOFFDS(3,1)=0;
                    end
                    
                    %if sumResponse>SNS.summary.ci_sumONDS(1) && sumResponse<SNS.summary.ci_sumONDS(2)
                    if sumResponse>SNS.summary.ci_sumONDS(1)
                        isONDS(4,1)=1;
                        isONOFFDS(4,1)=0;
                        %elseif sumResponse>SNS.summary.ci_sumONOFFDS(1) && sumResponse<SNS.summary.ci_sumONOFFDS(2)
                    elseif sumResponse<SNS.summary.ci_sumONOFFDS(2)
                        isONDS(4,1)=0;
                        isONOFFDS(4,1)=1;
                    else
                        isONDS(4,1)=0;
                        isONOFFDS(4,1)=0;
                    end
                    
                    %                      % excludes from analysis every cell with a weak
                    %                      % and noisy response (more than 4 detected peaks)
                    %                      if size(peakLocMax,1)>4
                    %                      isONDS(1:4,1)=0;
                    %                      isONOFFDS(1:4,1)=0;
                    %                      end
                    
                    
                    %remove the slope criterion
                    isONDS(3,1)=0;
                    isONOFFDS(3,1)=0;
                    %%%%%%%%%%%%%%
                    
                    if sum(isONDS)>sum(isONOFFDS)
                        %if sum(isONDS)>=2 && sum(isONOFFDS)<2
                        mergedmeta.meanRepDB.ROIdata{1,ROIn}.summary.ONorONOFForOFF=1;
                        mergedmeta.ONorONOFForOFF(ROIn,1)=1;
                    elseif sum(isONDS)==sum(isONOFFDS) && offPeak==0
                        mergedmeta.meanRepDB.ROIdata{1,ROIn}.summary.ONorONOFForOFF=1;
                        mergedmeta.ONorONOFForOFF(ROIn,1)=1;
                    elseif sum(isONDS)==sum(isONOFFDS) && offPeak==1
                        mergedmeta.meanRepDB.ROIdata{1,ROIn}.summary.ONorONOFForOFF=2;
                        mergedmeta.ONorONOFForOFF(ROIn,1)=2;
                    elseif sum(isONOFFDS)>sum(isONDS) && offPeak==1 && slopeOff<-0.01
                        %elseif sum(isONOFFDS)>=2 && sum(isONDS)<2 && offPeak==1
                        mergedmeta.meanRepDB.ROIdata{1,ROIn}.summary.ONorONOFForOFF=2;
                        mergedmeta.ONorONOFForOFF(ROIn,1)=2;
                    elseif sum(isONOFFDS)>sum(isONDS) && offPeak==0
                        %elseif sum(isONOFFDS)>=2 && sum(isONDS)<2 && offPeak==0
                        mergedmeta.meanRepDB.ROIdata{1,ROIn}.summary.ONorONOFForOFF=4;  % ONorONOFForOFF=4 denotes ON transient response
                        mergedmeta.ONorONOFForOFF(ROIn,1)=4;
                    elseif sum(isONOFFDS)>sum(isONDS) && offPeak==1 && slopeOff>=-0.01
                        %elseif sum(isONOFFDS)>=2 && sum(isONDS)<2 && offPeak==0
                        mergedmeta.meanRepDB.ROIdata{1,ROIn}.summary.ONorONOFForOFF=5;  % ONorONOFForOFF=5 denotes ON sustained OFF response
                        mergedmeta.ONorONOFForOFF(ROIn,1)=5;
                        %                     elseif sum(isONOFFDS)== sum(isONDS)
                        %                     %elseif sum(isONOFFDS)<2 && sum(isONDS)<2
                        %                         mergedmeta.meanRepDB.ROIdata{1,ROIn}.summary.ONorONOFForOFF=0;
                        %                         mergedmeta.ONorONOFForOFF(ROIn,1)=0;
                        %                     elseif sum(isONOFFDS)>=2 && sum(isONDS)>=2
                        %                         mergedmeta.meanRepDB.ROIdata{1,ROIn}.summary.ONorONOFForOFF=0;
                        %                         mergedmeta.ONorONOFForOFF(ROIn,1)=0;
                    end
                    
                    % excludes from analysis every cell with a weak
                    %                    and noisy response (more than 4 detected peaks)
                    if size(peakLocMax,1)>5
                        mergedmeta.meanRepDB.ROIdata{1,ROIn}.summary.ONorONOFForOFF=0;
                        mergedmeta.ONorONOFForOFF(ROIn,1)=0;
                    end
                    
                    mergedmeta.meanRepDB.ROIdata{ROIn}.summary.isONDS=isONDS;
                    mergedmeta.meanRepDB.ROIdata{ROIn}.summary.isONOFFDS=isONOFFDS;
                    mergedmeta.meanRepDB.ROIdata{ROIn}.summary.maxminRatio=maxminRatio;
                    mergedmeta.meanRepDB.ROIdata{ROIn}.summary.peakLocMaxfSec=peakLocMaxfSec;
                    mergedmeta.meanRepDB.ROIdata{ROIn}.summary.slope=slope;
                    mergedmeta.meanRepDB.ROIdata{ROIn}.summary.sumResponse=sumResponse;
                    mergedmeta.isONDS(ROIn,1:4)=isONDS';
                    mergedmeta.isONOFFDS(ROIn,1:4)=isONOFFDS';
                    
                    %                 else
                    %                     mergedmeta.meanRepDB.ROIdata{ROIn}.summary.ONorONOFForOFF=0;
                    %                     mergedmeta.ONorONOFForOFF(ROIn,1)=0;
                    %end
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                
                
                
                
                
                
                
                
                dispfigh=figure(2000+ROIn);clf;
                set(dispfigh,'position',[965 150 950 800]);
                for p=1:size(r,2)
                    subplot(ceil((size(r,2)+1)/3),3,p);
                    for nrep=1:repDB
                        for nruns=1:mergedmeta.metadataDB{nrep}.log.numberofruns
                            plot(mergedmeta.metadataDB{1,nrep}.ROIdata{1,ROIn}.ordfluotime{1,nruns},mergedmeta.metadataDB{1,nrep}.ROIdata{1,ROIn}.nororddff{1,nruns}(:,p),'Color',[0.7,0.7,0.7],'linewidth',1);
                            hold on
                        end
                    end
                    plot(mergedmeta.meanRepDB.ROIdata{1,1}.summary.ordfluotime,r(:,p),'r','linewidth',2);
                    title(strcat(num2str(a1(p),'%6.0f'), ' degrees'));
                    maxy=max(max(r(:,1:size(r,2))));miny=min(min(r(:,1:size(r,2))));
                    for nrep=1:repDB
                        for nruns=1:mergedmeta.metadataDB{nrep}.log.numberofruns
                            maxym(nrep*nruns)=max(max(mergedmeta.metadataDB{1,nrep}.ROIdata{1,ROIn}.nororddff{1,nruns}(:,1:size(r,2))));
                            minym(nrep*nruns)=min(min(mergedmeta.metadataDB{1,nrep}.ROIdata{1,ROIn}.nororddff{1,nruns}(:,1:size(r,2))));
                        end
                    end
                    maxym=max(maxym);
                    minym=min(minym);
                    minF=min([miny minym]);
                    maxF=max([maxy maxym]);
                    ylim([minF maxF]);
                    xlim([0 max(mergedmeta.meanRepDB.ROIdata{1,1}.summary.ordfluotime)]);
                    
                    
                    
z=750;      % at non-right angles, if the bar start to move within the screnn (rather than outside of it, as is the case for right angles).                     
screenWidthum=3500;  % 3 mm at the plane of the retina
pixelprojwidth=screenWidthum/600;   % 3 mm / 600 pixels
radDendTreeum=300;
corrSpeedFactor=1.2;
delay=0;
                              
speedterm=pixelprojwidth*corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed;
t=(screenWidthum-widthum)/2;
b=roiLocumX;
c=radDendTreeum;
k=roiLocumY;
e=widthum;
f=heightum;
g=sqrt((f-k)^2+(b-e/2)^2);
beta=acosd((b-e/2)/g);
h=screenWidthum;
m=sqrt(h^2+h^2)/2;  % half diagonal of the screen



switch p*45-45
%switch mergedmeta.metadataDB{1,1}.log.tested(p)
    case 0
        onset1=(delay+t+b-c)/speedterm;
        onset2=(delay+t+b)/speedterm;
    case 45
        x=g*sind(beta-45);
        n=m-x;
        onset1=(delay+n-z-c)/speedterm;
        onset2=(delay+n-z)/speedterm;
    case 180
        onset1=(delay+t+(e-b)-c)/speedterm;
        onset2=(delay+t+(e-b))/speedterm;
    case 225
        x=g*sind(beta-45);
        n=m+x;
        onset1=(delay+n-z-c)/speedterm;
        onset2=(delay+n-z)/speedterm;
    case 90
        onset1=(delay+t+k-c)/speedterm;
        onset2=(delay+t+k)/speedterm;
    case 135
        x=g*cosd(beta-45);
        n=m-x;
        onset1=(delay+n-z-c)/speedterm;
        onset2=(delay+n-z)/speedterm;
    case 270
        onset1=(delay+t+(e-f)+(f-k)-c)/speedterm;
        onset2=(delay+t+(e-f)+(f-k))/speedterm;
    case 315
        x=g*cosd(beta-45);
        n=m+x;
        onset1=(delay+n-z-c)/speedterm;
        onset2=(delay+n-z)/speedterm;
end

%                     switch mergedmeta.metadataDB{1,1}.log.tested(p)
%                         case 0
%                             onset1=(((screenWidthum-widthum)/2)+roiLocumX-radDendTreeum)/(pixelprojwidth*corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed);
%                             onset2=(((screenWidthum-widthum)/2)+roiLocumX)/(pixelprojwidth*corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed);
%                         case 45
%                             onset1=((((screenWidthum-widthum)/2)+roiLocumX-radDendTreeum)/cosd(45))/(pixelprojwidth*corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed);
%                             onset2=((((screenWidthum-widthum)/2)+roiLocumX)/cosd(45))/(pixelprojwidth*corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed);
%                         case 180
%                             onset1=(((screenWidthum-widthum)/2)+widthum-roiLocumX-radDendTreeum)/(pixelprojwidth*corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed);
%                             onset2=(((screenWidthum-widthum)/2)+widthum-roiLocumX)/(pixelprojwidth*corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed);
%                         case 225
%                             onset1=((((screenWidthum-widthum)/2)+widthum-roiLocumX-radDendTreeum)/cosd(45))/(pixelprojwidth*corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed);
%                             onset2=((((screenWidthum-widthum)/2)+widthum-roiLocumX/cosd(45)))/(pixelprojwidth*corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed);
%                         case 90
%                             onset1=(((screenWidthum-widthum)/2)+roiLocumY-radDendTreeum)/(pixelprojwidth*corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed);
%                             onset2=(((screenWidthum-widthum)/2)+roiLocumY)/(pixelprojwidth*corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed);
%                         case 135
%                             onset1=((((screenWidthum-widthum)/2)+widthum-roiLocumX-radDendTreeum)/cosd(45))/(pixelprojwidth*corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed);
%                             onset2=((((screenWidthum-widthum)/2)+widthum-roiLocumX/cosd(45)))/(pixelprojwidth*corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed);
%                         case 270
%                             onset1=(((screenWidthum-widthum)/2)+widthum-heightum+heightum-roiLocumY-radDendTreeum)/(pixelprojwidth*corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed);
%                             onset2=(((screenWidthum-widthum)/2)+widthum-heightum+heightum-roiLocumY)/(pixelprojwidth*corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed);
%                         case 315
%                             onset1=((((screenWidthum-widthum)/2)+widthum-roiLocumX-radDendTreeum)/cosd(45))/(pixelprojwidth*corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed);
%                             onset2=((((screenWidthum-widthum)/2)+widthum-roiLocumX/cosd(45)))/(pixelprojwidth*corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed);
%                     end


                    offset1=onset1+(mergedmeta.metadataDB{1,1}.log.width/(corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed));
                    offset2=onset2+(mergedmeta.metadataDB{1,1}.log.width/(corrSpeedFactor*60*mergedmeta.metadataDB{1,1}.log.speed));
                    
                    hold on;
                    plot([onset1 onset1],[0 10],'b');
                    hold on;
                    plot([onset2 onset2],[0 10],'b');
                    hold on;
                    plot([offset1 offset1],[0 10],'b');
                    hold on;
                    plot([offset2 offset2],[0 10],'b');
                end
                
                
                
                subplot(ceil((size(r,2)+1)/3),3,p+1);
                dpm=[d,d(1)]; apm=[a,a(1)];
                maxap=0;
                for nrep=1:repDB
                    for nruns=1:mergedmeta.metadataDB{nrep}.log.numberofruns
                        ap=[mergedmeta.metadataDB{1,nrep}.ROIdata{1,ROIn}.meanr{1,nruns},mergedmeta.metadataDB{1,nrep}.ROIdata{1,ROIn}.meanr{1,nruns}(1)];
                        maxap=max([maxap apm ap]);
                    end
                end
                
                for nrep=1:repDB
                    for nruns=1:mergedmeta.metadataDB{nrep}.log.numberofruns
                        ap=[mergedmeta.metadataDB{1,nrep}.ROIdata{1,ROIn}.meanr{1,nruns},mergedmeta.metadataDB{1,nrep}.ROIdata{1,ROIn}.meanr{1,nruns}(1)];
                        P=polarSS2p(deg2rad(dpm), ap/max(ap),1,'b',1);
                        hold on
                    end
                end
                P=polarSS2p(deg2rad(dpm), apm/max(apm),1,'b',2);
                hold on;
                [X,Y]=pol2cart(deg2rad(pdir), 1);
                compassSS2p(X,Y,1,'r',2);
                
                
                for nrep=1:repDB
                    for nruns=1:mergedmeta.metadataDB{nrep}.log.numberofruns
                        hold on
                        [X,Y]=pol2cart(deg2rad(mergedmeta.metadataDB{1,nrep}.ROIdata{1,ROIn}.pdir{1,nruns}), 1);
                        compassSS2p(X,Y,1,'r',1);
                    end
                end
                
                %title(strcat('n = ',num2str(nr,'%6.0f'),',pDir = ',num2str(pdir,'%6.0f'),',SDpDir = ',num2str(sd,'%6.0f'), ', DSI = ',num2str(DSI,'%6.3f'), ', DI = ',num2str(DI,'%6.3f')));
                annotation('textbox', [0.55, 0.27, 0.1, 0.1], 'string', ['n = ',num2str(mergedmeta.meanRepDB.ROIdata{ROIn}.summary.n,'%6.0f'),', pDir = ',num2str(pdir,'%6.0f'),', SDpDir = ',num2str(sd,'%6.0f'),', pOri = ',num2str(pori,'%6.0f'),', oriR2 = ',num2str(oriR2,'%6.3f'),', SDpOri = ',num2str(SDpOri,'%6.0f')],'LineStyle','none')
                annotation('textbox', [0.6, 0, 0.1, 0.1], 'string', ['DSI = ',num2str(DSI,'%6.3f'), ', DI = ',num2str(DI,'%6.3f'), ', OSI = ',num2str(OSI,'%6.3f'), ', globalOSI = ',num2str(globalOSI,'%6.3f')],'LineStyle','none')
                annotation('textbox', [0.2, 0, 0.7, 0.05], 'string', ['isDS = ',num2str(mergedmeta.meanRepDB.ROIdata{ROIn}.summary.isDS,'%d'), ', isOS = ',num2str(mergedmeta.meanRepDB.ROIdata{ROIn}.summary.isOS,'%d'), ', ONorONOFForOFF = ',num2str(mergedmeta.meanRepDB.ROIdata{ROIn}.summary.ONorONOFForOFF,'%d')],'LineStyle','none','FontSize',12,'Color','r');
                hold off;
                drawnow nocallbacks
                
            end
        end
        
        
        try mergedmeta.isONDS=[mergedmeta.isONDS; zeros(size(mergedmeta.isDS,1)-size(mergedmeta.isONDS,1),4)]; end;
        try mergedmeta.isONOFFDS=[mergedmeta.isONOFFDS; zeros(size(mergedmeta.isDS,1)-size(mergedmeta.isONOFFDS,1),4)]; end;
        
        
        
        if strcmp('DS',type)
            for ROIn=1:size(mergedmeta.metadataDG{1}.ROIdata,2)
                tabledata.data(ROIn,1)=mergedmeta.isDS(ROIn,1);
                tabledata.data(ROIn,2)=mergedmeta.ONorONOFForOFF(ROIn,1);
                if  isfield(mergedmeta.meanRepDB.ROIdata{ROIn}.summary,'isCART')
                    tabledata.data(ROIn,3)=mergedmeta.isCART(ROIn,1);
                else
                    tabledata.data(ROIn,3)=0;
                end
                if  isfield(mergedmeta.meanRepDB.ROIdata{ROIn}.summary,'isRBPMS')
                    tabledata.data(ROIn,4)=mergedmeta.isRBPMS(ROIn,1);
                else
                    tabledata.data(ROIn,4)=0;
                end
                if  isfield(mergedmeta.meanRepDB.ROIdata{ROIn}.summary,'isRetro')
                    tabledata.data(ROIn,5)=mergedmeta.isRetro(ROIn,1);
                else
                    tabledata.data(ROIn,5)=0;
                end
                tabledata.label={'isDS','ONorONOFForOFF','isCART','isRBPMS','isRetro'};
            end
        elseif strcmp('OS',type)
             for ROIn=1:size(mergedmeta.metadataDG{1}.ROIdata,2)
                tabledata.data(ROIn,1)=mergedmeta.isOS(ROIn,1);
                tabledata.data(ROIn,2)=mergedmeta.ONorONOFForOFF(ROIn,1);
                if  isfield(mergedmeta.meanRepDB.ROIdata{ROIn}.summary,'isCART')
                    tabledata.data(ROIn,3)=mergedmeta.isCART(ROIn,1);
                else
                    tabledata.data(ROIn,3)=0;
                end
                if  isfield(mergedmeta.meanRepDB.ROIdata{ROIn}.summary,'isRBPMS')
                    tabledata.data(ROIn,4)=mergedmeta.isRBPMS(ROIn,1);
                else
                    tabledata.data(ROIn,4)=0;
                end
                if  isfield(mergedmeta.meanRepDB.ROIdata{ROIn}.summary,'isRetro')
                    tabledata.data(ROIn,5)=mergedmeta.isRetro(ROIn,1);
                else
                    tabledata.data(ROIn,5)=0;
                end
                tabledata.label={'isOS','ONorONOFForOFF','isCART','isRBPMS','isRetro'};
            end
        end
        save('tabledata.mat', 'tabledata');
        
        f=cellTypeSelectionGui;
        waitfor(f);
        
        load tabledata
        
        if strcmp('DS',type)
            for ROIn=1:size(mergedmeta.metadataDG{1}.ROIdata,2)
                mergedmeta.meanRepDG.ROIdata{ROIn}.summary.isDS=tabledata.data(ROIn,1);
                mergedmeta.isDS(ROIn,1)=tabledata.data(ROIn,1);
                mergedmeta.meanRepDB.ROIdata{ROIn}.summary.ONorONOFForOFF=tabledata.data(ROIn,2);
                mergedmeta.ONorONOFForOFF(ROIn,1)=tabledata.data(ROIn,2);
                mergedmeta.meanRepDB.ROIdata{ROIn}.summary.isCART=tabledata.data(ROIn,3);
                mergedmeta.isCART(ROIn,1)=tabledata.data(ROIn,3);
                mergedmeta.meanRepDB.ROIdata{ROIn}.summary.isRBPMS=tabledata.data(ROIn,4);
                mergedmeta.isRBPMS(ROIn,1)=tabledata.data(ROIn,4);
                mergedmeta.meanRepDB.ROIdata{ROIn}.summary.isRetro=tabledata.data(ROIn,5);
                mergedmeta.isRetro(ROIn,1)=tabledata.data(ROIn,5);
            end
        elseif strcmp('OS',type)
            for ROIn=1:size(mergedmeta.metadataDG{1}.ROIdata,2)
                mergedmeta.meanRepDG.ROIdata{ROIn}.summary.isOS=tabledata.data(ROIn,1);
                mergedmeta.isOS(ROIn,1)=tabledata.data(ROIn,1);
                mergedmeta.meanRepDB.ROIdata{ROIn}.summary.ONorONOFForOFF=tabledata.data(ROIn,2);
                mergedmeta.ONorONOFForOFF(ROIn)=tabledata.data(ROIn,2);
                mergedmeta.meanRepDB.ROIdata{ROIn}.summary.isCART=tabledata.data(ROIn,3);
                mergedmeta.isCART(ROIn,1)=tabledata.data(ROIn,3);
                mergedmeta.meanRepDB.ROIdata{ROIn}.summary.isRBPMS=tabledata.data(ROIn,4);
                mergedmeta.isRBPMS(ROIn,1)=tabledata.data(ROIn,4);
                mergedmeta.meanRepDB.ROIdata{ROIn}.summary.isRetro=tabledata.data(ROIn,5);
                mergedmeta.isRetro(ROIn,1)=tabledata.data(ROIn,5);
            end
        end
        delete tabledata.mat
    end
    
    
    savdir=cd;
    save(([mergedmetaFileList{filen}]), 'mergedmeta');
    
%     clear all; clear java; 
javaaddpath(which('MatlabGarbageCollector.jar'))
jheapcl;
    
    
end
end


