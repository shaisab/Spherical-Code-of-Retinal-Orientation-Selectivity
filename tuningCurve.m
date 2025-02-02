
function []= tuningCurve()

figure;

filedir=cd;
[mergedmetaFileList]=getFileList(filedir,'mergedmeta',0,'anywhere');
if size(mergedmetaFileList,2)>0
    j=1;
    for fov=1:size(mergedmetaFileList,2)
        load(mergedmetaFileList{fov});
        
        for ROIn=1:size(mergedmeta.metadataDG{1}.ROIdata,2)
            
            tuning.pdir(j,1)=mergedmeta.meanRepDG.ROIdata{1,ROIn}.summary.pdir;
            tuning.DSI(j,1)=mergedmeta.meanRepDG.ROIdata{1,ROIn}.summary.DSI;
            tuning.isDS(j,1)=mergedmeta.isDS(ROIn);
            tuning.ONorONOFForOFF(j,1)=mergedmeta.ONorONOFForOFF(ROIn);
            try
                tuning.isCART(j,1)=mergedmeta.isCART(ROIn);
                tuning.isRBPMS(j,1)=mergedmeta.isRBPMS(ROIn);
                tuning.isRetro(j,1)=mergedmeta.isRetro(ROIn);
            end
            tuning.DSI(j,1)=mergedmeta.meanRepDG.ROIdata{1,ROIn}.summary.DSI;
            tuning.meanr(j,:)=mergedmeta.meanRepDG.ROIdata{1,ROIn}.summary.meanr;
            
            %%% calculates the preferred direction 'pdir' based on fit to a Von Mises function
            [~,p]=max(mergedmeta.meanRepDG.ROIdata{1,ROIn}.summary.meanr);
            tuning.pdirInd(j,1)=p;
            shift=4-tuning.pdirInd(j,1);       % shift the max response to be at the middle of the graph
            a=circshift(tuning.meanr(j,:),shift,2);
            d=mergedmeta.meanRepDG.ROIdata{1,ROIn}.summary.dir;
            
            
            
            
            if tuning.isDS(j)
                [dirFWHM,u,k,gof2,x,R,xq,Rq]= singleGaussianFit(ROIn,d,a);
                tuning.pdirFit(j)=u;
                tuning.dirR2Fit(j)=gof2.rsquare;
                tuning.dirFWHM(j)=dirFWHM;
                tuning.k(j)=k;
                tuning.gof{j}=gof2;
                tuning.x(:,j)=x;
                tuning.xq(:,j)=xq';
                tuning.R(:,j)=R;
                tuning.Rq(:,j)=Rq;
                
                j=j+1;
            end
        end
    end
    
    % shift x and R to the middle of the graph and plot them
    figure;
    for j=1:size(tuning.R,2)
        [~,p]=max(tuning.R(:,j));
        tuning.pdirInd(j,1)=p;
        shift=floor(size(tuning.R,1)/2)+1-tuning.pdirInd(j,1);       % shift the max response to be at the middle of the graph
        tuning.RShifted(:,j)=circshift(tuning.R(:,j),shift,1);
        plot(tuning.x,tuning.RShifted(:,j),'bo-');
        hold on
    end
    
    
    % shift xq and Rq to the middle of the graph and plot them
    figure;
    for j=1:size(tuning.Rq,2)
        [~,p]=max(tuning.Rq(:,j));
        tuning.pdirIndq(j,1)=p;
        shift=floor(size(tuning.Rq,1)/2)+1-tuning.pdirIndq(j,1);       % shift the max response to be at the middle of the graph
        tuning.RqShifted(:,j)=circshift(tuning.Rq(:,j),shift,1);
        plot(xq,tuning.RqShifted(:,j),'b');
        hold on
    end
    
    
    
    % calculate the mean interpolated tuning curve
    tuning.meanRqShiftedAll=mean(tuning.RqShifted,2);
    tuning.meanRqShiftedON=mean(tuning.RqShifted(:,tuning.ONorONOFForOFF==1),2);        % for ON DS
    tuning.meanRqShiftedONOFF=mean(tuning.RqShifted(:,tuning.ONorONOFForOFF==2),2);     % for ON-OFF DS
    tuning.xshifted=-180:1:180;
    figure;
    plot(tuning.xshifted',tuning.meanRqShiftedAll,'k');
    hold on
    plot(tuning.xshifted',tuning.meanRqShiftedON,'b');
    plot(tuning.xshifted',tuning.meanRqShiftedONOFF,'r');
    figure;
    plot(tuning.xshifted(181:end)',tuning.meanRqShiftedAll(181:end));
    %xtmp=(0:1:360);
    shadedErrorBar(0:round(size(tuning.RqShifted,1)/2)-1,mean(tuning.RqShifted(181:end,:),2),std(tuning.RqShifted(181:end,:),[],2),'k');
    xlim([0 180]); ylim([0 1]);
    xlabel('Angle from preferred direction (deg.)');
    ylabel('Response');
    
    
    tuning.dirFWHM_ON=tuning.dirFWHM(tuning.ONorONOFForOFF==1);
    tuning.dirFWHM_ONOFF=tuning.dirFWHM(tuning.ONorONOFForOFF==2);
    [tuning.ttest.h,tuning.ttest.p,tuning.ttest.ci,tuning.ttest.stats]=ttest2(tuning.dirFWHM_ON,tuning.dirFWHM_ONOFF,'Vartype','unequal')
    tuning.dirFWHM_ONmean=mean(tuning.dirFWHM_ON);
    tuning.dirFWHM_ONOFFmean=mean(tuning.dirFWHM_ONOFF);
    tuning.dirFWHM_ONstd=std(tuning.dirFWHM_ON);
    tuning.dirFWHM_ONOFFstd=std(tuning.dirFWHM_ONOFF);    

    tuning.meanx=0:round(size(tuning.RqShifted,1)/2)-1;
    tuning.meanR=mean(tuning.RqShifted(181:end,:),2);
    
%     % finds the gaussian that describes the mean tuning curve
%     [dirFWHM,u,k,gof2,x,R,~,~]= singleGaussianFit(1,xtmp,tuning.meanRqShifted');
%     tuning.meanpdirFit=u;
%     tuning.meandirR2Fit=gof2.rsquare;
%     tuning.meandirFWHM=dirFWHM;
%     tuning.meank=rad2deg(k);
%     tuning.meangof=gof2;
    
%     
%     figure;
%     shadedErrorBar(0:round(size(tuning.RqShifted,1)/2)-1,R(181:361),std(tuning.RqShifted(181:end,:),[],2),'k');
%     xlim([0 180]); ylim([0 1]);
%     xlabel('Angle from preferred direction (deg.)');
%     ylabel('Response');
 

    ind=strfind(filedir,'\');
    filename=filedir;
    filename(ind(end):length(filename))=[];
    filename(1:ind(end-1))=[];
    save(['tuning curves_',filename,'.mat'],'tuning');
    
end


