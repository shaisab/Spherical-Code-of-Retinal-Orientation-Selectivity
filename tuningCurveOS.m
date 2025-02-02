function []= tuningCurveOS()

% figure;

filedir=cd;
%% Load files
[FileName,PathName] = uigetfile('*.mat','Select pooledmap file for all clusters');

load(FileName);


tuning.dir=[0,45,90,135,180,225,270,315];
tuning.OSI=pooledmap.OSI;
tuning.symmetry_ratio=pooledmap.symmetry_ratio;
tuning.poriSz=pooledmap.poriSz;
tuning.meanr=pooledmap.meanr;
tuning.meanrs=pooledmap.meanrs;
try tuning.typeClust=pooledmap.typeClust; end

%ON SUS - 3, ON TRANS - 8, ON-OFF SUS - 4, ON-OFF TRANS - 6
            
for j=1:size(tuning.OSI,1)
    
%     [~,p]=max(tuning.meanr(j,:));
%     tuning.poriInd(j,1)=p;
%     shift=4-tuning.poriInd(j,1);       % shift the max response to be at the middle of the graph
%     a=circshift(tuning.meanr(j,:),shift,2);
    
    a=pooledmap.meanr(j,:);
%     ar=pooledmap.meanr(j,:);
    d=tuning.dir;
%     dr=pooledmap.theta(j,:);
%     pori=tuning.poriSz(j);
    [oriFWHM,u,k,gof2,x,xq,Rq]=doubleGaussianFitOSI(j,d,a);
%     if oriFWHM<45
%     plotPolarAndQuiverPooledMapOS(ar,dr,pori);
%     end
    tuning.poriFit(j)=u;
    tuning.oriR2Fit(j)=gof2.rsquare;
    tuning.oriFWHM(j)=oriFWHM;
    tuning.k(j)=k;
    tuning.gof{j}=gof2;
    tuning.x(:,j)=x;
    tuning.xq(:,j)=xq';
%     tuning.R(:,j)=R;
%     tuning.Rq(:,j)=Rq;
    tuning.Rq(:,j)=Rq/max(Rq);
end
       
   
    
%     shift x and R to the middle of the graph and plot them
%     figure;
%     for j=1:size(tuning.R,2)
%         [~,p]=max(tuning.R(:,j));
%         tuning.poriInd(j,1)=p;
%         shift=floor(size(tuning.R,1)/2)+1-tuning.poriInd(j,1);       % shift the max response to be at the middle of the graph
%         tuning.RShifted(:,j)=circshift(tuning.R(:,j),shift,1);
%         plot(tuning.x,tuning.RShifted(:,j),'bo-');
%         hold on
%     end
%     
%     
%     shift xq and Rq to the middle of the graph and plot them
%     figure;
    for j=1:size(tuning.Rq,2)
        [pks,locs] = findpeaks(tuning.Rq(:,j));
        tuning.poriIndq(j,1)=min(locs);
%         shift=floor(size(tuning.Rq,1)/2)+1-tuning.poriIndq(j,1);       % shift the max response to be at the middle of the graph

        shift(j)=size(tuning.Rq,1)+1-tuning.poriIndq(j,1);       % shift the max response to be at the beginning of the graph

        tuning.RqShifted(:,j)=circshift(tuning.Rq(:,j),shift(j),1);
%         plot(xq,tuning.RqShifted(:,j),'b');
%         plot(xq,tuning.Rq(:,j),'b');
%         hold on
    end
%     
%     
%     
%     calculate the mean interpolated tuning curve
    tuning.meanRqShiftedAll=mean(tuning.RqShifted,2);
%     tuning.xshifted=-180:1:135;
    tuning.xshifted=0:1:315;
    figure;
    plot(tuning.xshifted',tuning.meanRqShiftedAll,'k');

    hold on
    
[Max_meanRqShiftedAll,Ind_max] = max(tuning.meanRqShiftedAll);

    figure;
    tuning.meanRq=mean(tuning.Rq,2);
    plot(tuning.xq(:,1),tuning.meanRq,'k');
    
    hold on

%     figure;
%     plot(tuning.xshifted(181:end)',tuning.meanRqShiftedAll(181:end));
%     plot(tuning.xshifted(Ind_max:Ind_max+90)',tuning.meanRqShiftedAll(Ind_max:Ind_max+90));
%     xtmp=(0:1:360);
%     shadedErrorBar(0:round(size(tuning.RqShifted,1)/2)-1,mean(tuning.RqShifted(181:end,:),2),std(tuning.RqShifted(181:end,:),[],2),'k');
%     shadedErrorBar(0:round(size(tuning.RqShifted,1)/2)-2,mean(tuning.RqShifted(1:180,:),2),std(tuning.RqShifted(1:180,:),[],2),'k');
%     shadedErrorBar(0:90,mean(tuning.RqShifted(Ind_max:Ind_max+90,:),2),std(tuning.RqShifted(Ind_max:Ind_max+90,:),[],2),'k');
%     xlim([0 89]); ylim([0 1]);
%     xlabel('Angle from preferred direction (deg.)');
%     ylabel('Response');

    figure;
    plot(tuning.xshifted(Ind_max:Ind_max+90)',tuning.meanRqShiftedAll(Ind_max:Ind_max+90));
    shadedErrorBar(0:90,mean(tuning.RqShifted(Ind_max:Ind_max+90,:),2),std(tuning.RqShifted(Ind_max:Ind_max+90,:),[],2),'k');
    xlim([0 89]); ylim([0 1]);
    xlabel('Angle from preferred direction (deg.)');
    ylabel('Response');

%     figure;
%     plot(tuning.xq(1:91,1)',tuning.meanRq(1:91));
%     shadedErrorBar(0:90,mean(tuning.Rq(1:91,:),2),std(tuning.Rq(1:91,:),[],2),'k');
%     xlim([0 89]); ylim([0 1]);
%     xlabel('Angle from preferred direction (deg.)');
%     ylabel('Response');
    
    
    %ON SUS - 3, ON TRANS - 8, ON-OFF SUS - 4, ON-OFF TRANS - 6
    
%     tuning.oriFWHM_ON_SUS=tuning.oriFWHM(tuning.typeClust==3);
%     tuning.oriFWHM_ON_TRANS=tuning.oriFWHM(tuning.typeClust==8);
%     tuning.oriFWHM_ONOFF_SUS=tuning.oriFWHM(tuning.typeClust==4);
%     tuning.oriFWHM_ONOFF_TRANS=tuning.oriFWHM(tuning.typeClust==6);
%     
%     tuning.oriFWHM_ON_SUSmean=mean(tuning.oriFWHM_ON_SUS);
%     tuning.oriFWHM_ON_TRANSmean=mean(tuning.oriFWHM_ON_TRANS);
%     tuning.oriFWHM_ONOFF_SUSmean=mean(tuning.oriFWHM_ONOFF_SUS);
%     tuning.oriFWHM_ONOFF_TRANSmean=mean(tuning.oriFWHM_ONOFF_TRANS);
%     
%     tuning.oriFWHM_ON_SUSstd=std(tuning.oriFWHM_ON_SUS);
%     tuning.oriFWHM_ON_TRANSstd=std(tuning.oriFWHM_ON_TRANS);
%     tuning.oriFWHM_ONOFF_SUSstd=std(tuning.oriFWHM_ONOFF_SUS);
%     tuning.oriFWHM_ONOFF_TRANSstd=std(tuning.oriFWHM_ONOFF_TRANS);
%     
%     [p,tbl,stats] = anova1(tuning.oriFWHM,tuning.typeClust)
%     [c,m,h,gnames] = multcompare(stats)
       
%     tuning.meanx=0:round(size(tuning.RqShifted,1)/2)-1;
%     tuning.meanR=mean(tuning.RqShifted(181:end,:),2);
    
    tuning.meanx=0:90;
    tuning.meanR=mean(tuning.RqShifted(Ind_max:Ind_max+90,:),2);
    
%     % finds the gaussian that describes the mean tuning curve
%     [oriFWHM,u,k,gof2,x,R,~,~]= doubleGaussianFitOSI(1,xtmp,tuning.meanRq');
%     tuning.meanporiFit=u;
%     tuning.meanporiR2Fit=gof2.rsquare;
%     tuning.meanoriFWHM=oriFWHM;
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
