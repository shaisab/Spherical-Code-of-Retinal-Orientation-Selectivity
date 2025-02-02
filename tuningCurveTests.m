clear all
load('tuning curves_OS data_All_cells');

tuning.meanRqShifted_On_Sus=mean(tuning.RqShifted(:,tuning.typeClust==3),2);
tuning.meanRqShifted_On_Trans=mean(tuning.RqShifted(:,tuning.typeClust==8),2);
tuning.meanRqShifted_On_Off_Sus=mean(tuning.RqShifted(:,tuning.typeClust==4),2);
tuning.meanRqShifted_On_Off_Trans=mean(tuning.RqShifted(:,tuning.typeClust==6),2);
%     tuning.xshifted=-180:1:135;
    tuning.xshifted=0:1:315;
    
tuning.RqShifted_On_Sus=tuning.RqShifted(:,tuning.typeClust==3);
tuning.RqShifted_On_Trans=tuning.RqShifted(:,tuning.typeClust==8);
tuning.RqShifted_On_Off_Sus=tuning.RqShifted(:,tuning.typeClust==4);
tuning.RqShifted_On_Off_Trans=tuning.RqShifted(:,tuning.typeClust==6);

[Max_meanRqShiftedAll,Ind_max] = max(tuning.meanRqShiftedAll);

    figure;
    
    shadedErrorBar(0:90,mean(tuning.RqShifted_On_Sus(Ind_max:Ind_max+90,:),2),std(tuning.RqShifted_On_Trans(Ind_max:Ind_max+90,:),[],2),'-r',1); 
    
    hold on
    
    shadedErrorBar(0:90,mean(tuning.RqShifted_On_Trans(Ind_max:Ind_max+90,:),2),std(tuning.RqShifted_On_Trans(Ind_max:Ind_max+90,:),[],2),'-g',1);
    
    hold on
    
    shadedErrorBar(0:90,mean(tuning.RqShifted_On_Off_Sus(Ind_max:Ind_max+90,:),2),std(tuning.RqShifted_On_Off_Sus(Ind_max:Ind_max+90,:),[],2),'-b',1);
    
    hold on
    
    shadedErrorBar(0:90,mean(tuning.RqShifted_On_Off_Trans(Ind_max:Ind_max+90,:),2),std(tuning.RqShifted_On_Off_Trans(Ind_max:Ind_max+90,:),[],2),'-k',1);
    xlim([0 89]); ylim([0 1]);
    xlabel('Angle from preferred direction (deg.)');
    ylabel('Response');
    
    hold on
    
%     plot(tuning.xshifted(Ind_max:Ind_max+90)',tuning.meanRqShifted_On_Sus(Ind_max:Ind_max+90));
%     plot(tuning.xshifted(Ind_max:Ind_max+90)',tuning.meanRqShifted_On_Trans(Ind_max:Ind_max+90),'r');
%     plot(tuning.xshifted(Ind_max:Ind_max+90)',tuning.meanRqShifted_On_Off_Sus(Ind_max:Ind_max+90),'g');
%     plot(tuning.xshifted(Ind_max:Ind_max+90)',tuning.meanRqShifted_On_Off_Trans(Ind_max:Ind_max+90),'b');
    
legend('ON Sus','ON Sus','ON Trans','ON Trans','ON OFF Sus','ON OFF Sus','ON OFF Trans','ON OFF Trans')

    
    
    