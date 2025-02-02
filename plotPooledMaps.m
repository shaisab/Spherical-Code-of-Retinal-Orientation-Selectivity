



% current directory should be the 'All maps' directory


namesMap=getFileList(filedir,'map_',0,'anywhere');
for i=1:size(namesMap,2)   
    if strfind(namesMap{i},'.mat')
    Allmaps(i)=load(num2str(namesMap{i}));
    end
end

SzMapOutputDS=[];
for i=1:size(Allmaps,2)
SzMapOutputDS=[SzMapOutputDS;Allmaps(1,i).map.SzMapOutputDS]; 
end

ONDS=[];
for i=1:size(Allmaps,2)
ONDS=logical([ONDS;Allmaps(1,i).map.ONDS]); 
end

ONOFFDS=[];
for i=1:size(Allmaps,2)
ONOFFDS=logical([ONOFFDS;Allmaps(1,i).map.ONOFFDS]); 
end

retroONDS=[];
for i=1:size(Allmaps,2)
    if isfield(Allmaps(1,i).map,'isRetro')
        retroONDS=logical([retroONDS;Allmaps(1,i).map.retroONDS]);
    else
        retroONDS=logical([retroONDS;zeros(size(Allmaps(1,i).map.ONDS,1),1)]);
    end
end

retroONOFFDS=[];
for i=1:size(Allmaps,2)
    if isfield(Allmaps(1,i).map,'isRetro')
        retroONOFFDS=logical([retroONOFFDS;Allmaps(1,i).map.retroONOFFDS]);
    else
        retroONOFFDS=logical([retroONOFFDS;zeros(size(Allmaps(1,i).map.ONDS,1),1)]);
    end
end

retroNonDS=[];
for i=1:size(Allmaps,2)
    if isfield(Allmaps(1,i).map,'isRetro')
        retroNonDS=logical([retroNonDS;Allmaps(1,i).map.retroOther]);
    else
        retroNonDS=logical([retroNonDS;zeros(size(Allmaps(1,i).map.ONDS,1),1)]);
    end
end

% retroONtDS=[];
% for i=1:size(Allmaps,2)
%     if isfield(Allmaps(1,i).map,'isRetro')
%         retroONtDS=logical([retroONtDS;Allmaps(1,i).map.ONtDS]);
%     else
%         retroONtDS=logical([retroONtDS;zeros(size(Allmaps(1,i).map.ONDS,1),1)]);
%     end
% end

f1=figure;
set(f1,'position',[450 150 1500 500]);
subplot(1,4,1)
compassSS2p(SzMapOutputDS(ONDS,3),SzMapOutputDS(ONDS,4),1,'r',1);
annotation('textbox', [0.18, 0.85, 0.1, 0.1], 'string','ON DS','LineStyle','none','FontSize',14);
subplot(1,4,2)
compassSS2p(SzMapOutputDS(ONOFFDS,3),SzMapOutputDS(ONOFFDS,4),1,'g',1);
annotation('textbox', [0.37, 0.85, 0.1, 0.1], 'string','ON-OFF DS','LineStyle','none','FontSize',14);
subplot(1,4,3)
compassSS2p(SzMapOutputDS(retroONDS,3),SzMapOutputDS(retroONDS,4),1,'b',1);
annotation('textbox', [0.57, 0.85, 0.1, 0.1], 'string','retroON DS','LineStyle','none','FontSize',14);
subplot(1,4,4)
compassSS2p(SzMapOutputDS(retroONOFFDS,3),SzMapOutputDS(retroONOFFDS,4),1,'b',1);
annotation('textbox', [0.77, 0.85, 0.1, 0.1], 'string','retroONOFF DS','LineStyle','none','FontSize',14);
hold off


        percentONDS=100*(sum(ONDS)/length(ONDS));
        percentONOFFDS=100*(sum(ONOFFDS)/length(ONOFFDS));
        percentRetroONDS=100*(sum(retroONDS)/length(retroONDS));
        percentRetroONOFFDS=100*(sum(retroONOFFDS)/length(retroONOFFDS));
        percentRetroNonDS=100*(sum(retroNonDS)/length(retroNonDS));
       
        map.ONDS=ONDS;
        map.ONOFFDS=ONOFFDS;
        map.retroONDS=retroONDS;
        map.retroONOFFDS=retroONOFFDS;
        map.retroNonDS=retroNonDS;
       
        map.numONDS=sum(ONDS);
        map.numONOFFDS=sum(ONOFFDS);
        map.numRetroONDS=sum(retroONDS);
        map.numRetroONOFFDS=sum(retroONOFFDS);
        map.numRetroNonDS=sum(retroNonDS);
        
        map.SzMapOutputDS=SzMapOutputDS;

saveas(f1,['polarPlot_allMaps_'],'fig');
save(['map_All','.mat'],'map');




