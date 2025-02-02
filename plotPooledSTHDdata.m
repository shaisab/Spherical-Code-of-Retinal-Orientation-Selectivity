function []= plotPooledSTHDdata(typeStr)

% The current directory should be 'All maps'

[filedir,go]=bbdirselector('select data folder',cd);
if ~go disp('no folder selected'); return; end

namesSTH=getFileList(filedir,'STHData_',0,'anywhere');
for i=1:size(namesSTH,2)  
    if strfind(namesSTH{i},typeStr)
    AllSTH{i}=load(namesSTH{i});
    end
end

STHdata=[];
for i=1:size(AllSTH,2)
STHdata=[STHdata;AllSTH{i}]; 
end

Sn=STHdata(:,1);
THn=STHdata(:,2);
vecSn=STHdata(:,3);
vecTHn=STHdata(:,4);

f1=figure;
M=1;

switch typeStr
    case 'ONDS' 
quiver(rad2deg(Sn),rad2deg(THn),vecSn,vecTHn,.3,'r','linewidth',1);     
    case 'ONOFFDS' 
quiver(rad2deg(Sn),rad2deg(THn),vecSn,vecTHn,.3,'g','linewidth',1);          
    case 'OFFDS' 
quiver(rad2deg(Sn),rad2deg(THn),vecSn,vecTHn,.3,'b','linewidth',1);          
end
%axis([0,2*pi,0,2*pi]);
axis([0,360,0,360]);


f2=figure;
P=rad2deg(atan2(STHdata(:,4),STHdata(:,3)));
P1=[P(P>0); 360+P(P<0)];
pdir=rad2deg(STHdata(:,2))-P1;
pdir1=[pdir(pdir>0); 360+pdir(pdir<0)];
[X,Y]=pol2cart(deg2rad(pdir1), 1);
switch typeStr
    case 'ONDS'
        compassSS2p(X,Y,1,'r',3);
    case 'ONOFFDS'
        compassSS2p(X,Y,1,'g',3);
    case 'OFFDS'
        compassSS2p(X,Y,1,'b',3);
end

dlmwrite(['STHDataPooled_',typeStr],[Sn,THn,vecSn,vecTHn]);
saveas(f1,['STHDataPooled_',typeStr],'fig');
saveas(f2,['STHtoPolarDataPooled_',typeStr],'fig');
hold off;

end
