
%% invert polarity of Left global maps probed with rotation
[allPossibleFieldsName,PathName] = uigetfile('*.mat','Select a allPossibleFields Left global map probed with rotation');
load(allPossibleFieldsName);

outputRes=1;
figure;
azi=-180:outputRes:180;
elev=-90:outputRes:90;
contourf(azi,elev,allPossibleFields.meanOutput)

tmp=allPossibleFields.meanOutput;
tmp=flipud(tmp);                    % flip matrix up side down
elements=floor(size(tmp,2)/2);
tmp=circshift(tmp,[0 elements]);    % shift azimuth by 180 degrees

allPossibleFields.meanOutput=tmp;   % assign the inverted matrix to meanOutput 

f1=figure;
ax=contourf(azi,elev,tmp);
hold on
plot([-180 180],[0 0],'m','LineWidth',1);
plot([0 0],[-90 90],'m','LineWidth',1);
xlim([-180 180]);ylim([-90 90]);
grid on
ax=gca;
ax.XTick = -180:45:180;
ax.YTick = -90:45:90;
ax.XTickLabel = {'left','-135','-90','-45','ahead','45','90','135','right'};
ax.FontSize=12;
ax.YTickLabel = {'nadir','-45', 'horizon','45','zenith'};

%% save inverted file and figure
filename=allPossibleFieldsName;
ind=strfind(filename,'_');
filenameNew=filename(1:ind(13));
singname=filename(ind(13):end-4);
currenttime=datestr(now,'mmm_dd_yyyy_HH_MM');
save([filenameNew,'inverse polarity',singname,'.mat'],'allPossibleFields');
savefig(f1,[filenameNew,'inverse polarity',singname]);




