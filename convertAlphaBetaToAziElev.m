
[allPossibleFieldsName,PathName] = uigetfile('*.mat','Select a allPossibleFields translation or rotation in reitinal space');
load(allPossibleFieldsName);

filedir=cd;
names=getFileList(filedir,'retinal to global lookup table.mat',0,'anywhere');
load(num2str(names{1}));
elevR=real(elevR);
elevL=real(elevL);

resolution=5;
%resolution=allPossibleFields.metadata.resolution;
alphaold=0:resolution:360;
betaold=0:resolution:180;
[Aold,Bold]=meshgrid(alphaold,betaold);

f1=figure;
contourf(alphaold,betaold,allPossibleFields.meanOutput);    % plot original alpha-beta map

% interpolate alpha-beta map and put it is OutputNew
newRes=0.5;
alphanew=0:newRes:360;
betanew=0:newRes:180;
[Anew,Bnew]=meshgrid(alphanew,betanew);
OutputNew=griddata(Aold,Bold,allPossibleFields.meanOutput,Anew,Bnew);
clear allPossibleFields

GSR=zeros(size(OutputNew,1)*size(OutputNew,2),3);
GSL=zeros(size(OutputNew,1)*size(OutputNew,2),3);
k=1;
for i=1:size(OutputNew,1)
    for j=1:size(OutputNew,2)
GSR(k,:)=[aziR(i,j),elevR(i,j),OutputNew(i,j)];
GSL(k,:)=[aziL(i,j),elevL(i,j),OutputNew(i,j)];
k=k+1;
    end
end

%% plot azimuth-elevation map of Right eye before and after interpolation
f2=figure; 
scatter(GSR(:,1),GSR(:,2),500,GSR(:,3),'filled');
maxValR=max(max(GSR(:,3)));
maxAziR=GSR(GSR(:,3)==maxValR,1);
maxElevR=GSR(GSR(:,3)==maxValR,2);

outputRes=1;
alpha=0:outputRes:360;
beta=0:outputRes:180;
[A,B]=meshgrid(alpha,beta);
CR=griddata(real(GSR(:,1)),real(GSR(:,2)),real(GSR(:,3)),A,B);

%  for k = 1:size(GSR,1)
% x(k,1) = isreal(GSR(k,3));
%  end
%  find(x==0)
 
f3=figure;
contourf(A,B,CR);

% convert axis and martixes to new global space conventions
elements=floor(size(CR,2)/2);
CRdOutput=circshift(CR,[0 elements]);   %shift azimuth 0 to middle of graph
CRdOutput(:,elements+1)=[];
CRdOutput=[CRdOutput(:,end),CRdOutput];
CRdOutput=fliplr(CRdOutput);
f4=figure;
azi=-180:outputRes:180;
elev=-90:outputRes:90;
ax=contourf(azi,elev,CRdOutput); 
hold on
plot([-180 180],[0 0],'m','LineWidth',1);
plot([0 0],[-90 90],'m','LineWidth',1);
xlim([-180 180]);ylim([-90 90]);
grid on;
ax=gca;
ax.XTick = -180:45:180;
ax.YTick = -90:45:90;
ax.XTickLabel = {'left','-135','-90','-45','ahead','45','90','135','right'};
ax.FontSize=12;
ax.YTickLabel = {'nadir','-45', 'horizon','45','zenith'};

%% store data of Right eye is a structure
allPossibleFields.meanOutput=CRdOutput;
allPossibleFields.maxAziR=maxAziR;
allPossibleFields.maxElevR=maxElevR;

%% save allPossibleFields file for Right eye
filename=allPossibleFieldsName;
ind=strfind(filename,'_');
filenameNew=filename(1:ind(10));
currenttime=datestr(now,'mmm_dd_yyyy_HH_MM');
save([filenameNew,currenttime,'_global space','_Right','.mat'],'allPossibleFields');
savefig(f4,[filenameNew,currenttime,'_global space','_Right']);
clear allPossibleFields

%% plot azimuth-elevation map of Left eye before and after interpolation
f5=figure; 
scatter(GSL(:,1),GSL(:,2),500,GSL(:,3),'filled');
maxValL=max(max(GSL(:,3)));
maxAziL=GSL(GSL(:,3)==maxValL,1);
maxElevL=GSL(GSL(:,3)==maxValL,2);

CL=griddata(real(GSL(:,1)),real(GSL(:,2)),real(GSL(:,3)),A,B);
f6=figure;
contourf(A,B,CL);

% convert axis and martixes to new global space conventions
elements=floor(size(CL,2)/2);
CLdOutput=circshift(CL,[0 elements]);   %shift azimuth 0 to middle of graph
CLdOutput(:,elements+1)=[];
CLdOutput=[CLdOutput(:,end),CLdOutput];
CLdOutput=fliplr(CLdOutput);
f7=figure;
ax=contourf(azi,elev,CLdOutput); 
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

%% store data of Left eye is a structure
allPossibleFields.meanOutput=CLdOutput;
allPossibleFields.maxAziL=maxAziL;
allPossibleFields.maxElevL=maxElevL;

%% save allPossibleFields file for Left eye
filename=allPossibleFieldsName;
ind=strfind(filename,'_');
filenameNew=filename(1:ind(10));
currenttime=datestr(now,'mmm_dd_yyyy_HH_MM');
save([filenameNew,currenttime,'_global space','_Left','.mat'],'allPossibleFields');
savefig(f7,[filenameNew,currenttime,'_global space','_Left']);


