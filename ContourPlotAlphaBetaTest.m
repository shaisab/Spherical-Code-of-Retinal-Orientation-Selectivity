function []=ContourPlotAlphaBeta(resolution)

VizAngle=deg2rad(180);
discpoints=40;

phys=dlmread('phys');

R=phys(1);
M=phys(2);

% VizAngle should be in radians


        [FileName,PathName] = uigetfile('*.mat','Select a pooledMapAlphaCorr file');

load(FileName); 
%load('allPossibleFields_Dec_03_2014_40_DS_created_Apr_08_2015_22_05.mat')

c=10;
cutoff=(['c',num2str(c)]);
k=1;
h=1;
% for Alpha=0:5:355
% for beta=90:5:180
for Alpha=0:resolution:360
for beta=90:resolution:180
    mat180(k,h)=allPossibleFields.matchind180(k,h).(cutoff);
    mat90(k,h)=allPossibleFields.matchind90(k,h).(cutoff);
    k=k+1;
end
k=1;
h=h+1;
end


data=mat90;

%% Standard Contour
[a,b]=size(data);

Alpha1=0.6257+pi;
beta1=(pi/2-1.2624)+pi/2;
Alpha2=0.8848;
beta2=3.0135;
Alpha3=-1.0012+pi;
beta3=(pi/2-1.14502)+pi/2;

Adeg=Alpha;
Bdeg=beta;

[Adeg,Bdeg]=meshgrid(Adeg,Bdeg);
    
Alpha=linspace(0,2*pi,b+1);
Alpha=Alpha(2:b+1);

beta=linspace(pi-VizAngle/2,pi,a);


[A,B]=meshgrid(Alpha,beta);

figure;
hold on;
contourf(A,B,data,'EdgeColor','none','LineStyle','none');
plot(Alpha1,beta1,'ok',Alpha2,beta2,'ok',Alpha3,beta3,'ok');
hold off;
shading interp;

%% Contour on Sphere

X=R*cos(A).*sin(B);
Y=R*sin(A).*sin(B);
Z=-R*cos(B);

X(:,b+1)=X(:,1);
Y(:,b+1)=Y(:,1);
Z(:,b+1)=Z(:,1);
data(:,b+1)=data(:,1);

% x=X(1:10:size(X,1),1:10:size(X,2));
% y=Y(1:10:size(Y,1),1:10:size(Y,2));
% z=X(1:10:size(X,1),1:10:size(Z,2));

figure;
% mesh(x,y,z);
hold on;
surf(X,Y,Z,data,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
shading interp;
alpha(0.95);
axis equal

thetafine=linspace(0,2*pi,40);
sfine=linspace(pi/2,pi,10);

[THfine,Sfine]=meshgrid(thetafine,sfine);

Xfine=(R+0.005)*cos(THfine).*sin(Sfine);
Yfine=(R+0.005)*sin(THfine).*sin(Sfine);
Zfine=-(R+0.005)*cos(Sfine);

surf(Xfine,Yfine,Zfine,'FaceColor','none');

hold off;
%% Contour on Standard retina

%Mapping to S-TH coordinates
Sdata=-2*M/VizAngle*acos(-Z/R)+2*M*pi/VizAngle;
THdata=atan2(Y,X)+pi;


StSec1Data=dlmread(['StandardSector1_',num2str(discpoints)]);
StSec2Data=dlmread(['StandardSector2_',num2str(discpoints)]);
StSec3Data=dlmread(['StandardSector3_',num2str(discpoints)]);
StSec4Data=dlmread(['StandardSector4_',num2str(discpoints)]);

n=size(StSec1Data,1);


StS1=StSec1Data(1:n,2*n+1:3*n);
StTH1=StSec1Data(1:n,3*n+1:4*n);
StRHO1=StSec1Data(1:n,1:n);
StF1=StSec1Data(1:n,n+1:2*n);

StU1=StRHO1.*cos(StF1);
StV1=StRHO1.*sin(StF1);


StS2=StSec2Data(1:n,2*n+1:3*n);
StTH2=StSec2Data(1:n,3*n+1:4*n);
StRHO2=StSec2Data(1:n,1:n);
StF2=StSec2Data(1:n,n+1:2*n);

StU2=StRHO2.*cos(StF2);
StV2=StRHO2.*sin(StF2);


StS3=StSec3Data(1:n,2*n+1:3*n);
StTH3=StSec3Data(1:n,3*n+1:4*n);
StRHO3=StSec3Data(1:n,1:n);
StF3=StSec3Data(1:n,n+1:2*n);

StU3=StRHO3.*cos(StF3);
StV3=StRHO3.*sin(StF3);


StS4=StSec4Data(1:n,2*n+1:3*n);
StTH4=StSec4Data(1:n,3*n+1:4*n);
StRHO4=StSec4Data(1:n,1:n);
StF4=StSec4Data(1:n,n+1:2*n);

StU4=StRHO4.*cos(StF4);
StV4=StRHO4.*sin(StF4);



for i=1:a,
    for j=1:b,
        if Adeg(i,j)>=0 &&  Adeg(i,j)<90,
            StUdata(i,j)=griddata(StS1,StTH1,StU1,Sdata(i,j),THdata(i,j));
            StVdata(i,j)=griddata(StS1,StTH1,StV1,Sdata(i,j),THdata(i,j));
        elseif Adeg(i,j)>=90 &&  Adeg(i,j)<180,
          StUdata(i,j)=griddata(StS2,StTH2,StU2,Sdata(i,j),THdata(i,j));
          StVdata(i,j)=griddata(StS2,StTH2,StV2,Sdata(i,j),THdata(i,j));
        elseif Adeg(i,j)>=180 &&  Adeg(i,j)<270,
          StUdata(i,j)=griddata(StS3,StTH3,StU3,Sdata(i,j),THdata(i,j));
          StVdata(i,j)=griddata(StS3,StTH3,StV3,Sdata(i,j),THdata(i,j)); 
        else
            StUdata(i,j)=griddata(StS4,StTH4,StU4,Sdata(i,j),THdata(i,j));
            StVdata(i,j)=griddata(StS4,StTH4,StV4,Sdata(i,j),THdata(i,j));
        end
    end
end

StUdata(isnan(StUdata))=0;
StVdata(isnan(StVdata))=0;



StUdata1=StUdata;
StVdata1=StVdata;
data1=data;

StUdata1(~(THdata>=0 &  THdata<pi/2))=NaN;
StVdata1(~(THdata>=0 &  THdata<pi/2))=NaN;
data1(~(THdata>=0 &  THdata<pi/2))=NaN;

StUdata2=StUdata;
StVdata2=StVdata;
data2=data;

StUdata2(~(THdata>=pi/2 &  THdata<pi))=NaN;
StVdata2(~(THdata>=pi/2 &  THdata<pi))=NaN;
data2(~(THdata>=pi/2 &  THdata<pi))=NaN;

StUdata3=StUdata;
StVdata3=StVdata;
data3=data;

StUdata3(~(THdata>=pi &  THdata<3*pi/2))=NaN;
StVdata3(~(THdata>=pi &  THdata<3*pi/2))=NaN;
data3(~(THdata>=pi &  THdata<3*pi/2))=NaN;

StUdata4=StUdata;
StVdata4=StVdata;
data4=data;

StUdata4(~(THdata>=3*pi/2 &  THdata<2*pi))=NaN;
StVdata4(~(THdata>=3*pi/2 &  THdata<2*pi))=NaN;
data4(~(THdata>=3*pi/2 &  THdata<2*pi))=NaN;

figure;
hold on;
contourf(StUdata1,StVdata1,data1,'EdgeColor','none','LineStyle','none');
contourf(StUdata2,StVdata2,data2,'EdgeColor','none','LineStyle','none');
contourf(StUdata3,StVdata3,data3,'EdgeColor','none','LineStyle','none');
contourf(StUdata4,StVdata4,data4,'EdgeColor','none','LineStyle','none');
hold off;

shading interp;
