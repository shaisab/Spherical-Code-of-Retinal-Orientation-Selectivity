function [StU,StV,StVecU,StVecV] = XYZtoStandard(X,Y,Z,VecX,VecY,VecZ,discPoints)

%% Setting Parameters
[n,~]=size(X);


%% Data
phys=dlmread('phys');

R=phys(1);
M=phys(2); %This should be normalized to 1.

m1=phys(3);
m2=phys(4);
m3=phys(5);
m4=phys(6);

nu=phys(7);

a1=phys(8);
a2=phys(9);
a3=phys(10);
a4=phys(11);

ntemp=discPoints;

%% Converting to STH
S=R*acos((R-Z)./R);
TH=mod(atan2(Y,X),2*pi);

VecS=R./(Z.*(2*R-Z)).^(1/2).*VecZ;
VecTH=VecX.*(-Y./(X.^2+Y.^2))+VecY.*(X./(X.^2+Y.^2));

%% Sector 1 data
StSec1Data=dlmread(['StandardSector1_',num2str(discPoints)]);

StRHO=StSec1Data(1:ntemp,1:ntemp);
StF=StSec1Data(1:ntemp,ntemp+1:2*ntemp);
StS=StSec1Data(1:ntemp,2*ntemp+1:3*ntemp);
StTH=StSec1Data(1:ntemp,3*ntemp+1:4*ntemp);

%Calcution of derivatives
ds=StS(1,2)-StS(1,1);
DS=DiffX(ntemp)/ds;
StRHOS=(DS*StRHO')';
StFS=(DS*StF')';

dth=StTH(2,1)-StTH(1,1);
DTH=DiffX(ntemp)/dth;
StRHOTH=(DTH*StRHO);
StFTH=(DTH*StF);

StU=StRHO.*cos(StF);
StV=StRHO.*sin(StF);

StU1=StU;
StV1=StV;
StRHO1=StRHO;
StRHOS1=StRHOS;
StRHOTH1=StRHOTH;
StF1=StF;
StFS1=StFS;
StFTH1=StFTH;
StS1=StS;
StTH1=StTH;


%% Sector 2 data
StSec2Data=dlmread(['StandardSector2_',num2str(discPoints)]);

StRHO=StSec2Data(1:ntemp,1:ntemp);
StF=StSec2Data(1:ntemp,ntemp+1:2*ntemp);
StS=StSec2Data(1:ntemp,2*ntemp+1:3*ntemp);
StTH=StSec2Data(1:ntemp,3*ntemp+1:4*ntemp);

%Calcution of derivatives
ds=StS(1,2)-StS(1,1);
DS=DiffX(ntemp)/ds;
StRHOS=(DS*StRHO')';
StFS=(DS*StF')';

dth=StTH(2,1)-StTH(1,1);
DTH=DiffX(ntemp)/dth;
StRHOTH=(DTH*StRHO);
StFTH=(DTH*StF);

StU=StRHO.*cos(StF);
StV=StRHO.*sin(StF);

StU2=StU;
StV2=StV;
StRHO2=StRHO;
StRHOS2=StRHOS;
StRHOTH2=StRHOTH;
StF2=StF;
StFS2=StFS;
StFTH2=StFTH;
StS2=StS;
StTH2=StTH;

%% Sector 3 data
StSec3Data=dlmread(['StandardSector3_',num2str(discPoints)]);

StRHO=StSec3Data(1:ntemp,1:ntemp);
StF=StSec3Data(1:ntemp,ntemp+1:2*ntemp);
StS=StSec3Data(1:ntemp,2*ntemp+1:3*ntemp);
StTH=StSec3Data(1:ntemp,3*ntemp+1:4*ntemp);

%Calcution of derivatives
ds=StS(1,2)-StS(1,1);
DS=DiffX(ntemp)/ds;
StRHOS=(DS*StRHO')';
StFS=(DS*StF')';

dth=StTH(2,1)-StTH(1,1);
DTH=DiffX(ntemp)/dth;
StRHOTH=(DTH*StRHO);
StFTH=(DTH*StF);

StU=StRHO.*cos(StF);
StV=StRHO.*sin(StF);

StU3=StU;
StV3=StV;
StRHO3=StRHO;
StRHOS3=StRHOS;
StRHOTH3=StRHOTH;
StF3=StF;
StFS3=StFS;
StFTH3=StFTH;
StS3=StS;
StTH3=StTH;

%% Sector 4 data
StSec4Data=dlmread(['StandardSector4_',num2str(discPoints)]);

StRHO=StSec4Data(1:ntemp,1:ntemp);
StF=StSec4Data(1:ntemp,ntemp+1:2*ntemp);
StS=StSec4Data(1:ntemp,2*ntemp+1:3*ntemp);
StTH=StSec4Data(1:ntemp,3*ntemp+1:4*ntemp);

%Calcution of derivatives
ds=StS(1,2)-StS(1,1);
DS=DiffX(ntemp)/ds;
StRHOS=(DS*StRHO')';
StFS=(DS*StF')';

dth=StTH(2,1)-StTH(1,1);
DTH=DiffX(ntemp)/dth;
StRHOTH=(DTH*StRHO);
StFTH=(DTH*StF);

StU=StRHO.*cos(StF);
StV=StRHO.*sin(StF);


StU4=StU;
StV4=StV;
StRHO4=StRHO;
StRHOS4=StRHOS;
StRHOTH4=StRHOTH;
StF4=StF;
StFS4=StFS;
StFTH4=StFTH;
StS4=StS;
StTH4=StTH;

StU=zeros(n,1);
StV=zeros(n,1);
StRHO=zeros(n,1);
StRHOS=zeros(n,1);
StRHOTH=zeros(n,1);
StF=zeros(n,1);
StFS=zeros(n,1);
StFTH=zeros(n,1);

for i=1:n
    
    if TH(i)>=0 & TH(i)<pi/2,
        %Calculation of closest grid point.
        StU(i)=interp2(StS1,StTH1,StU1,S(i),TH(i));
        StV(i)=interp2(StS1,StTH1,StV1,S(i),TH(i));
        
        StRHO(i)=interp2(StS1,StTH1,StRHO1,S(i),TH(i));
        StRHOS(i)=interp2(StS1,StTH1,StRHOS1,S(i),TH(i));
        StRHOTH(i)=interp2(StS1,StTH1,StRHOTH1,S(i),TH(i));
        
        StF(i)=interp2(StS1,StTH1,StF1,S(i),TH(i));
        StFS(i)=interp2(StS1,StTH1,StFS1,S(i),TH(i));
        StFTH(i)=interp2(StS1,StTH1,StFTH1,S(i),TH(i));
        
    elseif TH(i)>=pi/2 & TH(i)<pi,
        %Calculation of closest grid point.
        StU(i)=interp2(StS2,StTH2,StU2,S(i),TH(i));
        StV(i)=interp2(StS2,StTH2,StV2,S(i),TH(i));
        
        StRHO(i)=interp2(StS2,StTH2,StRHO2,S(i),TH(i));
        StRHOS(i)=interp2(StS2,StTH2,StRHOS2,S(i),TH(i));
        StRHOTH(i)=interp2(StS2,StTH2,StRHOTH2,S(i),TH(i));
        
        StF(i)=interp2(StS2,StTH2,StF2,S(i),TH(i));
        StFS(i)=interp2(StS2,StTH2,StFS2,S(i),TH(i));
        StFTH(i)=interp2(StS2,StTH2,StFTH2,S(i),TH(i));
        
    elseif TH(i)>=pi & TH(i)<3*pi/2,
        %Calculation of closest grid point.
        StU(i)=interp2(StS3,StTH3,StU3,S(i),TH(i));
        StV(i)=interp2(StS3,StTH3,StV3,S(i),TH(i));
        
        StRHO(i)=interp2(StS3,StTH3,StRHO3,S(i),TH(i));
        StRHOS(i)=interp2(StS3,StTH3,StRHOS3,S(i),TH(i));
        StRHOTH(i)=interp2(StS3,StTH3,StRHOTH3,S(i),TH(i));
        
        StF(i)=interp2(StS3,StTH3,StF3,S(i),TH(i));
        StFS(i)=interp2(StS3,StTH3,StFS3,S(i),TH(i));
        StFTH(i)=interp2(StS3,StTH3,StFTH3,S(i),TH(i));
        
    else
        %Calculation of closest grid point.
        StU(i)=interp2(StS4,StTH4,StU4,S(i),TH(i));
        StV(i)=interp2(StS4,StTH4,StV4,S(i),TH(i));
        
        StRHO(i)=interp2(StS4,StTH4,StRHO4,S(i),TH(i));
        StRHOS(i)=interp2(StS4,StTH4,StRHOS4,S(i),TH(i));
        StRHOTH(i)=interp2(StS4,StTH4,StRHOTH4,S(i),TH(i));
        
        StF(i)=interp2(StS4,StTH4,StF4,S(i),TH(i));
        StFS(i)=interp2(StS4,StTH4,StFS4,S(i),TH(i));
        StFTH(i)=interp2(StS4,StTH4,StFTH4,S(i),TH(i));
    end
    
    StVecU(i)=VecS(i)*(StRHOS(i).*cos(StF(i))-StRHO(i).*StFS(i).*sin(StF(i)))+VecTH(i).*(StRHOTH(i).*cos(StF(i))-StRHO(i).*StFTH(i).*sin(StF(i)));
    StVecV(i)=VecS(i)*(StRHOS(i).*sin(StF(i))+StRHO(i).*StFS(i).*cos(StF(i)))+VecTH(i).*(StRHOTH(i).*sin(StF(i))+StRHO(i).*StFTH(i).*cos(StF(i)));
    
end


%% plot standard retina surface %%%%%%%%
f1=figure;
StSec1Data=dlmread(['StandardSector1_',num2str(discPoints)]);
StSec2Data=dlmread(['StandardSector2_',num2str(discPoints)]);
StSec3Data=dlmread(['StandardSector3_',num2str(discPoints)]);
StSec4Data=dlmread(['StandardSector4_',num2str(discPoints)]);

n=size(StSec1Data,1);

StRHO1=StSec1Data(1:n,1:n);
StF1=StSec1Data(1:n,n+1:2*n);
StU1=StRHO1.*cos(StF1);
StV1=StRHO1.*sin(StF1);

StRHO2=StSec2Data(1:n,1:n);
StF2=StSec2Data(1:n,n+1:2*n);
StU2=StRHO2.*cos(StF2);
StV2=StRHO2.*sin(StF2);

StRHO3=StSec3Data(1:n,1:n);
StF3=StSec3Data(1:n,n+1:2*n);
StU3=StRHO3.*cos(StF3);
StV3=StRHO3.*sin(StF3);

StRHO4=StSec4Data(1:n,1:n);
StF4=StSec4Data(1:n,n+1:2*n);
StU4=StRHO4.*cos(StF4);
StV4=StRHO4.*sin(StF4);

hold on;
s1=surf(StU1,StV1,zeros(size(StU1,1),size(StU1,1)),'FaceColor','k','EdgeColor', 'none');
s2=surf(StU2,StV2,zeros(size(StU2,1),size(StU2,1)),'FaceColor','k','EdgeColor', 'none');
s3=surf(StU3,StV3,zeros(size(StU3,1),size(StU3,1)),'FaceColor','k','EdgeColor', 'none');
s4=surf(StU4,StV4,zeros(size(StU4,1),size(StU4,1)), 'FaceColor','k','EdgeColor', 'none');
alpha(.15);
uistack(s1,'bottom');uistack(s2,'bottom');uistack(s3,'bottom');uistack(s4,'bottom');
hold on;


