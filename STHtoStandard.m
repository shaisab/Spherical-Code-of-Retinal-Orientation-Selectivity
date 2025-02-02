function [StU,StV,StVecU,StVecV] = STHtoStandard(S,TH,VecS,VecTH,discPoints)

%% Setting Parameters
[n,~]=size(S);

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

StU=[];
StV=[];
StRHO=[];
StRHOS=[];
StRHOTH=[];
StF=[];
StFS=[];
StFTH=[];

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
        StV(i)=interp3(StS3,StTH3,StV3,S(i),TH(i));
        
        StRHO(i)=interp3(StS3,StTH3,StRHO3,S(i),TH(i));
        StRHOS(i)=interp3(StS3,StTH3,StRHOS3,S(i),TH(i));
        StRHOTH(i)=interp3(StS3,StTH3,StRHOTH3,S(i),TH(i));
        
        StF(i)=interp3(StS3,StTH3,StF3,S(i),TH(i));
        StFS(i)=interp3(StS3,StTH3,StFS3,S(i),TH(i));
        StFTH(i)=interp3(StS3,StTH3,StFTH3,S(i),TH(i));
        
    else
        %Calculation of closest grid point.
        StU(i)=interp2(StS4,StTH4,StU4,S(i),TH(i));
        StV(i)=interp4(StS4,StTH4,StV4,S(i),TH(i));
        
        StRHO(i)=interp4(StS4,StTH4,StRHO4,S(i),TH(i));
        StRHOS(i)=interp4(StS4,StTH4,StRHOS4,S(i),TH(i));
        StRHOTH(i)=interp4(StS4,StTH4,StRHOTH4,S(i),TH(i));
        
        StF(i)=interp4(StS4,StTH4,StF4,S(i),TH(i));
        StFS(i)=interp4(StS4,StTH4,StFS4,S(i),TH(i));
        StFTH(i)=interp4(StS4,StTH4,StFTH4,S(i),TH(i));
    end
    
    StVecU(i)=VecS(i)*(StRHOS(i).*cos(StF(i))-StRHO(i).*StFS(i).*sin(StF(i)))+VecTH(i).*(StRHOTH(i).*cos(StF(i))-StRHO(i).*StFTH(i).*sin(StF(i)));
    StVecV(i)=VecS(i)*(StRHOS(i).*sin(StF(i))+StRHO(i).*StFS(i).*cos(StF(i)))+VecTH(i).*(StRHOTH(i).*sin(StF(i))+StRHO(i).*StFTH(i).*cos(StF(i)));
    
end


