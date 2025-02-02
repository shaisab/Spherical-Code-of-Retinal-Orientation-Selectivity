
function [] = plotintersectionLine(discPoints)

%% calculate equation of saccular plane (values taken from worksheet 'all mice saculae' in the file 'LandmarkData-Saccular-Macula-28-9-15' 
%SMs-SMi (SM1)
%normal=[-0.144,0.315,0.938];        % vector of the sacular macula
normal=[-0.119,-0.286,0.951];        % vector of the sacular macula Symmetric
[alpharot,betarot]=AB(normal);
alpha1=rad2deg(alpharot);
beta1=rad2deg(betarot);

%SMf-SMi (SM2)
%normal=[0.220,-0.680,0.699];        % vector of the sacular macula
normal=[0.246,-0.695,0.675];        % vector of the sacular macula Symmetric
[alpharot,betarot]=AB(normal);
alpha3=rad2deg(alpharot);
beta3=rad2deg(betarot);

%SMs-SMf (SM3)
%normal=[-0.391,-0.030,0.920];       % vector of the sacular macula
normal=[-0.381,0.121,0.917];       % vector of the sacular macula Symmetric
[alpharot,betarot]=AB(normal);
alpha2=rad2deg(alpharot);
beta2=rad2deg(betarot);



alphaAvg=deg2rad(mean([alpha1,alpha2,alpha3]));      % this is the average alpha for the three vectors of the saccular macula.
alphaAvgUtr=alphaAvg+deg2rad(93);                         % this is the average alpha based on Honda data. 
%It should be plus 93 because the normal vector of the utricle is farther away from the saccule vector relative to the nasal direction (rectus muscle)

%% calculate all sectors of standard retina
%Standard Sector 1
StSec1Data=dlmread(['StandardSector1_',num2str(discPoints)]);

StRHO1=StSec1Data(1:discPoints,1:discPoints);
StF1=StSec1Data(1:discPoints,discPoints+1:2*discPoints);
StS1=StSec1Data(1:discPoints,2*discPoints+1:3*discPoints);
StTH1=StSec1Data(1:discPoints,3*discPoints+1:4*discPoints);

StU1=StRHO1.*cos(StF1);
StV1=StRHO1.*sin(StF1);

%Calcution of derivatives
ds=StS1(1,2)-StS1(1,1);
DS=DiffX(discPoints)/ds;
StRHOS1=(DS*StRHO1')';
StFS1=(DS*StF1')';

dth=StTH1(2,1)-StTH1(1,1);
DTH=DiffX(discPoints)/dth;
StRHOTH1=(DTH*StRHO1);
StFTH1=(DTH*StF1);

%Standard Sector 2
StSec2Data=dlmread(['StandardSector2_',num2str(discPoints)]);

StRHO2=StSec2Data(1:discPoints,1:discPoints);
StF2=StSec2Data(1:discPoints,discPoints+1:2*discPoints);
StS2=StSec2Data(1:discPoints,2*discPoints+1:3*discPoints);
StTH2=StSec2Data(1:discPoints,3*discPoints+1:4*discPoints);

StU2=StRHO2.*cos(StF2);
StV2=StRHO2.*sin(StF2);

%Calcution of derivatives
ds=StS2(1,2)-StS2(1,1);
DS=DiffX(discPoints)/ds;
StRHOS2=(DS*StRHO2')';
StFS2=(DS*StF2')';

dth=StTH2(2,1)-StTH2(1,1);
DTH=DiffX(discPoints)/dth;
StRHOTH2=(DTH*StRHO2);
StFTH2=(DTH*StF2);


%Standard Sector 3
StSec3Data=dlmread(['StandardSector3_',num2str(discPoints)]);

StRHO3=StSec3Data(1:discPoints,1:discPoints);
StF3=StSec3Data(1:discPoints,discPoints+1:2*discPoints);
StS3=StSec3Data(1:discPoints,2*discPoints+1:3*discPoints);
StTH3=StSec3Data(1:discPoints,3*discPoints+1:4*discPoints);

StU3=StRHO3.*cos(StF3);
StV3=StRHO3.*sin(StF3);

%Calcution of derivatives
ds=StS3(1,2)-StS3(1,1);
DS=DiffX(discPoints)/ds;
StRHOS3=(DS*StRHO3')';
StFS3=(DS*StF3')';

dth=StTH3(2,1)-StTH3(1,1);
DTH=DiffX(discPoints)/dth;
StRHOTH3=(DTH*StRHO3);
StFTH3=(DTH*StF3);

%Standard Sector 4
StSec4Data=dlmread(['StandardSector4_',num2str(discPoints)]);

StRHO4=StSec4Data(1:discPoints,1:discPoints);
StF4=StSec4Data(1:discPoints,discPoints+1:2*discPoints);
StS4=StSec4Data(1:discPoints,2*discPoints+1:3*discPoints);
StTH4=StSec4Data(1:discPoints,3*discPoints+1:4*discPoints);

StU4=StRHO4.*cos(StF4);
StV4=StRHO4.*sin(StF4);

%Calcution of derivatives
ds=StS4(1,2)-StS4(1,1);
DS=DiffX(discPoints)/ds;
StRHOS4=(DS*StRHO4')';
StFS4=(DS*StF4')';

dth=StTH4(2,1)-StTH4(1,1);
DTH=DiffX(discPoints)/dth;
StRHOTH4=(DTH*StRHO4);
StFTH4=(DTH*StF4);

%% calculate the intersection line of the saccular macula plane and the standard retina

sret=linspace(0,1,100);            
for i=1:100,
        %Calculation of closest grid point sector 2.
        Uavg1(i)=interp2(1*StS2,StTH2,StU2,sret(i),alphaAvg);
        Vavg1(i)=interp2(1*StS2,StTH2,StV2,sret(i),alphaAvg); 
        %Calculation of closest grid point sector 4.
        Uavg2(i)=interp2(1*StS4,StTH4,StU4,sret(i),alphaAvg+pi);
        Vavg2(i)=interp2(1*StS4,StTH4,StV4,sret(i),alphaAvg+pi);              
end
hold on;
plot(Uavg1,Vavg1,'r-',Uavg2,Vavg2,'r-','LineWidth',3)

%% calculate the intersection line of the utricular macula plane and the standard retina

sret=linspace(0,1,100);            
for i=1:100,
        %Calculation of closest grid point sector 1.
        UavgUtr1(i)=interp2(1*StS1,StTH1,StU1,sret(i),mod(alphaAvgUtr+pi,2*pi));
        VavgUtr1(i)=interp2(1*StS1,StTH1,StV1,sret(i),mod(alphaAvgUtr+pi,2*pi)); 
        %Calculation of closest grid point sector 3.
        UavgUtr2(i)=interp2(1*StS3,StTH3,StU3,sret(i),alphaAvgUtr);
        VavgUtr2(i)=interp2(1*StS3,StTH3,StV3,sret(i),alphaAvgUtr);              
end
hold on;
plot(UavgUtr1,VavgUtr1,'b-',UavgUtr2,VavgUtr2,'b-','LineWidth',3)
end