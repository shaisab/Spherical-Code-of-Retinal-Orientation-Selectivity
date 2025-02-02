function [Uout,Vout] = StandardRemove(Uin,Vin)

discPoints=40;
error=.025;

%% Standard Retina Data
StSec1Data=dlmread(['StandardSector1_',num2str(discPoints)]);
StSec2Data=dlmread(['StandardSector2_',num2str(discPoints)]);
StSec3Data=dlmread(['StandardSector3_',num2str(discPoints)]);
StSec4Data=dlmread(['StandardSector4_',num2str(discPoints)]);

n=size(StSec1Data,1);

StRHO1=StSec1Data(1:n,1:n);
StF1=StSec1Data(1:n,n+1:2*n);
StS1=StSec1Data(1:n,2*n+1:3*n);

StTH1=StSec1Data(1:n,3*n+1:4*n);

StU1=StRHO1.*cos(StF1);
StV1=StRHO1.*sin(StF1);

%Calcution of derivatives
ds=StS1(1,2)-StS1(1,1);
DS=DiffX(n)/ds;
StRHO1S=(DS*StRHO1')';
StF1S=(DS*StF1')';

dth=StTH1(2,1)-StTH1(1,1);
DTH=DiffX(n)/dth;
StRHO1TH=(DTH*StRHO1);
StF1TH=(DTH*StF1);

StRHO2=StSec2Data(1:n,1:n);
StF2=StSec2Data(1:n,n+1:2*n);
StS2=StSec2Data(1:n,2*n+1:3*n);
StTH2=StSec2Data(1:n,3*n+1:4*n);

StU2=StRHO2.*cos(StF2);
StV2=StRHO2.*sin(StF2);

%Calcution of derivatives
ds=StS2(1,2)-StS2(1,1);
DS=DiffX(n)/ds;
StRHO2S=(DS*StRHO2')';
StF2S=(DS*StF2')';

dth=StTH2(2,1)-StTH2(1,1);
DTH=DiffX(n)/dth;
StRHO2TH=(DTH*StRHO2);
StF2TH=(DTH*StF2);

StRHO3=StSec3Data(1:n,1:n);
StF3=StSec3Data(1:n,n+1:2*n);
StS3=StSec3Data(1:n,2*n+1:3*n);
StTH3=StSec3Data(1:n,3*n+1:4*n);

StU3=StRHO3.*cos(StF3);
StV3=StRHO3.*sin(StF3);

%Calcution of derivatives
ds=StS3(1,2)-StS3(1,1);
DS=DiffX(n)/ds;
StRHO3S=(DS*StRHO3')';
StF3S=(DS*StF3')';

dth=StTH3(2,1)-StTH3(1,1);
DTH=DiffX(n)/dth;
StRHO3TH=(DTH*StRHO3);
StF3TH=(DTH*StF3);

StRHO4=StSec4Data(1:n,1:n);
StF4=StSec4Data(1:n,n+1:2*n);
StS4=StSec4Data(1:n,2*n+1:3*n);
StTH4=StSec4Data(1:n,3*n+1:4*n);

StU4=StRHO4.*cos(StF4);
StV4=StRHO4.*sin(StF4);

%Calcution of derivatives
ds=StS4(1,2)-StS4(1,1);
DS=DiffX(n)/ds;
StRHO4S=(DS*StRHO4')';
StF4S=(DS*StF4')';

dth=StTH4(2,1)-StTH4(1,1);
DTH=DiffX(n)/dth;
StRHO4TH=(DTH*StRHO4);
StF4TH=(DTH*StF4);

if Uin >= 0 && Vin >= 0
    dist=min(min((Uin-StU1).^2+(Vin-StV1).^2).^.5);
    if dist>error
        Uout=NaN;
        Vout=NaN;
    else
        Uout=Uin;
        Vout=Vin;
    end
elseif Uin <= 0 && Vin >= 0
    dist=min(min((Uin-StU2).^2+(Vin-StV2).^2).^.5);
    if dist>error
        Uout=NaN;
        Vout=NaN;
    else
        Uout=Uin;
        Vout=Vin;
    end
elseif Uin <= 0 && Vin <= 0
    dist=min(min((Uin-StU3).^2+(Vin-StV3).^2).^.5);
    if dist>error
        Uout=NaN;
        Vout=NaN;
    else
        Uout=Uin;
        Vout=Vin;
    end
elseif Uin >= 0 && Vin <= 0
    dist=min(min((Uin-StU4).^2+(Vin-StV4).^2).^.5);
    if dist>error
        Uout=NaN;
        Vout=NaN;
    else
        Uout=Uin;
        Vout=Vin;
    end
end

% dist
% figure;
% hold on;
% surf(StU1,StV1,zeros(40,40));
% surf(StU2,StV2,zeros(40,40));
% surf(StU3,StV3,zeros(40,40));
% surf(StU4,StV4,zeros(40,40));
% plot(Uout,Vout,'ro');
