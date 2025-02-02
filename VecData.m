

function []=VecData(alpha, beta, polarity, discPoints, typeStr)

% typestr can be either 'ONDS', 'ONOFFDS', 'OFFDS', or 'retroONDS'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This script compares the vector fields obtained from the elastic model
%   with that from Shai's data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Physical Parameters.
% This data is imported from data exported from the mapping script. These
% values correspond to measured values on the flattened retina.

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Construction of Vector Fields
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% STH Vector Fields
n=discPoints;
s=linspace(0,M,n);
th=linspace(0,2*pi,n);

[S,TH]=meshgrid(s,th);
for i=1:n
    for j=1:n,
        [t(i,j),phi(i,j)]=InverseSTH(alpha,beta,R,S(i,j),TH(i,j));       
            vecS(i,j)=-cos(S(i,j)/R)*cos(beta)*sin(phi(i,j))*cos(TH(i,j)-alpha)+...
            cos(S(i,j)/R)*cos(phi(i,j))*(sin(TH(i,j)-alpha))-...
            sin(S(i,j)/R)*sin(phi(i,j))*sin(beta);
            vecS(i,j)=R*sin(t(i,j)/R)*vecS(i,j);
            vecTH(i,j)=(cos(beta)*sin(phi(i,j))*sin(TH(i,j)-alpha)+...
            cos(phi(i,j))*cos(alpha-TH(i,j)));
            vecTH(i,j)=R*sin(t(i,j)/R)*vecTH(i,j);
            vecTH(i,j)=1/(R*sin(S(i,j)/R))*vecTH(i,j);
            
    end
end

vecStemp=vecS./(vecS.^2+vecTH.^2).^(1/2);
vecTHtemp=vecTH./(vecS.^2+vecTH.^2).^(1/2);

vecS=vecStemp./(vecStemp.^2+vecTHtemp.^2);
vecTH=vecTHtemp./(vecStemp.^2+vecTHtemp.^2);

figure(2);

hold on;

quiver(S,TH/(2*pi),vecS,vecTH,1,'linewidth',1);
axis([0,1,0,1]);

hold off;

%% Plotting Vector Field in Generalized (S,TH) coordinates
figure(3);
hold on;

X=S.*cos(TH);
Y=S.*sin(TH);

vecX=cos(TH).*vecS-S.*sin(TH).*vecTH;
vecY=sin(TH).*vecS+S.*cos(TH).*vecTH;

a=vecX./(vecX.^2+vecY.^2).^(1/2);
b=vecY./(vecX.^2+vecY.^2).^(1/2);

vecX=a;
vecY=b;

quiver(X,Y,vecX,vecY);
plot(R*beta*cos(alpha),R*beta*sin(alpha),'o');
axis([-1,1,-1,1]);
axis square;

xlabel('x^{\prime}');
ylabel('y^{\prime}');

hold off;

Xprime=X;
Yprime=Y;

vecXprime=vecX;
vecYprime=vecY;

%% Plotting Vector Fields on the sphere
figure(4);
hold on;
X=R*sin(S/R).*cos(TH);
Y=R*sin(S/R).*sin(TH);
Z=R-R*cos(S/R);

vecX=vecS.*cos(S/R).*cos(TH)-R*vecTH.*sin(S/R).*sin(TH);
vecY=vecS.*cos(S/R).*sin(TH)+R*vecTH.*sin(S/R).*cos(TH);
vecZ=vecS.*sin(S/R);

% Normalization of vector field
veca=vecX./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);
vecb=vecY./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);
vecc=vecZ./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);

vecX=veca;
vecY=vecb;
vecZ=vecc;

ntemp=discPoints;

s=linspace(0,M,ntemp);
th=linspace(a1,a2,ntemp);
[S,TH]=meshgrid(s,th);

surf(R*sin(S/R).*cos(TH),R*sin(S/R).*sin(TH),R-R*cos(S/R),S,'EdgeColor','none','FaceAlpha','interp','AlphaDataMapping','scaled',...
    'AlphaData',S);
colormap('bone');
alphamap('default')
view(104,6);

quiver3(X,Y,Z,vecX,vecY,vecZ);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

axis([-R,R,-R,R,0,2*R]);
axis square;

plot3(R*cos(alpha)*sin(beta),R*sin(alpha)*sin(beta),R-R*cos(beta),'o');
xlabel('x');
ylabel('y');
zlabel('z');
hold off;

%% Sector 1
%This part of the code computes the vector fields in the sectors. We
%read in the numerical solution computed using the script Mapping.m to
%compute the value of rho. The vector fields are then found by computing
%the pushforward in the coordinates (u,v).

%Reading in data;
Sec1Data=dlmread('Sector1'); %radial stretching factor.
n=size(Sec1Data,1);
RHO1=Sec1Data(1:n,1:n);
F1=Sec1Data(1:n,n+1:2*n);

%Record the number of points that were discretized over.

% Construction of domain in s th coordinates
s=linspace(0,M,n+1);
s=linspace(s(2),M,n);
th=linspace(a1,a2,n);


[S1,TH1]=meshgrid(s,th);
%Calculation of the vector field of rotations in the coordinates s and th.
%The variables (t,phi) denote geodesic polar coordinate around the axis of
%rotation. 


t=zeros(n,n);
phi=zeros(n,n);
vecS1=zeros(n,n);
vecTH1=zeros(n,n);

%This loop constructs the components of the vector field of rotations
%in the the s, th coordinates. The function InverseSTH determines the
%coordinates of (t,phi) given coordiantes s and th and is necessary to
%calculate the components of the vector field. 
for i=1:n
    for j=1:n,
        [t(i,j),phi(i,j)]=InverseSTH(alpha,beta,R,S1(i,j),TH1(i,j)); 
        
            vecS1(i,j)=-cos(S1(i,j)/R)*cos(beta)*sin(phi(i,j))*cos(TH1(i,j)-alpha)+...
                cos(S1(i,j)/R)*cos(phi(i,j))*(sin(TH1(i,j)-alpha))-...
                sin(S1(i,j)/R)*sin(phi(i,j))*sin(beta);
            
            vecS1(i,j)=R*sin(t(i,j)/R)*vecS1(i,j);
            
            vecTH1(i,j)=(cos(beta)*sin(phi(i,j))*sin(TH1(i,j)-alpha)+...
                cos(phi(i,j))*cos(alpha-TH1(i,j)));
        
            vecTH1(i,j)=R*sin(t(i,j)/R)*vecTH1(i,j);
            vecTH1(i,j)=1/(R*sin(S1(i,j)/R))*vecTH1(i,j);
            
    end
end
vecTH1=vecTH1.*polarity;

%Calcution of derivatives
ds=s(2)-s(1);
DS=DiffX(n)/ds;
RHOS1=(DS*RHO1')';
FS1=(DS*F1')';

dth=th(2)-th(1);
DTH=DiffX(n)/dth;
RHOTH1=(DTH*RHO1);
FTH1=(DTH*F1);

%Calculating the vector fields from the formula for the push forward of the
%derivative.
vecU1=vecS1.*(RHOS1.*cos(F1)-RHO1.*FS1.*sin(F1))+vecTH1.*(RHOTH1.*cos(F1)-RHO1.*FTH1.*sin(F1));
vecV1=vecS1.*(RHOS1.*sin(F1)+RHO1.*FS1.*cos(F1))+vecTH1.*(RHOTH1.*sin(F1)+RHO1.*FTH1.*cos(F1));

%Plotting vector field in the sector. We plot the physical retina so we can
%see where the vectors are located on the physical retina.
figure(1);
hold on;

C=abs(R.^2*sin(S1./R).^2-S1.^2);
surf(RHO1.*cos(F1),RHO1.*sin(F1),zeros(n,n),C);
view(0,90);
caxis([0,1]);
colorbar;
colormap('bone');
shading interp;

U1=RHO1.*cos(F1);
V1=RHO1.*sin(F1);
quiver(U1,V1,vecU1,vecV1,'b');
axis square;
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

plot(beta*R*cos(alpha),beta*R*sin(alpha),'o');
xlabel('u');
ylabel('v');

hold off;

%% Sector 2
%This part of the code computes the vector fields in the sectors. We
%read in the numerical solution computed using the script Mapping.m to
%compute the value of rho. The vector fields are then found by computing
%the pushforward in the coordinates (u,v).

%Reading in data;
Sec2Data=dlmread('Sector2'); %radial stretching factor.
n=size(Sec2Data,1);
RHO2=Sec2Data(1:n,1:n);
F2=Sec2Data(1:n,n+1:2*n);

%Record the number of points that were discretized over.

% Construction of domain in s th coordinates
s=linspace(0,M,n+1);
s=linspace(s(2),M,n);
th=linspace(a2,a3,n);


[S2,TH2]=meshgrid(s,th);
%Calculation of the vector field of rotations in the coordinates s and th.
%The variables (t,phi) denote geodesic polar coordinate around the axis of
%rotation. 


t=zeros(n,n);
phi=zeros(n,n);
vecS2=zeros(n,n);
vecTH2=zeros(n,n);

%This loop constructs the components of the vector field of rotations
%in the the s, th coordinates. The function InverseSTH determines the
%coordinates of (t,phi) given coordiantes s and th and is necessary to
%calculate the components of the vector field. 
for i=1:n
    for j=1:n,
        [t(i,j),phi(i,j)]=InverseSTH(alpha,beta,R,S2(i,j),TH2(i,j)); 
        
            vecS2(i,j)=-cos(S2(i,j)/R)*cos(beta)*sin(phi(i,j))*cos(TH2(i,j)-alpha)+...
                cos(S2(i,j)/R)*cos(phi(i,j))*(sin(TH2(i,j)-alpha))-...
                sin(S2(i,j)/R)*sin(phi(i,j))*sin(beta);
            
            vecS2(i,j)=R*sin(t(i,j)/R)*vecS2(i,j);
            
            vecTH2(i,j)=(cos(beta)*sin(phi(i,j))*sin(TH2(i,j)-alpha)+...
                cos(phi(i,j))*cos(alpha-TH2(i,j)));
        
            vecTH2(i,j)=R*sin(t(i,j)/R)*vecTH2(i,j);
            vecTH2(i,j)=1/(R*sin(S2(i,j)/R))*vecTH2(i,j);
            
    end
end
vecTH2=vecTH2.*polarity;

%Calcution of derivatives
ds=s(2)-s(1);
DS=DiffX(n)/ds;
RHOS2=(DS*RHO2')';
FS2=(DS*F2')';

dth=th(2)-th(1);
DTH=DiffX(n)/dth;
RHOTH2=(DTH*RHO2);
FTH2=(DTH*F2);

%Calculating the vector fields from the formula for the push forward of the
%derivative.
vecU2=vecS2.*(RHOS2.*cos(F2)-RHO2.*FS2.*sin(F2))+vecTH2.*(RHOTH2.*cos(F2)-RHO2.*FTH2.*sin(F2));
vecV2=vecS2.*(RHOS2.*sin(F2)+RHO2.*FS2.*cos(F2))+vecTH2.*(RHOTH2.*sin(F2)+RHO2.*FTH2.*cos(F2));

%Plotting vector field in the sector. We plot the physical retina so we can
%see where the vectors are located on the physical retina.
figure(1);
hold on;

C=abs(R.^2*sin(S2./R).^2-S2.^2);
surf(RHO2.*cos(F2),RHO2.*sin(F2),zeros(n,n),C);
view(0,90);
caxis([0,1]);
colorbar;
colormap('bone');
shading interp;

U2=RHO2.*cos(F2);
V2=RHO2.*sin(F2);
quiver(U2,V2,vecU2,vecV2,'b');
axis square;
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

plot(beta*R*cos(alpha),beta*R*sin(alpha),'o');
xlabel('u');
ylabel('v');

hold off;


%% Sector 3
%This part of the code computes the vector fields in the sectors. We
%read in the numerical solution computed using the script Mapping.m to
%compute the value of rho. The vector fields are then found by computing
%the pushforward in the coordinates (u,v).

%Reading in data;
Sec3Data=dlmread('Sector3'); %radial stretching factor.
n=size(Sec3Data,1);
RHO3=Sec3Data(1:n,1:n);
F3=Sec3Data(1:n,n+1:2*n);

%Record the number of points that were discretized over.

% Construction of domain in s th coordinates
s=linspace(0,M,n+1);
s=linspace(s(2),M,n);
th=linspace(a3,a4,n);


[S3,TH3]=meshgrid(s,th);
% Calculation of the vector field of rotations in the coordinates s and th.
% The variables (t,phi) denote geodesic polar coordinate around the axis of
% rotation. 


t=zeros(n,n);
phi=zeros(n,n);
vecS3=zeros(n,n);
vecTH3=zeros(n,n);

%This loop constructs the components of the vector field of rotations
%in the the s, th coordinates. The function InverseSTH determines the
%coordinates of (t,phi) given coordiantes s and th and is necessary to
%calculate the components of the vector field. 
for i=1:n
    for j=1:n,
        [t(i,j),phi(i,j)]=InverseSTH(alpha,beta,R,S3(i,j),TH3(i,j)); 
        
            vecS3(i,j)=-cos(S3(i,j)/R)*cos(beta)*sin(phi(i,j))*cos(TH3(i,j)-alpha)+...
                cos(S3(i,j)/R)*cos(phi(i,j))*(sin(TH3(i,j)-alpha))-...
                sin(S3(i,j)/R)*sin(phi(i,j))*sin(beta);
            
            vecS3(i,j)=R*sin(t(i,j)/R)*vecS3(i,j);
            
            vecTH3(i,j)=(cos(beta)*sin(phi(i,j))*sin(TH3(i,j)-alpha)+...
                cos(phi(i,j))*cos(alpha-TH3(i,j)));
        
            vecTH3(i,j)=R*sin(t(i,j)/R)*vecTH3(i,j);
            vecTH3(i,j)=1/(R*sin(S3(i,j)/R))*vecTH3(i,j);
            
    end
end
vecTH3=vecTH3.*polarity;

%Calcution of derivatives
ds=s(2)-s(1);
DS=DiffX(n)/ds;
RHOS3=(DS*RHO3')';
FS3=(DS*F3')';

dth=th(2)-th(1);
DTH=DiffX(n)/dth;
RHOTH3=(DTH*RHO3);
FTH3=(DTH*F3);

%Calculating the vector fields from the formula for the push forward of the
%derivative.
vecU3=vecS3.*(RHOS3.*cos(F3)-RHO3.*FS3.*sin(F3))+vecTH3.*(RHOTH3.*cos(F3)-RHO3.*FTH3.*sin(F3));
vecV3=vecS3.*(RHOS3.*sin(F3)+RHO3.*FS3.*cos(F3))+vecTH3.*(RHOTH3.*sin(F3)+RHO3.*FTH3.*cos(F3));

%Plotting vector field in the sector. We plot the physical retina so we can
%see where the vectors are located on the physical retina.
figure(1);
hold on;

C=abs(R.^3*sin(S3./R).^2-S3.^2);
surf(RHO3.*cos(F3),RHO3.*sin(F3),zeros(n,n),C);
view(0,90);
caxis([0,1]);
colorbar;
colormap('bone');
shading interp;

U3=RHO3.*cos(F3);
V3=RHO3.*sin(F3);
quiver(U3,V3,vecU3,vecV3,'b');
axis square;
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

plot(beta*R*cos(alpha),beta*R*sin(alpha),'o');
xlabel('u');
ylabel('v');

hold off;

%% Sector 4
%This part of the code computes the vector fields in the sectors. We
%read in the numerical solution computed using the script Mapping.m to
%compute the value of rho. The vector fields are then found by computing
%the pushforward in the coordinates (u,v).

%Reading in data;
Sec4Data=dlmread('Sector4'); %radial stretching factor.
n=size(Sec4Data,1);
RHO4=Sec4Data(1:n,1:n);
F4=Sec4Data(1:n,n+1:2*n);

%Record the number of points that were discretized over.

% Construction of domain in s th coordinates
s=linspace(0,M,n+1);
s=linspace(s(2),M,n);
th=linspace(a4,a1+2*pi,n);


[S4,TH4]=meshgrid(s,th);
%Calculation of the vector field of rotations in the coordinates s and th.
%The variables (t,phi) denote geodesic polar coordinate around the axis of
%rotation. 


t=zeros(n,n);
phi=zeros(n,n);
vecS=zeros(n,n);
vecTH=zeros(n,n);

%This loop constructs the components of the vector field of rotations
%in the the s, th coordinates. The function InverseSTH determines the
%coordinates of (t,phi) given coordiantes s and th and is necessary to
%calculate the components of the vector field. 
for i=1:n
    for j=1:n,
        [t(i,j),phi(i,j)]=InverseSTH(alpha,beta,R,S4(i,j),TH4(i,j)); 
        
            vecS4(i,j)=-cos(S4(i,j)/R)*cos(beta)*sin(phi(i,j))*cos(TH4(i,j)-alpha)+...
                cos(S4(i,j)/R)*cos(phi(i,j))*(sin(TH4(i,j)-alpha))-...
                sin(S4(i,j)/R)*sin(phi(i,j))*sin(beta);
            
            vecS4(i,j)=R*sin(t(i,j)/R)*vecS4(i,j);
            
            vecTH4(i,j)=(cos(beta)*sin(phi(i,j))*sin(TH4(i,j)-alpha)+...
                cos(phi(i,j))*cos(alpha-TH4(i,j)));
        
            vecTH4(i,j)=R*sin(t(i,j)/R)*vecTH4(i,j);
            vecTH4(i,j)=1/(R*sin(S4(i,j)/R))*vecTH4(i,j);
            
    end
end
vecTH4=vecTH4.*polarity;

%Calcution of derivatives
ds=s(2)-s(1);
DS=DiffX(n)/ds;
RHOS4=(DS*RHO4')';
FS4=(DS*F4')';

dth=th(2)-th(1);
DTH=DiffX(n)/dth;
RHOTH4=(DTH*RHO4);
FTH4=(DTH*F4);

%Calculating the vector fields from the formula for the push forward of the
%derivative.
vecU4=vecS4.*(RHOS4.*cos(F4)-RHO4.*FS4.*sin(F4))+vecTH4.*(RHOTH4.*cos(F4)-RHO4.*FTH4.*sin(F4));
vecV4=vecS4.*(RHOS4.*sin(F4)+RHO4.*FS4.*cos(F4))+vecTH4.*(RHOTH4.*sin(F4)+RHO4.*FTH4.*cos(F4));

%Plotting vector field in the sector. We plot the physical retina so we can
%see where the vectors are located on the physical retina.
figure(1);
hold on;

C=abs(R.^3*sin(S4./R).^2-S4.^2);
surf(RHO4.*cos(F4),RHO4.*sin(F4),zeros(n,n),C);
view(0,90);
caxis([0,1]);
colorbar;
colormap('bone');
shading interp;

U4=RHO4.*cos(F4);
V4=RHO4.*sin(F4);
quiver(U4,V4,vecU4,vecV4,'b');
axis square;
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

plot(beta*R*cos(alpha),beta*R*sin(alpha),'o');
xlabel('u');
ylabel('v');

hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Import In Data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DirectionSelective Cells
[filedir,go]=bbdirselector('select data folder',cd);
if ~go disp('no folder selected'); return; end

namesMap=getFileList(filedir,'map_',0,'anywhere');
for i=1:size(namesMap,2)   
    if strfind(namesMap{i},'.mat')
    load(num2str(namesMap{i}));
    end
end

%load('map_May_27_2014_2.mat');
data=map.SzMapOutputDS;

U=data(:,1);
V=data(:,2);
vecU=data(:,3);
vecV=data(:,4);
Sec=data(:,5);

switch typeStr
    case 'ONDS' 
type=map.ONDS;
    case 'ONOFFDS' 
type=map.ONOFFDS;
    case 'OFFDS' 
type=map.OFFDS;
    case 'retroONDS' 
type=map.retroONDS;
end

n=size(U,1);

figure(1);
hold on;
% quiver(U,V,vecU,vecV);
% 
% hold off;

for i=1:n,
    
    if type(i)==1
        j=j+1;   
    if Sec(i)==1,
        
        
        
        temp=(U(i)-U1).^2+(V(i)-V1).^2;
        [minval,ind]=min(temp(:));
        [a,b]=ind2sub([size(temp,1),size(temp,2)],ind);
        
        
        [g1,g2] = InverseVec(vecU(i),vecV(i),RHO1(a,b),F1(a,b),RHOS1(a,b),FS1(a,b),RHOTH1(a,b),FTH1(a,b));
        Sn(i)=S1(a,b);
        THn(i)=TH1(a,b);
        vecSn(i)=g1;
        vecTHn(i)=g2;
        
        
        
    elseif Sec(i)==2,
        
        temp=(U(i)-U1).^2+(V(i)-V1).^2;
        [minval,ind]=min(temp(:));
        [a,b]=ind2sub([size(temp,1),size(temp,2)],ind);
        
        [g1,g2] = InverseVec(vecU(i),vecV(i),RHO2(a,b),F2(a,b),RHOS2(a,b),FS2(a,b),RHOTH2(a,b),FTH2(a,b));
        Sn(i)=S2(a,b);
        THn(i)=TH1(a,b);
        vecSn(i)=g1;
        vecTHn(i)=g2;
        
        
    elseif Sec(i)==3,
        temp=(U(i)-U1).^2+(V(i)-V1).^2;
        [minval,ind]=min(temp(:));
        [a,b]=ind2sub([size(temp,1),size(temp,2)],ind);
        
        [g1,g2] = InverseVec(vecU(i),vecV(i),RHO3(a,b),F3(a,b),RHOS3(a,b),FS3(a,b),RHOTH3(a,b),FTH3(a,b));
        Sn(i)=S3(a,b);
        THn(i)=TH3(a,b);
        vecSn(i)=g1;
        vecTHn(i)=g2;
    else
        temp=(U(i)-U1).^2+(V(i)-V1).^2;
        [minval,ind]=min(temp(:));
        [a,b]=ind2sub([size(temp,1),size(temp,2)],ind);
        
        [g1,g2] = InverseVec(vecU(i),vecV(i),RHO4(a,b),F4(a,b),RHOS4(a,b),FS4(a,b),RHOTH4(a,b),FTH4(a,b));
        Sn(i)=S4(a,b);
        THn(i)=TH4(a,b);
        vecSn(i)=g1;
        vecTHn(i)=g2;
        
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Plot of Data in STH Coordinates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


vecSntemp=vecSn./(vecSn.^2+vecTHn.^2).^(1/2);
vecTHntemp=vecTHn./(vecSn.^2+vecTHn.^2).^(1/2);

vecSn=vecSntemp;
vecTHn=vecTHntemp;

figure(2);
hold on;
quiver(Sn,THn/(2*pi),vecSn,vecTHn,1,'linewidth',1);

Sn=Sn';
THn=THn';
vecSn=vecSn';
vecTHn=vecTHn';

dlmwrite('STHData',[Sn,THn,vecSn,vecTHn]);


hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Plot of Data in Generalized STH Coordinates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3);

hold on;

Xtemp=Sn.*cos(THn);
Ytemp=Sn.*sin(THn);

vecXtemp=cos(THn).*vecSn-Sn.*sin(THn).*vecTHn;
vecYtemp=sin(THn).*vecSn+Sn.*cos(THn).*vecTHn;

a=vecXtemp./(vecXtemp.^2+vecYtemp.^2).^(1/2);
b=vecYtemp./(vecXtemp.^2+vecYtemp.^2).^(1/2);

vecXtemp=a;
vecYtemp=b;

quiver(Xtemp,Ytemp,vecXtemp,vecYtemp);
plot(R*beta*cos(alpha),R*beta*sin(alpha),'o');
axis([-1,1,-1,1]);
axis square;

xlabel('x^{\prime}');
ylabel('y^{\prime}');

hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Plot of Data on Retina
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4);

hold on;

Xtemp=R.*cos(THn).*sin(Sn/R);
Ytemp=R.*sin(THn).*sin(Sn/R);
Ztemp=R-R.*cos(Sn/R);


vecXtemp=cos(THn).*cos(Sn)-R*sin(THn).*sin(Sn);
vecYtemp=sin(THn).*cos(Sn)+R*cos(THn).*sin(Sn);
vecZtemp=sin(Sn);

a=vecXtemp./(vecXtemp.^2+vecYtemp.^2+vecZtemp.^2).^(1/2);
b=vecYtemp./(vecXtemp.^2+vecYtemp.^2+vecZtemp.^2).^(1/2);
c=vecZtemp./(vecXtemp.^2+vecYtemp.^2+vecZtemp.^2).^(1/2);
vecXtemp=a;
vecYtemp=b;
vecZtemp=c;

quiver3(Xtemp,Ytemp,Ztemp,vecXtemp,vecYtemp,vecZtemp);

hold off;
end


