function []= VectorFieldGeneralizedForGui(alpha, beta, polarity, discPoints,fieldType)


% VectorField Script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This script computes the vector field on the flattened state for a
%   vector field corresponding to an axis of rotation about a point with
%   azimuthal angle alpha and altitude angle beta as measured from the south
%   pole. 
%
%   Input Parameters:
%   1. R - Radius of the hemisphere. 
%   2. M - Length of meridean connecting the south pole of the retina to
%   the boundary of S. This should be normalized to 1.
%   3. mi - Length of meridan connecting the south pole to cut i.
%   4. nu - Poisson ratio of the material.
%   5. ai - angular position of the cuts.
%   6. n - number of points used in discretization.
%   7. alpha - azimuthal angle of axis of rotation.
%   8. beta - altitude angle of axis of rotation. 
%
%   Independent Variables:
%   1. s - The arclength coordinate along meridians.
%   2. th - The azimuthal angle along the sphere.
%
%   Dependent Variables:
%   1. vecU - component of flattened vector field in the U-direction.
%   2. vecV - component of flattened vector field in the V-direction.
%   3. vecS - component of flattened vector field in the S-direction.
%   4. vecTH - component of flattened vector field in the TH-direction.
%
%   Ouput:
%   1. vecU - component of flattened vector field in the U-direction. We
%   will output this as a function of polar coordinates.
%   2. vecV - component of flattened vector field in the V-direction. We
%   will output this as a function of polar coordinates.
%   3. vecS - component of flattened vector field in the S-direction.
%   4. vecTH - component of flattened vector field in the TH-direction.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all;
figure;
% Input Parameters

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

%% Computing coordinates of center of ration.

surfcolor='white';
vColor='r';

%% Sector 1
%This part of the code computes the vector fields in the sectors. We
%read in the numerical solution computed using the script Mapping.m to
%compute the value of rho. The vector fields are then found by computing
%the pushforward in the coordinates (u,v).

num2str(discPoints)

%Reading in data;
Sec1Data=dlmread(['Sector1_',num2str(discPoints)]); %radial stretching factor.
n=size(Sec1Data,1); %Record the number of points that were discretized over.
RHO=Sec1Data(1:n,1:n);
F=Sec1Data(1:n,n+1:2*n);


% Construction of domain in s-th coordinates
s=linspace(0,M,n+1);
s=linspace(s(2),M,n);
th=linspace(a1,a2,n);
[S,TH]=meshgrid(s,th);

%Construction of domain in visual optic space. This corresponds to the
%the coordinates on the hemisphere about the optical axis which the eye can
%see.

SV=-pi*R/(2*M)*S+pi*R;
THV=TH+pi;

%Calculation of the vector field of rotations in the coordinates svis and
%thvis. The variables (t,phi) denote geodesic polar coordinate around the axis of
%rotation. 

t=zeros(n,n);
phi=zeros(n,n);
vecS=zeros(n,n);
vecTH=zeros(n,n);


%This loop constructs the components of the vector field of rotations
%in the visual space of the eye.  The function InverseSTH determines the
%coordinates of (t,phi) given coordiantes svis and thvis and is necessary to
%calculate the components of the vector field. 

for i=1:n
    for j=1:n,
        [t(i,j),phi(i,j)]=InverseSTH(alpha,beta,R,SV(i,j),THV(i,j)); 
        
            vecSV(i,j)=-cos(SV(i,j)/R)*cos(beta)*sin(phi(i,j))*cos(THV(i,j)-alpha)+...
                cos(SV(i,j)/R)*cos(phi(i,j))*(sin(THV(i,j)-alpha))-...
                sin(SV(i,j)/R)*sin(phi(i,j))*sin(beta);
            
            vecSV(i,j)=R*sin(t(i,j)/R)*vecSV(i,j);
            
            vecTHV(i,j)=(cos(beta)*sin(phi(i,j))*sin(THV(i,j)-alpha)+...
                cos(phi(i,j))*cos(alpha-THV(i,j)));
        
            vecTHV(i,j)=R*sin(t(i,j)/R)*vecTHV(i,j);
            vecTHV(i,j)=1/(R*sin(SV(i,j)/R))*vecTHV(i,j);
            
    end
end

vecTHV=vecTHV.*polarity;

%Conversion of the vector fields in visual space to vector fields on the
%retina.

vecS=-2*M/(pi*R)*vecSV;
vecTH=vecTHV;

%Calcution of derivatives
ds=s(2)-s(1);
DS=DiffX(n)/ds;
RHOS=(DS*RHO')';
FS=(DS*F')';

dth=th(2)-th(1);
DTH=DiffX(n)/dth;
RHOTH=(DTH*RHO);
FTH=(DTH*F);

%Calculating the vector fields from the formula for the push forward of the
%derivative.
vecU=vecS.*(RHOS.*cos(F)-RHO.*FS.*sin(F))+vecTH.*(RHOTH.*cos(F)-RHO.*FTH.*sin(F));
vecV=vecS.*(RHOS.*sin(F)+RHO.*FS.*cos(F))+vecTH.*(RHOTH.*sin(F)+RHO.*FTH.*cos(F));

%Plotting vector field in the sector. We plot the physical retina so we can
%see where the vectors are located on the physical retina.
figure(1);
hold on;

C=abs(R.^2*sin(S./R).^2-S.^2);
surf(RHO.*cos(F),RHO.*sin(F),zeros(n,n),C);
view(0,90);
colormap(jet);
caxis([0,1]);
colorbar;
colormap(surfcolor);
shading interp;

U=RHO.*cos(F);
V=RHO.*sin(F);

quiver(U,V,vecU,vecV,vColor);
axis square;
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
plot(beta*R*cos(alpha),beta*R*sin(alpha),'o');
xlabel('u');
ylabel('v');

hold off;

%Plotting vector field in (S,TH) coordiantes;
figure(2);
hold on;
vecStemp=vecS./(vecS.^2+vecTH.^2).^(1/2);
vecTHtemp=vecTH./(vecS.^2+vecTH.^2).^(1/2);

vecS=vecStemp;
vecTH=vecTHtemp;

quiver(S,TH/(2*pi),vecS,vecTH,vColor);
plot(beta*R,alpha/(2*pi),'o');
axis([0,1,0,1]);
axis square;

xlabel('s');
ylabel('\theta');

set(gcf,'units','normalized','outerposition',[0 0 1 1]);

hold off;

%Plotting Vector Field in Generalized (S,TH) coordinates
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

quiver(X,Y,vecX,vecY,vColor);
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

%Plotting Vector Fields on the sphere
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

%ntemp=40;

s=linspace(0,M,ntemp);
th=linspace(a1,a2,ntemp);
[S,TH]=meshgrid(s,th);

surf(R*sin(S/R).*cos(TH),R*sin(S/R).*sin(TH),R-R*cos(S/R),S,'EdgeColor','none','FaceAlpha','interp','AlphaDataMapping','scaled',...
    'AlphaData',S);
colormap(surfcolor);
alphamap('default')
view(104,6);

quiver3(X,Y,Z,vecX,vecY,vecZ,vColor);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

axis([-R,R,-R,R,0,2*R]);
axis square;

plot3(R*cos(alpha)*sin(beta),R*sin(alpha)*sin(beta),R-R*cos(beta),'o');
xlabel('x');
ylabel('y');
zlabel('z');
hold off;

%% Plotting vector field in standard retina coordinates

StSec1Data=dlmread(['StandardSector1_',num2str(discPoints)]);

n=size(StSec1Data,1);

StRHO=StSec1Data(1:n,1:n);
StF=StSec1Data(1:n,n+1:2*n);
StS=StSec1Data(1:n,2*n+1:3*n);
StTH=StSec1Data(1:n,3*n+1:4*n);

StU=StRHO.*cos(StF);
StV=StRHO.*sin(StF);

%Calcution of derivatives
ds=StS(1,2)-StS(1,1);
DS=DiffX(n)/ds;
StRHOS=(DS*StRHO')';
StFS=(DS*StF')';

dth=StTH(2,1)-StTH(1,1);
DTH=DiffX(n)/dth;
StRHOTH=(DTH*StRHO);
StFTH=(DTH*StF);

StVecU=vecS.*(StRHOS.*cos(StF)-StRHO.*StFS.*sin(StF))+vecTH.*(StRHOTH.*cos(StF)-StRHO.*StFTH.*sin(StF));
StVecV=vecS.*(StRHOS.*sin(StF)+StRHO.*StFS.*cos(StF))+vecTH.*(StRHOTH.*sin(StF)+StRHO.*StFTH.*cos(StF));

figure(5);
hold on;
temp=(StVecU.^2+StVecV.^2).^(1/2);

StVecU=StVecU./temp;
StVecV=StVecV./temp;

quiver(StU,StV,StVecU,StVecV,'r');
axis square;
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

%Writing Vector Field Data.
dlmwrite(['Sec1VecField_',num2str(discPoints),'_',fieldType],[U,V,vecU,vecV,RHO,F,RHOS,RHOTH,FS,FTH,S,TH,vecS,vecTH,Xprime,Yprime,vecXprime,vecYprime,X,Y,Z,vecX,vecY,vecZ]);


% Sector 2
%This part of the code computes the vector fields in the sectors. We
%read in the numerical solution computed using the script Mapping.m to
%compute the value of rho. The vector fields are then found by computing
%the pushforward in the coordinates (u,v).

%Reading in data;
Sec2Data=dlmread(['Sector2_',num2str(discPoints)]); %radial stretching factor.
n=size(Sec2Data,1);
RHO=Sec2Data(1:n,1:n);
F=Sec2Data(1:n,n+1:2*n);

 %Record the number of points that were discretized over.

% Construction of domain in s-th coordinates
s=linspace(0,M,n+1);
s=linspace(s(2),M,n);
th=linspace(a2,a3,n);
[S,TH]=meshgrid(s,th);

%Construction of domain in visual optic space. This corresponds to the
%the coordinates on the hemisphere about the optical axis which the eye can
%see.

SV=-pi*R/(2*M)*S+pi*R;
THV=TH+pi;

%Calculation of the vector field of rotations in the coordinates svis and
%thvis. The variables (t,phi) denote geodesic polar coordinate around the axis of
%rotation. 

t=zeros(n,n);
phi=zeros(n,n);
vecS=zeros(n,n);
vecTH=zeros(n,n);


%This loop constructs the components of the vector field of rotations
%in the visual space of the eye.  The function InverseSTH determines the
%coordinates of (t,phi) given coordiantes svis and thvis and is necessary to
%calculate the components of the vector field. 

for i=1:n
    for j=1:n,
        [t(i,j),phi(i,j)]=InverseSTH(alpha,beta,R,SV(i,j),THV(i,j)); 
        
            vecSV(i,j)=-cos(SV(i,j)/R)*cos(beta)*sin(phi(i,j))*cos(THV(i,j)-alpha)+...
                cos(SV(i,j)/R)*cos(phi(i,j))*(sin(THV(i,j)-alpha))-...
                sin(SV(i,j)/R)*sin(phi(i,j))*sin(beta);
            
            vecSV(i,j)=R*sin(t(i,j)/R)*vecSV(i,j);
            
            vecTHV(i,j)=(cos(beta)*sin(phi(i,j))*sin(THV(i,j)-alpha)+...
                cos(phi(i,j))*cos(alpha-THV(i,j)));
        
            vecTHV(i,j)=R*sin(t(i,j)/R)*vecTHV(i,j);
            vecTHV(i,j)=1/(R*sin(SV(i,j)/R))*vecTHV(i,j);
            
    end
end

vecTHV=vecTHV.*polarity;

%Conversion of the vector fields in visual space to vector fields on the
%retina.

vecS=-2*M/(pi*R)*vecSV;
vecTH=vecTHV;

%Calcution of derivatives
ds=s(2)-s(1);
DS=DiffX(n)/ds;
RHOS=(DS*RHO')';
FS=(DS*F')';

dth=th(2)-th(1);
DTH=DiffX(n)/dth;
RHOTH=(DTH*RHO);
FTH=(DTH*F);

%Calculating the vector fields from the formula for the push forward of the
%derivative.
vecU=vecS.*(RHOS.*cos(F)-RHO.*FS.*sin(F))+vecTH.*(RHOTH.*cos(F)-RHO.*FTH.*sin(F));
vecV=vecS.*(RHOS.*sin(F)+RHO.*FS.*cos(F))+vecTH.*(RHOTH.*sin(F)+RHO.*FTH.*cos(F));

%Plotting vector field in the sector. We plot the physical retina so we can
%see where the vectors are located on the physical retina.
figure(1);
hold on;

C=abs(R.^2*sin(S./R).^2-S.^2);
surf(RHO.*cos(F),RHO.*sin(F),zeros(n,n),C);
view(0,90);
colormap(jet);
caxis([0,1]);
colorbar;
colormap(surfcolor);
shading interp;

U=RHO.*cos(F);
V=RHO.*sin(F);
quiver(U,V,vecU,vecV,vColor);
axis square;
axis ([-1,1,-1,1]);

%Plotting vector field in (S,TH) coordiantes;
figure(2);
hold on;
vecStemp=vecS./(vecS.^2+vecTH.^2).^(1/2);
vecTHtemp=vecTH./(vecS.^2+vecTH.^2).^(1/2);

vecS=vecStemp;
vecTH=vecTHtemp;

quiver(S,TH/(2*pi),vecS,vecTH,vColor);
plot(beta*R,alpha/(2*pi),'o');
axis([0,1,0,1]);
axis square;

set(gcf,'units','normalized','outerposition',[0 0 1 1]);

hold off;

%Plotting Vector Field in Generalized S TH coordinates
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

quiver(X,Y,vecX,vecY,vColor);
plot(R*beta*cos(alpha),R*beta*sin(alpha),'o');
axis([-1,1,-1,1]);
axis square;

xlabel('x');
ylabel('y');

hold off;

set(gcf,'units','normalized','outerposition',[0 0 1 1]);

hold off;

Xprime=X;
Yprime=Y;

vecXprime=vecX;
vecYprime=vecY;
%Plotting Vector Fields on the sphere
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

%ntemp=40;

s=linspace(0,M,ntemp);
th=linspace(a2,a3,ntemp);
[S,TH]=meshgrid(s,th);

surf(R*sin(S/R).*cos(TH),R*sin(S/R).*sin(TH),R-R*cos(S/R),S,'EdgeColor','none','FaceAlpha','interp','AlphaDataMapping','scaled',...
    'AlphaData',S);
colormap(surfcolor);
alphamap('default')
view(104,6);

quiver3(X,Y,Z,vecX,vecY,vecZ,vColor);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

hold off;

%% Plotting vector field in standard retina coordinates
StSec2Data=dlmread(['StandardSector2_',num2str(discPoints)]);

StRHO=StSec2Data(1:n,1:n);
StF=StSec2Data(1:n,n+1:2*n);
StS=StSec2Data(1:n,2*n+1:3*n);
StTH=StSec2Data(1:n,3*n+1:4*n);

StU=StRHO.*cos(StF);
StV=StRHO.*sin(StF);

%Calcution of derivatives
ds=StS(1,2)-StS(1,1);
DS=DiffX(n)/ds;
StRHOS=(DS*StRHO')';
StFS=(DS*StF')';

dth=StTH(2,1)-StTH(1,1);
DTH=DiffX(n)/dth;
StRHOTH=(DTH*StRHO);
StFTH=(DTH*StF);

StVecU=vecS.*(StRHOS.*cos(StF)-StRHO.*StFS.*sin(StF))+vecTH.*(StRHOTH.*cos(StF)-StRHO.*StFTH.*sin(StF));
StVecV=vecS.*(StRHOS.*sin(StF)+StRHO.*StFS.*cos(StF))+vecTH.*(StRHOTH.*sin(StF)+StRHO.*StFTH.*cos(StF));

figure(5);
temp=(StVecU.^2+StVecV.^2).^(1/2);

StVecU=StVecU./temp;
StVecV=StVecV./temp;

quiver(StU,StV,StVecU,StVecV,'r');


%Writing Vector Field Data.
dlmwrite(['Sec2VecField_',num2str(discPoints),'_',fieldType],[U,V,vecU,vecV,RHO,F,RHOS,RHOTH,FS,FTH,S,TH,vecS,vecTH,Xprime,Yprime,vecXprime,vecYprime,X,Y,Z,vecX,vecY,vecZ]);


% Sector 3
%This part of the code computes the vector fields in the sectors. We
%read in the numerical solution computed using the script Mapping.m to
%compute the value of rho. The vector fields are then found by computing
%the pushforward in the coordinates (u,v).

%Reading in data;
Sec3Data=dlmread(['Sector3_',num2str(discPoints)]); %radial stretching factor.
n=size(Sec3Data,1);
RHO=Sec3Data(1:n,1:n);
F=Sec3Data(1:n,n+1:2*n);

 %Record the number of points that were discretized over.

% Construction of domain in s-th coordinates
s=linspace(0,M,n+1);
s=linspace(s(2),M,n);
th=linspace(a3,a4,n);
[S,TH]=meshgrid(s,th);

%Construction of domain in visual optic space. This corresponds to the
%the coordinates on the hemisphere about the optical axis which the eye can
%see.

SV=-pi*R/(2*M)*S+pi*R;
THV=TH+pi;

%Calculation of the vector field of rotations in the coordinates svis and
%thvis. The variables (t,phi) denote geodesic polar coordinate around the axis of
%rotation. 

t=zeros(n,n);
phi=zeros(n,n);
vecS=zeros(n,n);
vecTH=zeros(n,n);


%This loop constructs the components of the vector field of rotations
%in the visual space of the eye.  The function InverseSTH determines the
%coordinates of (t,phi) given coordiantes svis and thvis and is necessary to
%calculate the components of the vector field. 

for i=1:n
    for j=1:n,
        [t(i,j),phi(i,j)]=InverseSTH(alpha,beta,R,SV(i,j),THV(i,j)); 
        
            vecSV(i,j)=-cos(SV(i,j)/R)*cos(beta)*sin(phi(i,j))*cos(THV(i,j)-alpha)+...
                cos(SV(i,j)/R)*cos(phi(i,j))*(sin(THV(i,j)-alpha))-...
                sin(SV(i,j)/R)*sin(phi(i,j))*sin(beta);
            
            vecSV(i,j)=R*sin(t(i,j)/R)*vecSV(i,j);
            
            vecTHV(i,j)=(cos(beta)*sin(phi(i,j))*sin(THV(i,j)-alpha)+...
                cos(phi(i,j))*cos(alpha-THV(i,j)));
        
            vecTHV(i,j)=R*sin(t(i,j)/R)*vecTHV(i,j);
            vecTHV(i,j)=1/(R*sin(SV(i,j)/R))*vecTHV(i,j);
            
    end
end

vecTHV=vecTHV.*polarity;

%Conversion of the vector fields in visual space to vector fields on the
%retina.

vecS=-2*M/(pi*R)*vecSV;
vecTH=vecTHV;

%Calcution of derivatives
ds=s(2)-s(1);
DS=DiffX(n)/ds;
RHOS=(DS*RHO')';
FS=(DS*F')';

dth=th(2)-th(1);
DTH=DiffX(n)/dth;
RHOTH=(DTH*RHO);
FTH=(DTH*F);

%Calculating the vector fields from the formula for the push forward of the
%derivative.
vecU=vecS.*(RHOS.*cos(F)-RHO.*FS.*sin(F))+vecTH.*(RHOTH.*cos(F)-RHO.*FTH.*sin(F));
vecV=vecS.*(RHOS.*sin(F)+RHO.*FS.*cos(F))+vecTH.*(RHOTH.*sin(F)+RHO.*FTH.*cos(F));

%Plotting vector field in the sector. We plot the physical retina so we can
%see where the vectors are located on the physical retina.
figure(1);
hold on;

C=abs(R.^2*sin(S./R).^2-S.^2);
surf(RHO.*cos(F),RHO.*sin(F),zeros(n,n),C);
view(0,90);
colormap(jet);
caxis([0,1]);
colorbar;
colormap('jet');
shading interp;

U=RHO.*cos(F);
V=RHO.*sin(F);
quiver(U,V,vecU,vecV,vColor);
axis square;

xlabel('u');
ylabel('v');
hold off;

%Plotting vector field in (S,TH) coordiantes;
figure(2);
hold on;
vecStemp=vecS./(vecS.^2+vecTH.^2).^(1/2);
vecTHtemp=vecTH./(vecS.^2+vecTH.^2).^(1/2);

vecS=vecStemp;
vecTH=vecTHtemp;

quiver(S,TH/(2*pi),vecS,vecTH,vColor);
plot(beta*R,alpha/(2*pi),'o');
axis([0,1,0,1]);
axis square;

set(gcf,'units','normalized','outerposition',[0 0 1 1]);

hold off;

%Plotting Vector Field in Generalized S TH coordinates
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

quiver(X,Y,vecX,vecY,vColor);
plot(R*beta*cos(alpha),R*beta*sin(alpha),'o');
axis([-1,1,-1,1]);
axis square;

xlabel('x');
ylabel('y');

hold off;

set(gcf,'units','normalized','outerposition',[0 0 1 1]);

hold off;

Xprime=X;
Yprime=Y;

vecXprime=vecX;
vecYprime=vecY;

%Plotting Vector Fields on the sphere
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

%ntemp=40;

s=linspace(0,M,ntemp);
th=linspace(a3,a4,ntemp);
[S,TH]=meshgrid(s,th);

surf(R*sin(S/R).*cos(TH),R*sin(S/R).*sin(TH),R-R*cos(S/R),S,'EdgeColor','none','FaceAlpha','interp','AlphaDataMapping','scaled',...
    'AlphaData',S);
colormap(surfcolor);
alphamap('default')
view(104,6);

quiver3(X,Y,Z,vecX,vecY,vecZ,vColor);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

hold off;

%% Plotting vector field in standard retina coordinates
StSec3Data=dlmread(['StandardSector3_',num2str(discPoints)]);

StRHO=StSec3Data(1:n,1:n);
StF=StSec3Data(1:n,n+1:2*n);
StS=StSec3Data(1:n,2*n+1:3*n);
StTH=StSec3Data(1:n,3*n+1:4*n);

StU=StRHO.*cos(StF);
StV=StRHO.*sin(StF);

%Calcution of derivatives
ds=StS(1,2)-StS(1,1);
DS=DiffX(n)/ds;
StRHOS=(DS*StRHO')';
StFS=(DS*StF')';

dth=StTH(2,1)-StTH(1,1);
DTH=DiffX(n)/dth;
StRHOTH=(DTH*StRHO);
StFTH=(DTH*StF);

StVecU=vecS.*(StRHOS.*cos(StF)-StRHO.*StFS.*sin(StF))+vecTH.*(StRHOTH.*cos(StF)-StRHO.*StFTH.*sin(StF));
StVecV=vecS.*(StRHOS.*sin(StF)+StRHO.*StFS.*cos(StF))+vecTH.*(StRHOTH.*sin(StF)+StRHO.*StFTH.*cos(StF));

figure(5);
temp=(StVecU.^2+StVecV.^2).^(1/2);

StVecU=StVecU./temp;
StVecV=StVecV./temp;

quiver(StU,StV,StVecU,StVecV,'r');

%Writing Vector Field Data.
dlmwrite(['Sec3VecField_',num2str(discPoints),'_',fieldType],[U,V,vecU,vecV,RHO,F,RHOS,RHOTH,FS,FTH,S,TH,vecS,vecTH,Xprime,Yprime,vecXprime,vecYprime,X,Y,Z,vecX,vecY,vecZ]);


% Sector 4
%This part of the code computes the vector fields in the sectors. We
%read in the numerical solution computed using the script Mapping.m to
%compute the value of rho. The vector fields are then found by computing
%the pushforward in the coordinates (u,v).

%Reading in data;
Sec4Data=dlmread(['Sector4_',num2str(discPoints)]); %radial stretching factor.
n=size(Sec4Data,1);
RHO=Sec4Data(1:n,1:n);
F=Sec4Data(1:n,n+1:2*n);

 %Record the number of points that were discretized over.

% Construction of domain in s-th coordinates
s=linspace(0,M,n+1);
s=linspace(s(2),M,n);
th=linspace(a4,a1+2*pi,n);
[S,TH]=meshgrid(s,th);

%Construction of domain in visual optic space. This corresponds to the
%the coordinates on the hemisphere about the optical axis which the eye can
%see.

SV=-pi*R/(2*M)*S+pi*R;
THV=TH+pi;

%Calculation of the vector field of rotations in the coordinates svis and
%thvis. The variables (t,phi) denote geodesic polar coordinate around the axis of
%rotation. 

t=zeros(n,n);
phi=zeros(n,n);
vecS=zeros(n,n);
vecTH=zeros(n,n);


%This loop constructs the components of the vector field of rotations
%in the visual space of the eye.  The function InverseSTH determines the
%coordinates of (t,phi) given coordiantes svis and thvis and is necessary to
%calculate the components of the vector field. 

for i=1:n
    for j=1:n,
        [t(i,j),phi(i,j)]=InverseSTH(alpha,beta,R,SV(i,j),THV(i,j)); 
        
            vecSV(i,j)=-cos(SV(i,j)/R)*cos(beta)*sin(phi(i,j))*cos(THV(i,j)-alpha)+...
                cos(SV(i,j)/R)*cos(phi(i,j))*(sin(THV(i,j)-alpha))-...
                sin(SV(i,j)/R)*sin(phi(i,j))*sin(beta);
            
            vecSV(i,j)=R*sin(t(i,j)/R)*vecSV(i,j);
            
            vecTHV(i,j)=(cos(beta)*sin(phi(i,j))*sin(THV(i,j)-alpha)+...
                cos(phi(i,j))*cos(alpha-THV(i,j)));
        
            vecTHV(i,j)=R*sin(t(i,j)/R)*vecTHV(i,j);
            vecTHV(i,j)=1/(R*sin(SV(i,j)/R))*vecTHV(i,j);
            
    end
end

vecTHV=vecTHV.*polarity;

%Conversion of the vector fields in visual space to vector fields on the
%retina.

vecS=-2*M/(pi*R)*vecSV;
vecTH=vecTHV;

%Calcution of derivatives
ds=s(2)-s(1);
DS=DiffX(n)/ds;
RHOS=(DS*RHO')';
FS=(DS*F')';

dth=th(2)-th(1);
DTH=DiffX(n)/dth;
RHOTH=(DTH*RHO);
FTH=(DTH*F);

%Calculating the vector fields from the formula for the push forward of the
%derivative.
vecU=vecS.*(RHOS.*cos(F)-RHO.*FS.*sin(F))+vecTH.*(RHOTH.*cos(F)-RHO.*FTH.*sin(F));
vecV=vecS.*(RHOS.*sin(F)+RHO.*FS.*cos(F))+vecTH.*(RHOTH.*sin(F)+RHO.*FTH.*cos(F));

%Plotting vector field in the sector. We plot the physical retina so we can
%see where the vectors are located on the physical retina.
figure(1);
hold on;

C=abs(R.^2*sin(S./R).^2-S.^2);
surf(RHO.*cos(F),RHO.*sin(F),zeros(n,n),C);
view(0,90);
colormap(jet);
caxis([0,1]);
colorbar;
colormap(surfcolor);
shading interp;

U=RHO.*cos(F);
V=RHO.*sin(F);
quiver(U,V,vecU,vecV,vColor);
axis square;
xlabel('u');
ylabel('v');

hold off;

%Plotting vector field in (S,TH) coordiantes;
figure(2);
hold on;
vecStemp=vecS./(vecS.^2+vecTH.^2).^(1/2);
vecTHtemp=vecTH./(vecS.^2+vecTH.^2).^(1/2);

vecS=vecStemp;
vecTH=vecTHtemp;

quiver(S,TH/(2*pi),vecS,vecTH,vColor);
plot(beta*R,alpha/(2*pi),'o');
axis([0,1,0,1]);
axis square;

set(gcf,'units','normalized','outerposition',[0 0 1 1]);

hold off;

%Plotting Vector Field in Generalized S TH coordinates
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

quiver(X,Y,vecX,vecY,vColor);
plot(R*beta*cos(alpha),R*beta*sin(alpha),vColor);
axis([-1,1,-1,1]);
axis square;

xlabel('x');
ylabel('y');

hold off;

set(gcf,'units','normalized','outerposition',[0 0 1 1]);

hold off;

Xprime=X;
Yprime=Y;

vecXprime=vecX;
vecYprime=vecY;

%Plotting Vector Fields on the sphere
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

%ntemp=40;

s=linspace(0,M,ntemp);
th=linspace(a4,a1+2*pi,ntemp);
[S,TH]=meshgrid(s,th);

surf(R*sin(S/R).*cos(TH),R*sin(S/R).*sin(TH),R-R*cos(S/R),S,'EdgeColor','none','FaceAlpha','interp','AlphaDataMapping','scaled',...
    'AlphaData',S);
colormap(surfcolor);
alphamap('default')
view(104,6);

quiver3(X,Y,Z,vecX,vecY,vecZ,vColor);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

hold off;

%% Plotting vector field in standard retina coordinates
StSec4Data=dlmread(['StandardSector4_',num2str(discPoints)]);

StRHO=StSec4Data(1:n,1:n);
StF=StSec4Data(1:n,n+1:2*n);
StS=StSec4Data(1:n,2*n+1:3*n);
StTH=StSec4Data(1:n,3*n+1:4*n);

StU=StRHO.*cos(StF);
StV=StRHO.*sin(StF);

%Calcution of derivatives
ds=StS(1,2)-StS(1,1);
DS=DiffX(n)/ds;
StRHOS=(DS*StRHO')';
StFS=(DS*StF')';

dth=StTH(2,1)-StTH(1,1);
DTH=DiffX(n)/dth;
StRHOTH=(DTH*StRHO);
StFTH=(DTH*StF);

StVecU=vecS.*(StRHOS.*cos(StF)-StRHO.*StFS.*sin(StF))+vecTH.*(StRHOTH.*cos(StF)-StRHO.*StFTH.*sin(StF));
StVecV=vecS.*(StRHOS.*sin(StF)+StRHO.*StFS.*cos(StF))+vecTH.*(StRHOTH.*sin(StF)+StRHO.*StFTH.*cos(StF));

figure(5);
temp=(StVecU.^2+StVecV.^2).^(1/2);

StVecU=StVecU./temp;
StVecV=StVecV./temp;

quiver(StU,StV,StVecU,StVecV,'r');
hold off;


%Writing Vector Field Data.
dlmwrite(['Sec4VecField_',num2str(discPoints),'_',fieldType],[U,V,vecU,vecV,RHO,F,RHOS,RHOTH,FS,FTH,S,TH,vecS,vecTH,Xprime,Yprime,vecXprime,vecYprime,X,Y,Z,vecX,vecY,vecZ]);

end

