function []= VectorFieldGeneralizedForGuiAlphaCorr(alpha, beta, polarity, discPoints,fieldType,VizAngle)


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
%   9. VisAngle - Angle of the visual field.
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
%figure;
%% Input Parameters

phys=dlmread('physAlphaCorr');

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

VizAngle=deg2rad(VizAngle);

%% Computing coordinates of center of ration.

surfcolor='white';
vColor='r';

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

%% Sector 1
%This part of the code computes the vector fields in the sectors. We
%read in the numerical solution computed using the script Mapping.m to
%compute the value of rho. The vector fields are then found by computing
%the pushforward in the coordinates (u,v).

num2str(discPoints)

%Reading in data;
Sec1Data=dlmread(['Sector1AlphaCorr_',num2str(discPoints)]); %radial stretching factor.
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

SV=-VizAngle*R/(2*M)*S+pi*R;
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

vecS=-2*M/(VizAngle*R)*vecSV;
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
% figure(1);
% hold on;

C=abs(R.^2*sin(S./R).^2-S.^2);
% surf(RHO.*cos(F),RHO.*sin(F),zeros(n,n),C);
% view(0,90);
% colormap(jet);
% caxis([0,1]);
% colorbar;
% colormap(surfcolor);
% shading interp;

U=RHO.*cos(F);
V=RHO.*sin(F);

% quiver(U,V,vecU,vecV,vColor);
% axis square;
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% plot(beta*R*cos(alpha),beta*R*sin(alpha),'o');
% xlabel('u');
% ylabel('v');
% 
% drawnow expose
% hold off;

%Plotting vector field in (S,TH) coordiantes;
% figure(2);
% hold on;
vecStemp=vecS./(vecS.^2+vecTH.^2).^(1/2);
vecTHtemp=vecTH./(vecS.^2+vecTH.^2).^(1/2);

vecS=vecStemp;
vecTH=vecTHtemp;

% quiver(S,TH/(2*pi),vecS,vecTH,vColor);
% plot(beta*R,alpha/(2*pi),'o');
% axis([0,1,0,1]);
% axis square;
% 
% xlabel('s');
% ylabel('\theta');
% 
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% 
% drawnow expose
% hold off;

%Plotting Vector Field in Generalized (S,TH) coordinates
% figure(3);
% hold on;
X=S.*cos(TH);
Y=S.*sin(TH);

vecX=cos(TH).*vecS-S.*sin(TH).*vecTH;
vecY=sin(TH).*vecS+S.*cos(TH).*vecTH;

a=vecX./(vecX.^2+vecY.^2).^(1/2);
b=vecY./(vecX.^2+vecY.^2).^(1/2);

vecX=a;
vecY=b;

% quiver(X,Y,vecX,vecY,vColor);
% plot(R*beta*cos(alpha),R*beta*sin(alpha),'o');
% axis([-1,1,-1,1]);
% axis square;
% 
% xlabel('x^{\prime}');
% ylabel('y^{\prime}');
% 
% drawnow expose
% hold off;

Xprime=X;
Yprime=Y;

vecXprime=vecX;
vecYprime=vecY;

%Plotting Vector Fields on the sphere
% figure(4);
% hold on;
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

% surf(R*sin(S/R).*cos(TH),R*sin(S/R).*sin(TH),R-R*cos(S/R),S,'EdgeColor','none','FaceAlpha','interp','AlphaDataMapping','scaled',...
%     'AlphaData',S);
% colormap(surfcolor);
% alphamap('default')
% view(104,6);
% 
% quiver3(X,Y,Z,vecX,vecY,vecZ,vColor);
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% 
% axis([-R,R,-R,R,0,2*R]);
% axis square;
% 
% plot3(R*cos(alpha)*sin(beta),R*sin(alpha)*sin(beta),R-R*cos(beta),'o');
% xlabel('x');
% ylabel('y');
% zlabel('z');
% drawnow expose
% hold off;

%% Plotting vector field in standard retina coordinates

for i=1:n,
    for j=1:n,
    angle=mod(TH(i,j),2*pi);
    if angle<pi/2,
        
        %Calculation of closest grid point in S,TH coordinates.
        temp=(S(i,j)-StS1).^2+(angle-StTH1).^2;
        [ind1,ind2] =find(temp==min(min(temp)));
        a(i)=ind1(1);
        b(i)=ind2(1);
        
        USt(i,j)=StRHO1(a(i),b(i))*cos(StF1(a(i),b(i)));
        VSt(i,j)=StRHO1(a(i),b(i))*sin(StF1(a(i),b(i)));
        
        VecUSt(i,j)=vecS(i,j).*(StRHO1S(a(i),b(i)).*cos(StF1(a(i),b(i)))-StRHO1(a(i),b(i)).*StF1S(a(i),b(i)).*sin(StF1(a(i),b(i))))+vecTH(i,j).*(StRHO1TH(a(i),b(i)).*cos(StF1(a(i),b(i)))-StRHO1(a(i),b(i)).*StF1TH(a(i),b(i)).*sin(StF1(a(i),b(i))));
        VecVSt(i,j)=vecS(i,j).*(StRHO1S(a(i),b(i)).*sin(StF1(a(i),b(i)))+StRHO1(a(i),b(i)).*StF1S(a(i),b(i)).*cos(StF1(a(i),b(i))))+vecTH(i,j).*(StRHO1TH(a(i),b(i)).*sin(StF1(a(i),b(i)))+StRHO1(a(i),b(i)).*StF1TH(a(i),b(i)).*cos(StF1(a(i),b(i))));
        
    elseif angle<pi,
        %Calculation of closest grid point in S,TH coordinates.
        temp=(S(i,j)-StS2).^2+(angle-StTH2).^2;
        [ind1,ind2] =find(temp==min(min(temp)));
        a(i)=ind1(1);
        b(i)=ind2(1);
        
        USt(i,j)=StRHO2(a(i),b(i))*cos(StF2(a(i),b(i)));
        VSt(i,j)=StRHO2(a(i),b(i))*sin(StF2(a(i),b(i)));
        
        VecUSt(i,j)=vecS(i,j).*(StRHO2S(a(i),b(i)).*cos(StF2(a(i),b(i)))-StRHO2(a(i),b(i)).*StF2S(a(i),b(i)).*sin(StF2(a(i),b(i))))+vecTH(i,j).*(StRHO2TH(a(i),b(i)).*cos(StF2(a(i),b(i)))-StRHO2(a(i),b(i)).*StF2TH(a(i),b(i)).*sin(StF2(a(i),b(i))));
        VecVSt(i,j)=vecS(i,j).*(StRHO2S(a(i),b(i)).*sin(StF2(a(i),b(i)))+StRHO2(a(i),b(i)).*StF2S(a(i),b(i)).*cos(StF2(a(i),b(i))))+vecTH(i,j).*(StRHO2TH(a(i),b(i)).*sin(StF2(a(i),b(i)))+StRHO2(a(i),b(i)).*StF2TH(a(i),b(i)).*cos(StF2(a(i),b(i))));
        
        
    elseif angle<3*pi/2,
        %Calculation of closest grid point in S,TH coordinates.
        temp=(S(i,j)-StS3).^2+(angle-StTH3).^2;
        [ind1,ind2] =find(temp==min(min(temp)));
        a(i)=ind1(1);
        b(i)=ind2(1);
        
        USt(i,j)=StRHO3(a(i),b(i))*cos(StF3(a(i),b(i)));
        VSt(i,j)=StRHO3(a(i),b(i))*sin(StF3(a(i),b(i)));
        
        VecUSt(i,j)=vecS(i,j).*(StRHO3S(a(i),b(i)).*cos(StF3(a(i),b(i)))-StRHO3(a(i),b(i)).*StF3S(a(i),b(i)).*sin(StF3(a(i),b(i))))+vecTH(i,j).*(StRHO3TH(a(i),b(i)).*cos(StF3(a(i),b(i)))-StRHO3(a(i),b(i)).*StF3TH(a(i),b(i)).*sin(StF3(a(i),b(i))));
        VecVSt(i,j)=vecS(i,j).*(StRHO3S(a(i),b(i)).*sin(StF3(a(i),b(i)))+StRHO3(a(i),b(i)).*StF3S(a(i),b(i)).*cos(StF3(a(i),b(i))))+vecTH(i,j).*(StRHO3TH(a(i),b(i)).*sin(StF3(a(i),b(i)))+StRHO3(a(i),b(i)).*StF3TH(a(i),b(i)).*cos(StF3(a(i),b(i))));
        
        
    else
        %Calculation of closest grid point in S,TH coordinates.
        temp=(S(i,j)-StS4).^2+(angle-StTH4).^2;
        [ind1,ind2] =find(temp==min(min(temp)));
        a(i)=ind1(1);
        b(i)=ind2(1);
        
        USt(i,j)=StRHO4(a(i),b(i))*cos(StF4(a(i),b(i)));
        VSt(i,j)=StRHO4(a(i),b(i))*sin(StF4(a(i),b(i)));
        
        VecUSt(i,j)=vecS(i,j).*(StRHO4S(a(i),b(i)).*cos(StF4(a(i),b(i)))-StRHO4(a(i),b(i)).*StF4S(a(i),b(i)).*sin(StF4(a(i),b(i))))+vecTH(i,j).*(StRHO4TH(a(i),b(i)).*cos(StF4(a(i),b(i)))-StRHO4(a(i),b(i)).*StF4TH(a(i),b(i)).*sin(StF4(a(i),b(i))));
        VecVSt(i,j)=vecS(i,j).*(StRHO4S(a(i),b(i)).*sin(StF4(a(i),b(i)))+StRHO4(a(i),b(i)).*StF4S(a(i),b(i)).*cos(StF4(a(i),b(i))))+vecTH(i,j).*(StRHO4TH(a(i),b(i)).*sin(StF4(a(i),b(i)))+StRHO4(a(i),b(i)).*StF4TH(a(i),b(i)).*cos(StF4(a(i),b(i))));
        
    end
    end
end

% figure(5);
% hold on;
temp=(VecUSt.^2+VecVSt.^2).^(1/2);
VecUSt=VecUSt./temp;
VecVSt=VecVSt./temp;

% scale=.4;
% quiver(USt,VSt,VecUSt,VecVSt,scale);
% hold off;

dlmwrite(['Sec1VecFieldAlphaCorr_',num2str(discPoints),'_',fieldType],[U,V,vecU,vecV,RHO,F,RHOS,RHOTH,FS,FTH,S,TH,vecS,vecTH,Xprime,Yprime,vecXprime,vecYprime,X,Y,Z,vecX,vecY,vecZ]);


% Sector 2
%This part of the code computes the vector fields in the sectors. We
%read in the numerical solution computed using the script Mapping.m to
%compute the value of rho. The vector fields are then found by computing
%the pushforward in the coordinates (u,v).

%Reading in data;
Sec2Data=dlmread(['Sector2AlphaCorr_',num2str(discPoints)]); %radial stretching factor.
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

SV=-VizAngle*R/(2*M)*S+pi*R;
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

vecS=-2*M/(VizAngle*R)*vecSV;
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
% figure(1);
% hold on;

C=abs(R.^2*sin(S./R).^2-S.^2);
% surf(RHO.*cos(F),RHO.*sin(F),zeros(n,n),C);
% view(0,90);
% colormap(jet);
% caxis([0,1]);
% colorbar;
% colormap(surfcolor);
% shading interp;

U=RHO.*cos(F);
V=RHO.*sin(F);
% quiver(U,V,vecU,vecV,vColor);
% axis square;
% axis ([-1,1,-1,1]);

%Plotting vector field in (S,TH) coordiantes;
% figure(2);
% hold on;
vecStemp=vecS./(vecS.^2+vecTH.^2).^(1/2);
vecTHtemp=vecTH./(vecS.^2+vecTH.^2).^(1/2);

vecS=vecStemp;
vecTH=vecTHtemp;

% quiver(S,TH/(2*pi),vecS,vecTH,vColor);
% plot(beta*R,alpha/(2*pi),'o');
% axis([0,1,0,1]);
% axis square;
% 
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% 
% drawnow expose
% hold off;

%Plotting Vector Field in Generalized S TH coordinates
% figure(3);
% hold on;
X=S.*cos(TH);
Y=S.*sin(TH);

vecX=cos(TH).*vecS-S.*sin(TH).*vecTH;
vecY=sin(TH).*vecS+S.*cos(TH).*vecTH;

a=vecX./(vecX.^2+vecY.^2).^(1/2);
b=vecY./(vecX.^2+vecY.^2).^(1/2);

vecX=a;
vecY=b;

% quiver(X,Y,vecX,vecY,vColor);
% plot(R*beta*cos(alpha),R*beta*sin(alpha),'o');
% axis([-1,1,-1,1]);
% axis square;
% 
% xlabel('x');
% ylabel('y');
% 
% drawnow expose
% hold off;
% 
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% 
% hold off;

Xprime=X;
Yprime=Y;

vecXprime=vecX;
vecYprime=vecY;
%Plotting Vector Fields on the sphere
% figure(4);
% hold on;
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

% surf(R*sin(S/R).*cos(TH),R*sin(S/R).*sin(TH),R-R*cos(S/R),S,'EdgeColor','none','FaceAlpha','interp','AlphaDataMapping','scaled',...
%     'AlphaData',S);
% colormap(surfcolor);
% alphamap('default')
% view(104,6);
% 
% quiver3(X,Y,Z,vecX,vecY,vecZ,vColor);
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% 
% drawnow expose
% hold off;

%% Plotting vector field in standard retina coordinates

for i=1:n,
    for j=1:n,
    angle=mod(TH(i,j),2*pi);
    if angle<pi/2,
        
        %Calculation of closest grid point in S,TH coordinates.
        temp=(S(i,j)-StS1).^2+(angle-StTH1).^2;
        [ind1,ind2] =find(temp==min(min(temp)));
        a(i)=ind1(1);
        b(i)=ind2(1);
        
        USt(i,j)=StRHO1(a(i),b(i))*cos(StF1(a(i),b(i)));
        VSt(i,j)=StRHO1(a(i),b(i))*sin(StF1(a(i),b(i)));
        
        VecUSt(i,j)=vecS(i,j).*(StRHO1S(a(i),b(i)).*cos(StF1(a(i),b(i)))-StRHO1(a(i),b(i)).*StF1S(a(i),b(i)).*sin(StF1(a(i),b(i))))+vecTH(i,j).*(StRHO1TH(a(i),b(i)).*cos(StF1(a(i),b(i)))-StRHO1(a(i),b(i)).*StF1TH(a(i),b(i)).*sin(StF1(a(i),b(i))));
        VecVSt(i,j)=vecS(i,j).*(StRHO1S(a(i),b(i)).*sin(StF1(a(i),b(i)))+StRHO1(a(i),b(i)).*StF1S(a(i),b(i)).*cos(StF1(a(i),b(i))))+vecTH(i,j).*(StRHO1TH(a(i),b(i)).*sin(StF1(a(i),b(i)))+StRHO1(a(i),b(i)).*StF1TH(a(i),b(i)).*cos(StF1(a(i),b(i))));
        
    elseif angle<pi,
        %Calculation of closest grid point in S,TH coordinates.
        temp=(S(i,j)-StS2).^2+(angle-StTH2).^2;
        [ind1,ind2] =find(temp==min(min(temp)));
        a(i)=ind1(1);
        b(i)=ind2(1);
        
        USt(i,j)=StRHO2(a(i),b(i))*cos(StF2(a(i),b(i)));
        VSt(i,j)=StRHO2(a(i),b(i))*sin(StF2(a(i),b(i)));
        
        VecUSt(i,j)=vecS(i,j).*(StRHO2S(a(i),b(i)).*cos(StF2(a(i),b(i)))-StRHO2(a(i),b(i)).*StF2S(a(i),b(i)).*sin(StF2(a(i),b(i))))+vecTH(i,j).*(StRHO2TH(a(i),b(i)).*cos(StF2(a(i),b(i)))-StRHO2(a(i),b(i)).*StF2TH(a(i),b(i)).*sin(StF2(a(i),b(i))));
        VecVSt(i,j)=vecS(i,j).*(StRHO2S(a(i),b(i)).*sin(StF2(a(i),b(i)))+StRHO2(a(i),b(i)).*StF2S(a(i),b(i)).*cos(StF2(a(i),b(i))))+vecTH(i,j).*(StRHO2TH(a(i),b(i)).*sin(StF2(a(i),b(i)))+StRHO2(a(i),b(i)).*StF2TH(a(i),b(i)).*cos(StF2(a(i),b(i))));
        
        
    elseif angle<3*pi/2,
        %Calculation of closest grid point in S,TH coordinates.
        temp=(S(i,j)-StS3).^2+(angle-StTH3).^2;
        [ind1,ind2] =find(temp==min(min(temp)));
        a(i)=ind1(1);
        b(i)=ind2(1);
        
        USt(i,j)=StRHO3(a(i),b(i))*cos(StF3(a(i),b(i)));
        VSt(i,j)=StRHO3(a(i),b(i))*sin(StF3(a(i),b(i)));
        
        VecUSt(i,j)=vecS(i,j).*(StRHO3S(a(i),b(i)).*cos(StF3(a(i),b(i)))-StRHO3(a(i),b(i)).*StF3S(a(i),b(i)).*sin(StF3(a(i),b(i))))+vecTH(i,j).*(StRHO3TH(a(i),b(i)).*cos(StF3(a(i),b(i)))-StRHO3(a(i),b(i)).*StF3TH(a(i),b(i)).*sin(StF3(a(i),b(i))));
        VecVSt(i,j)=vecS(i,j).*(StRHO3S(a(i),b(i)).*sin(StF3(a(i),b(i)))+StRHO3(a(i),b(i)).*StF3S(a(i),b(i)).*cos(StF3(a(i),b(i))))+vecTH(i,j).*(StRHO3TH(a(i),b(i)).*sin(StF3(a(i),b(i)))+StRHO3(a(i),b(i)).*StF3TH(a(i),b(i)).*cos(StF3(a(i),b(i))));
        
        
    else
        %Calculation of closest grid point in S,TH coordinates.
        temp=(S(i,j)-StS4).^2+(angle-StTH4).^2;
        [ind1,ind2] =find(temp==min(min(temp)));
        a(i)=ind1(1);
        b(i)=ind2(1);
        
        USt(i,j)=StRHO4(a(i),b(i))*cos(StF4(a(i),b(i)));
        VSt(i,j)=StRHO4(a(i),b(i))*sin(StF4(a(i),b(i)));
        
        VecUSt(i,j)=vecS(i,j).*(StRHO4S(a(i),b(i)).*cos(StF4(a(i),b(i)))-StRHO4(a(i),b(i)).*StF4S(a(i),b(i)).*sin(StF4(a(i),b(i))))+vecTH(i,j).*(StRHO4TH(a(i),b(i)).*cos(StF4(a(i),b(i)))-StRHO4(a(i),b(i)).*StF4TH(a(i),b(i)).*sin(StF4(a(i),b(i))));
        VecVSt(i,j)=vecS(i,j).*(StRHO4S(a(i),b(i)).*sin(StF4(a(i),b(i)))+StRHO4(a(i),b(i)).*StF4S(a(i),b(i)).*cos(StF4(a(i),b(i))))+vecTH(i,j).*(StRHO4TH(a(i),b(i)).*sin(StF4(a(i),b(i)))+StRHO4(a(i),b(i)).*StF4TH(a(i),b(i)).*cos(StF4(a(i),b(i))));
        
    end
    end
end

% figure(5);
% hold on;
temp=(VecUSt.^2+VecVSt.^2).^(1/2);
VecUSt=VecUSt./temp;
VecVSt=VecVSt./temp;

% scale=.4;
% quiver(USt,VSt,VecUSt,VecVSt,scale);
% hold off;


%Writing Vector Field Data.
dlmwrite(['Sec2VecFieldAlphaCorr_',num2str(discPoints),'_',fieldType],[U,V,vecU,vecV,RHO,F,RHOS,RHOTH,FS,FTH,S,TH,vecS,vecTH,Xprime,Yprime,vecXprime,vecYprime,X,Y,Z,vecX,vecY,vecZ]);


% Sector 3
%This part of the code computes the vector fields in the sectors. We
%read in the numerical solution computed using the script Mapping.m to
%compute the value of rho. The vector fields are then found by computing
%the pushforward in the coordinates (u,v).

%Reading in data;
Sec3Data=dlmread(['Sector3AlphaCorr_',num2str(discPoints)]); %radial stretching factor.
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

SV=-VizAngle*R/(2*M)*S+pi*R;
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

vecS=-2*M/(VizAngle*R)*vecSV;
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
% figure(1);
% hold on;

C=abs(R.^2*sin(S./R).^2-S.^2);
% surf(RHO.*cos(F),RHO.*sin(F),zeros(n,n),C);
% view(0,90);
% colormap(jet);
% caxis([0,1]);
% colorbar;
% colormap('jet');
% shading interp;

U=RHO.*cos(F);
V=RHO.*sin(F);
% quiver(U,V,vecU,vecV,vColor);
% axis square;
% 
% xlabel('u');
% ylabel('v');
% drawnow expose
% hold off;

%Plotting vector field in (S,TH) coordiantes;
% figure(2);
% hold on;
vecStemp=vecS./(vecS.^2+vecTH.^2).^(1/2);
vecTHtemp=vecTH./(vecS.^2+vecTH.^2).^(1/2);

vecS=vecStemp;
vecTH=vecTHtemp;

% quiver(S,TH/(2*pi),vecS,vecTH,vColor);
% plot(beta*R,alpha/(2*pi),'o');
% axis([0,1,0,1]);
% axis square;
% 
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% 
% drawnow expose
% hold off;

%Plotting Vector Field in Generalized S TH coordinates
% figure(3);
% hold on;
X=S.*cos(TH);
Y=S.*sin(TH);

vecX=cos(TH).*vecS-S.*sin(TH).*vecTH;
vecY=sin(TH).*vecS+S.*cos(TH).*vecTH;

a=vecX./(vecX.^2+vecY.^2).^(1/2);
b=vecY./(vecX.^2+vecY.^2).^(1/2);

vecX=a;
vecY=b;

% quiver(X,Y,vecX,vecY,vColor);
% plot(R*beta*cos(alpha),R*beta*sin(alpha),'o');
% axis([-1,1,-1,1]);
% axis square;
% 
% xlabel('x');
% ylabel('y');
% 
% drawnow expose
% hold off;
% 
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% 
% hold off;

Xprime=X;
Yprime=Y;

vecXprime=vecX;
vecYprime=vecY;

%Plotting Vector Fields on the sphere
% figure(4);
% hold on;
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

% surf(R*sin(S/R).*cos(TH),R*sin(S/R).*sin(TH),R-R*cos(S/R),S,'EdgeColor','none','FaceAlpha','interp','AlphaDataMapping','scaled',...
%     'AlphaData',S);
% colormap(surfcolor);
% alphamap('default')
% view(104,6);
% 
% quiver3(X,Y,Z,vecX,vecY,vecZ,vColor);
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% 
% drawnow expose
% hold off;

%% Plotting vector field in standard retina coordinates

for i=1:n,
    for j=1:n,
    angle=mod(TH(i,j),2*pi);
    if angle<pi/2,
        
        %Calculation of closest grid point in S,TH coordinates.
        temp=(S(i,j)-StS1).^2+(angle-StTH1).^2;
        [ind1,ind2] =find(temp==min(min(temp)));
        a(i)=ind1(1);
        b(i)=ind2(1);
        
        USt(i,j)=StRHO1(a(i),b(i))*cos(StF1(a(i),b(i)));
        VSt(i,j)=StRHO1(a(i),b(i))*sin(StF1(a(i),b(i)));
        
        VecUSt(i,j)=vecS(i,j).*(StRHO1S(a(i),b(i)).*cos(StF1(a(i),b(i)))-StRHO1(a(i),b(i)).*StF1S(a(i),b(i)).*sin(StF1(a(i),b(i))))+vecTH(i,j).*(StRHO1TH(a(i),b(i)).*cos(StF1(a(i),b(i)))-StRHO1(a(i),b(i)).*StF1TH(a(i),b(i)).*sin(StF1(a(i),b(i))));
        VecVSt(i,j)=vecS(i,j).*(StRHO1S(a(i),b(i)).*sin(StF1(a(i),b(i)))+StRHO1(a(i),b(i)).*StF1S(a(i),b(i)).*cos(StF1(a(i),b(i))))+vecTH(i,j).*(StRHO1TH(a(i),b(i)).*sin(StF1(a(i),b(i)))+StRHO1(a(i),b(i)).*StF1TH(a(i),b(i)).*cos(StF1(a(i),b(i))));
        
    elseif angle<pi,
        %Calculation of closest grid point in S,TH coordinates.
        temp=(S(i,j)-StS2).^2+(angle-StTH2).^2;
        [ind1,ind2] =find(temp==min(min(temp)));
        a(i)=ind1(1);
        b(i)=ind2(1);
        
        USt(i,j)=StRHO2(a(i),b(i))*cos(StF2(a(i),b(i)));
        VSt(i,j)=StRHO2(a(i),b(i))*sin(StF2(a(i),b(i)));
        
        VecUSt(i,j)=vecS(i,j).*(StRHO2S(a(i),b(i)).*cos(StF2(a(i),b(i)))-StRHO2(a(i),b(i)).*StF2S(a(i),b(i)).*sin(StF2(a(i),b(i))))+vecTH(i,j).*(StRHO2TH(a(i),b(i)).*cos(StF2(a(i),b(i)))-StRHO2(a(i),b(i)).*StF2TH(a(i),b(i)).*sin(StF2(a(i),b(i))));
        VecVSt(i,j)=vecS(i,j).*(StRHO2S(a(i),b(i)).*sin(StF2(a(i),b(i)))+StRHO2(a(i),b(i)).*StF2S(a(i),b(i)).*cos(StF2(a(i),b(i))))+vecTH(i,j).*(StRHO2TH(a(i),b(i)).*sin(StF2(a(i),b(i)))+StRHO2(a(i),b(i)).*StF2TH(a(i),b(i)).*cos(StF2(a(i),b(i))));
        
        
    elseif angle<3*pi/2,
        %Calculation of closest grid point in S,TH coordinates.
        temp=(S(i,j)-StS3).^2+(angle-StTH3).^2;
        [ind1,ind2] =find(temp==min(min(temp)));
        a(i)=ind1(1);
        b(i)=ind2(1);
        
        USt(i,j)=StRHO3(a(i),b(i))*cos(StF3(a(i),b(i)));
        VSt(i,j)=StRHO3(a(i),b(i))*sin(StF3(a(i),b(i)));
        
        VecUSt(i,j)=vecS(i,j).*(StRHO3S(a(i),b(i)).*cos(StF3(a(i),b(i)))-StRHO3(a(i),b(i)).*StF3S(a(i),b(i)).*sin(StF3(a(i),b(i))))+vecTH(i,j).*(StRHO3TH(a(i),b(i)).*cos(StF3(a(i),b(i)))-StRHO3(a(i),b(i)).*StF3TH(a(i),b(i)).*sin(StF3(a(i),b(i))));
        VecVSt(i,j)=vecS(i,j).*(StRHO3S(a(i),b(i)).*sin(StF3(a(i),b(i)))+StRHO3(a(i),b(i)).*StF3S(a(i),b(i)).*cos(StF3(a(i),b(i))))+vecTH(i,j).*(StRHO3TH(a(i),b(i)).*sin(StF3(a(i),b(i)))+StRHO3(a(i),b(i)).*StF3TH(a(i),b(i)).*cos(StF3(a(i),b(i))));
        
        
    else
        %Calculation of closest grid point in S,TH coordinates.
        temp=(S(i,j)-StS4).^2+(angle-StTH4).^2;
        [ind1,ind2] =find(temp==min(min(temp)));
        a(i)=ind1(1);
        b(i)=ind2(1);
        
        USt(i,j)=StRHO4(a(i),b(i))*cos(StF4(a(i),b(i)));
        VSt(i,j)=StRHO4(a(i),b(i))*sin(StF4(a(i),b(i)));
        
        VecUSt(i,j)=vecS(i,j).*(StRHO4S(a(i),b(i)).*cos(StF4(a(i),b(i)))-StRHO4(a(i),b(i)).*StF4S(a(i),b(i)).*sin(StF4(a(i),b(i))))+vecTH(i,j).*(StRHO4TH(a(i),b(i)).*cos(StF4(a(i),b(i)))-StRHO4(a(i),b(i)).*StF4TH(a(i),b(i)).*sin(StF4(a(i),b(i))));
        VecVSt(i,j)=vecS(i,j).*(StRHO4S(a(i),b(i)).*sin(StF4(a(i),b(i)))+StRHO4(a(i),b(i)).*StF4S(a(i),b(i)).*cos(StF4(a(i),b(i))))+vecTH(i,j).*(StRHO4TH(a(i),b(i)).*sin(StF4(a(i),b(i)))+StRHO4(a(i),b(i)).*StF4TH(a(i),b(i)).*cos(StF4(a(i),b(i))));
        
    end
    end
end

% figure(5);
% hold on;
temp=(VecUSt.^2+VecVSt.^2).^(1/2);
VecUSt=VecUSt./temp;
VecVSt=VecVSt./temp;

% scale=.4;
% quiver(USt,VSt,VecUSt,VecVSt,scale);
% hold off;

%Writing Vector Field Data.
dlmwrite(['Sec3VecFieldAlphaCorr_',num2str(discPoints),'_',fieldType],[U,V,vecU,vecV,RHO,F,RHOS,RHOTH,FS,FTH,S,TH,vecS,vecTH,Xprime,Yprime,vecXprime,vecYprime,X,Y,Z,vecX,vecY,vecZ]);


% Sector 4
%This part of the code computes the vector fields in the sectors. We
%read in the numerical solution computed using the script Mapping.m to
%compute the value of rho. The vector fields are then found by computing
%the pushforward in the coordinates (u,v).

%Reading in data;
Sec4Data=dlmread(['Sector4AlphaCorr_',num2str(discPoints)]); %radial stretching factor.
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

SV=-VizAngle*R/(2*M)*S+pi*R;
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

vecS=-2*M/(VizAngle*R)*vecSV;
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
% figure(1);
% hold on;

C=abs(R.^2*sin(S./R).^2-S.^2);
% surf(RHO.*cos(F),RHO.*sin(F),zeros(n,n),C);
% view(0,90);
% colormap(jet);
% caxis([0,1]);
% colorbar;
% colormap(surfcolor);
% shading interp;

U=RHO.*cos(F);
V=RHO.*sin(F);
% quiver(U,V,vecU,vecV,vColor);
% axis square;
% xlabel('u');
% ylabel('v');
% 
% drawnow expose
% hold off;

%Plotting vector field in (S,TH) coordiantes;
% figure(2);
% hold on;
vecStemp=vecS./(vecS.^2+vecTH.^2).^(1/2);
vecTHtemp=vecTH./(vecS.^2+vecTH.^2).^(1/2);

vecS=vecStemp;
vecTH=vecTHtemp;

% quiver(S,TH/(2*pi),vecS,vecTH,vColor);
% plot(beta*R,alpha/(2*pi),'o');
% axis([0,1,0,1]);
% axis square;
% 
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% 
% drawnow expose
% hold off;

%Plotting Vector Field in Generalized S TH coordinates
% figure(3);
% hold on;
X=S.*cos(TH);
Y=S.*sin(TH);

vecX=cos(TH).*vecS-S.*sin(TH).*vecTH;
vecY=sin(TH).*vecS+S.*cos(TH).*vecTH;

a=vecX./(vecX.^2+vecY.^2).^(1/2);
b=vecY./(vecX.^2+vecY.^2).^(1/2);

vecX=a;
vecY=b;

% quiver(X,Y,vecX,vecY,vColor);
% plot(R*beta*cos(alpha),R*beta*sin(alpha),vColor);
% axis([-1,1,-1,1]);
% axis square;
% 
% xlabel('x');
% ylabel('y');
% 
% drawnow expose
% hold off;
% 
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% 
% hold off;

Xprime=X;
Yprime=Y;

vecXprime=vecX;
vecYprime=vecY;

%Plotting Vector Fields on the sphere
% figure(4);
% hold on;
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

% surf(R*sin(S/R).*cos(TH),R*sin(S/R).*sin(TH),R-R*cos(S/R),S,'EdgeColor','none','FaceAlpha','interp','AlphaDataMapping','scaled',...
%     'AlphaData',S);
% colormap(surfcolor);
% alphamap('default')
% view(104,6);
% 
% quiver3(X,Y,Z,vecX,vecY,vecZ,vColor);
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% 
% drawnow expose
% hold off;

%% Plotting vector field in standard retina coordinates

for i=1:n,
    for j=1:n,
    angle=mod(TH(i,j),2*pi);
    if angle<pi/2,
        
        %Calculation of closest grid point in S,TH coordinates.
        temp=(S(i,j)-StS1).^2+(angle-StTH1).^2;
        [ind1,ind2] =find(temp==min(min(temp)));
        a(i)=ind1(1);
        b(i)=ind2(1);
        
        USt(i,j)=StRHO1(a(i),b(i))*cos(StF1(a(i),b(i)));
        VSt(i,j)=StRHO1(a(i),b(i))*sin(StF1(a(i),b(i)));
        
        VecUSt(i,j)=vecS(i,j).*(StRHO1S(a(i),b(i)).*cos(StF1(a(i),b(i)))-StRHO1(a(i),b(i)).*StF1S(a(i),b(i)).*sin(StF1(a(i),b(i))))+vecTH(i,j).*(StRHO1TH(a(i),b(i)).*cos(StF1(a(i),b(i)))-StRHO1(a(i),b(i)).*StF1TH(a(i),b(i)).*sin(StF1(a(i),b(i))));
        VecVSt(i,j)=vecS(i,j).*(StRHO1S(a(i),b(i)).*sin(StF1(a(i),b(i)))+StRHO1(a(i),b(i)).*StF1S(a(i),b(i)).*cos(StF1(a(i),b(i))))+vecTH(i,j).*(StRHO1TH(a(i),b(i)).*sin(StF1(a(i),b(i)))+StRHO1(a(i),b(i)).*StF1TH(a(i),b(i)).*cos(StF1(a(i),b(i))));
        
    elseif angle<pi,
        %Calculation of closest grid point in S,TH coordinates.
        temp=(S(i,j)-StS2).^2+(angle-StTH2).^2;
        [ind1,ind2] =find(temp==min(min(temp)));
        a(i)=ind1(1);
        b(i)=ind2(1);
        
        USt(i,j)=StRHO2(a(i),b(i))*cos(StF2(a(i),b(i)));
        VSt(i,j)=StRHO2(a(i),b(i))*sin(StF2(a(i),b(i)));
        
        VecUSt(i,j)=vecS(i,j).*(StRHO2S(a(i),b(i)).*cos(StF2(a(i),b(i)))-StRHO2(a(i),b(i)).*StF2S(a(i),b(i)).*sin(StF2(a(i),b(i))))+vecTH(i,j).*(StRHO2TH(a(i),b(i)).*cos(StF2(a(i),b(i)))-StRHO2(a(i),b(i)).*StF2TH(a(i),b(i)).*sin(StF2(a(i),b(i))));
        VecVSt(i,j)=vecS(i,j).*(StRHO2S(a(i),b(i)).*sin(StF2(a(i),b(i)))+StRHO2(a(i),b(i)).*StF2S(a(i),b(i)).*cos(StF2(a(i),b(i))))+vecTH(i,j).*(StRHO2TH(a(i),b(i)).*sin(StF2(a(i),b(i)))+StRHO2(a(i),b(i)).*StF2TH(a(i),b(i)).*cos(StF2(a(i),b(i))));
        
        
    elseif angle<3*pi/2,
        %Calculation of closest grid point in S,TH coordinates.
        temp=(S(i,j)-StS3).^2+(angle-StTH3).^2;
        [ind1,ind2] =find(temp==min(min(temp)));
        a(i)=ind1(1);
        b(i)=ind2(1);
        
        USt(i,j)=StRHO3(a(i),b(i))*cos(StF3(a(i),b(i)));
        VSt(i,j)=StRHO3(a(i),b(i))*sin(StF3(a(i),b(i)));
        
        VecUSt(i,j)=vecS(i,j).*(StRHO3S(a(i),b(i)).*cos(StF3(a(i),b(i)))-StRHO3(a(i),b(i)).*StF3S(a(i),b(i)).*sin(StF3(a(i),b(i))))+vecTH(i,j).*(StRHO3TH(a(i),b(i)).*cos(StF3(a(i),b(i)))-StRHO3(a(i),b(i)).*StF3TH(a(i),b(i)).*sin(StF3(a(i),b(i))));
        VecVSt(i,j)=vecS(i,j).*(StRHO3S(a(i),b(i)).*sin(StF3(a(i),b(i)))+StRHO3(a(i),b(i)).*StF3S(a(i),b(i)).*cos(StF3(a(i),b(i))))+vecTH(i,j).*(StRHO3TH(a(i),b(i)).*sin(StF3(a(i),b(i)))+StRHO3(a(i),b(i)).*StF3TH(a(i),b(i)).*cos(StF3(a(i),b(i))));
        
        
    else
        %Calculation of closest grid point in S,TH coordinates.
        temp=(S(i,j)-StS4).^2+(angle-StTH4).^2;
        [ind1,ind2] =find(temp==min(min(temp)));
        a(i)=ind1(1);
        b(i)=ind2(1);
        
        USt(i,j)=StRHO4(a(i),b(i))*cos(StF4(a(i),b(i)));
        VSt(i,j)=StRHO4(a(i),b(i))*sin(StF4(a(i),b(i)));
        
        VecUSt(i,j)=vecS(i,j).*(StRHO4S(a(i),b(i)).*cos(StF4(a(i),b(i)))-StRHO4(a(i),b(i)).*StF4S(a(i),b(i)).*sin(StF4(a(i),b(i))))+vecTH(i,j).*(StRHO4TH(a(i),b(i)).*cos(StF4(a(i),b(i)))-StRHO4(a(i),b(i)).*StF4TH(a(i),b(i)).*sin(StF4(a(i),b(i))));
        VecVSt(i,j)=vecS(i,j).*(StRHO4S(a(i),b(i)).*sin(StF4(a(i),b(i)))+StRHO4(a(i),b(i)).*StF4S(a(i),b(i)).*cos(StF4(a(i),b(i))))+vecTH(i,j).*(StRHO4TH(a(i),b(i)).*sin(StF4(a(i),b(i)))+StRHO4(a(i),b(i)).*StF4TH(a(i),b(i)).*cos(StF4(a(i),b(i))));
        
    end
    end
end

% figure(5);
% hold on;
temp=(VecUSt.^2+VecVSt.^2).^(1/2);
VecUSt=VecUSt./temp;
VecVSt=VecVSt./temp;

% scale=.4;
% quiver(USt,VSt,VecUSt,VecVSt,scale);
% hold off;


%Writing Vector Field Data.
dlmwrite(['Sec4VecFieldAlphaCorr_',num2str(discPoints),'_',fieldType],[U,V,vecU,vecV,RHO,F,RHOS,RHOTH,FS,FTH,S,TH,vecS,vecTH,Xprime,Yprime,vecXprime,vecYprime,X,Y,Z,vecX,vecY,vecZ]);


%% Standard Vector Field No Gaps

% Read in data
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

% Sector 1

% Construction of domain in s-th coordinates
s=linspace(0,M,n+1);
s=linspace(s(2),M,n);
th=linspace(0,pi/2,n);
[S,TH]=meshgrid(s,th);

%Construction of domain in visual optic space. This corresponds to the
%the coordinates on the hemisphere about the optical axis which the eye can
%see.

SV=-VizAngle*R/(2*M)*S+pi*R;
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

vecS=-2*M/(VizAngle*R)*vecSV;
vecTH=vecTHV;

USt=StRHO1.*cos(StF1);
        VSt=StRHO1.*sin(StF1);
        
        VecUSt=vecS.*(StRHO1S.*cos(StF1)-StRHO1.*StF1S.*sin(StF1))+vecTH.*(StRHO1TH.*cos(StF1)-StRHO1.*StF1TH.*sin(StF1));
        VecVSt=vecS.*(StRHO1S.*sin(StF1)+StRHO1.*StF1S.*cos(StF1))+vecTH.*(StRHO1TH.*sin(StF1)+StRHO1.*StF1TH.*cos(StF1));
        
        figure(6);
hold on;
temp=(VecUSt.^2+VecVSt.^2).^(1/2);
VecUSt=VecUSt./temp;
VecVSt=VecVSt./temp;

scale=.4;
quiver(USt,VSt,VecUSt,VecVSt,scale);
hold off;

% Sector 2

% Construction of domain in s-th coordinates
s=linspace(0,M,n+1);
s=linspace(s(2),M,n);
th=linspace(pi/2,pi,n);
[S,TH]=meshgrid(s,th);

%Construction of domain in visual optic space. This corresponds to the
%the coordinates on the hemisphere about the optical axis which the eye can
%see.

SV=-VizAngle*R/(2*M)*S+pi*R;
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

vecS=-2*M/(VizAngle*R)*vecSV;
vecTH=vecTHV;

USt=StRHO2.*cos(StF2);
        VSt=StRHO2.*sin(StF2);
        
        VecUSt=vecS.*(StRHO2S.*cos(StF2)-StRHO2.*StF2S.*sin(StF2))+vecTH.*(StRHO2TH.*cos(StF2)-StRHO2.*StF2TH.*sin(StF2));
        VecVSt=vecS.*(StRHO2S.*sin(StF2)+StRHO2.*StF2S.*cos(StF2))+vecTH.*(StRHO2TH.*sin(StF2)+StRHO2.*StF2TH.*cos(StF2));
        
        figure(6);
hold on;
temp=(VecUSt.^2+VecVSt.^2).^(1/2);
VecUSt=VecUSt./temp;
VecVSt=VecVSt./temp;

scale=.4;
quiver(USt,VSt,VecUSt,VecVSt,scale);
hold off;

% Sector 3

% Construction of domain in s-th coordinates
s=linspace(0,M,n+1);
s=linspace(s(2),M,n);
th=linspace(pi,3*pi/2,n);
[S,TH]=meshgrid(s,th);

%Construction of domain in visual optic space. This corresponds to the
%the coordinates on the hemisphere about the optical axis which the eye can
%see.

SV=-VizAngle*R/(2*M)*S+pi*R;
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

vecS=-2*M/(VizAngle*R)*vecSV;
vecTH=vecTHV;

USt=StRHO3.*cos(StF3);
        VSt=StRHO3.*sin(StF3);
        
        VecUSt=vecS.*(StRHO3S.*cos(StF3)-StRHO3.*StF3S.*sin(StF3))+vecTH.*(StRHO3TH.*cos(StF3)-StRHO3.*StF3TH.*sin(StF3));
        VecVSt=vecS.*(StRHO3S.*sin(StF3)+StRHO3.*StF3S.*cos(StF3))+vecTH.*(StRHO3TH.*sin(StF3)+StRHO3.*StF3TH.*cos(StF3));
        
        figure(6);
hold on;
temp=(VecUSt.^2+VecVSt.^2).^(1/2);
VecUSt=VecUSt./temp;
VecVSt=VecVSt./temp;

scale=.4;
quiver(USt,VSt,VecUSt,VecVSt,scale);
hold off;

% Sector 4

% Construction of domain in s-th coordinates
s=linspace(0,M,n+1);
s=linspace(s(2),M,n);
th=linspace(3*pi/2,2*pi,n);
[S,TH]=meshgrid(s,th);

%Construction of domain in visual optic space. This corresponds to the
%the coordinates on the hemisphere about the optical axis which the eye can
%see.

SV=-VizAngle*R/(2*M)*S+pi*R;
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

vecS=-2*M/(VizAngle*R)*vecSV;
vecTH=vecTHV;

USt=StRHO4.*cos(StF4);
        VSt=StRHO4.*sin(StF4);
        
        VecUSt=vecS.*(StRHO4S.*cos(StF4)-StRHO4.*StF4S.*sin(StF4))+vecTH.*(StRHO4TH.*cos(StF4)-StRHO4.*StF4TH.*sin(StF4));
        VecVSt=vecS.*(StRHO4S.*sin(StF4)+StRHO4.*StF4S.*cos(StF4))+vecTH.*(StRHO4TH.*sin(StF4)+StRHO4.*StF4TH.*cos(StF4));
        
        figure(6);
hold on;
temp=(VecUSt.^2+VecVSt.^2).^(1/2);
VecUSt=VecUSt./temp;
VecVSt=VecVSt./temp;

scale=.4;
quiver(USt,VSt,VecUSt,VecVSt,scale);
hold off;
end

