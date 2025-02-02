%% VectorField Script

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

clear all;
figure;
%% Input Parameters

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


%% Computing coordinates of center of ration.
%Currently, we set this each time we run the script by modifiying the
%vector normal. The code can easily be modified to accept a normal as
%input.

%%%%%%%%%%%%%%%%%%%     Nathan's data     %%%%%%%%%%%%%%%%%%%%
%normal=[0.513,0.238,-0.822];        %ASC
%normal=[-0.513,-0.238,0.822];        %ASC flipped
%normal=[0.696,-0.483,0.526];        %PSC
%normal=[-0.696,0.483,-0.526];        %flipped PSC
%normal=[-0.407,-0.778,-0.470];       %LSC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% our made up data %%%%%%%%%%%%%%%%%%%%%
normal=[-0.7868,-0.5919,-0.1612];     %ASC
%normal=[0.7290,-0.6736,0.0958];        %PSC
%normal=[-0.2547,0.0903,0.9585];       %LSC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%normal=[1/sqrt(3),1/sqrt(3),1/sqrt(3)];


[alpha,beta]=NormToVecField(normal);

surfcolor='white';
vColor='r';

%%%%%%%%%%%%%%%%%%%% fixed directions %%%%%%%%%%%%%%%%%%%%%%%%%
% error=0;
% fixedLSC=-1.7+error;
% 
% %alpha=fixedLSC;                  %fixed horizontal
% alpha=deg2rad(120)+fixedLSC;      %fixed down-backward
% %alpha=-deg2rad(120)+fixedLSC;    %fixed up-forward
% beta=pi/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Sector 1
%This part of the code computes the vector fields in the sectors. We
%read in the numerical solution computed using the script Mapping.m to
%compute the value of rho. The vector fields are then found by computing
%the pushforward in the coordinates (u,v).

%Reading in data;
Sec1Data=dlmread('Sector1'); %radial stretching factor.
n=size(Sec1Data,1);
RHO=Sec1Data(1:n,1:n);
F=Sec1Data(1:n,n+1:2*n);

 %Record the number of points that were discretized over.

% Construction of domain in s th coordinates
s=linspace(0,M,n+1);
s=linspace(s(2),M,n);
th=linspace(a1,a2,n);


[S,TH]=meshgrid(s,th);
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

%% Normalization of vector field
veca=vecX./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);
vecb=vecY./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);
vecc=vecZ./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);

vecX=veca;
vecY=vecb;
vecZ=vecc;

ntemp=40;

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

%Writing Vector Field Data.
dlmwrite('Sec1VecField',[U,V,vecU,vecV,RHO,F,RHOS,RHOTH,FS,FTH,S,TH,vecS,vecTH,Xprime,Yprime,vecXprime,vecYprime,X,Y,Z,vecX,vecY,vecZ]);

%% Sector 2
%This part of the code computes the vector fields in the sectors. We
%read in the numerical solution computed using the script Mapping.m to
%compute the value of rho. The vector fields are then found by computing
%the pushforward in the coordinates (u,v).

%Reading in data;
Sec2Data=dlmread('Sector2'); %radial stretching factor.
n=size(Sec2Data,1);
RHO=Sec2Data(1:n,1:n);
F=Sec2Data(1:n,n+1:2*n);

 %Record the number of points that were discretized over.

% Construction of domain in s th coordinates
s=linspace(0,M,n+1);
s=linspace(s(2),M,n);
th=linspace(a2,a3,n);

[S,TH]=meshgrid(s,th);
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

%% Normalization of vector field
veca=vecX./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);
vecb=vecY./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);
vecc=vecZ./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);

vecX=veca;
vecY=vecb;
vecZ=vecc;

ntemp=40;

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

%Writing Vector Field Data.
dlmwrite('Sec2VecField',[U,V,vecU,vecV,RHO,F,RHOS,RHOTH,FS,FTH,S,TH,vecS,vecTH,Xprime,Yprime,vecXprime,vecYprime,X,Y,Z,vecX,vecY,vecZ]);


%% Sector 3
%This part of the code computes the vector fields in the sectors. We
%read in the numerical solution computed using the script Mapping.m to
%compute the value of rho. The vector fields are then found by computing
%the pushforward in the coordinates (u,v).

%Reading in data;
Sec3Data=dlmread('Sector3'); %radial stretching factor.
n=size(Sec3Data,1);
RHO=Sec3Data(1:n,1:n);
F=Sec3Data(1:n,n+1:2*n);

 %Record the number of points that were discretized over.

% Construction of domain in s th coordinates
s=linspace(0,M,n+1);
s=linspace(s(2),M,n);
th=linspace(a3,a4,n);

[S,TH]=meshgrid(s,th);
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

%% Normalization of vector field
veca=vecX./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);
vecb=vecY./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);
vecc=vecZ./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);

vecX=veca;
vecY=vecb;
vecZ=vecc;

ntemp=40;

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

%Writing Vector Field Data.
dlmwrite('Sec3VecField',[U,V,vecU,vecV,RHO,F,RHOS,RHOTH,FS,FTH,S,TH,vecS,vecTH,Xprime,Yprime,vecXprime,vecYprime,X,Y,Z,vecX,vecY,vecZ]);


%% Sector 4
%This part of the code computes the vector fields in the sectors. We
%read in the numerical solution computed using the script Mapping.m to
%compute the value of rho. The vector fields are then found by computing
%the pushforward in the coordinates (u,v).

%Reading in data;
Sec4Data=dlmread('Sector4'); %radial stretching factor.
n=size(Sec4Data,1);
RHO=Sec4Data(1:n,1:n);
F=Sec4Data(1:n,n+1:2*n);

 %Record the number of points that were discretized over.

% Construction of domain in s th coordinates
s=linspace(0,M,n+1);
s=linspace(s(2),M,n);
th=linspace(a4,a1+2*pi,n);

[S,TH]=meshgrid(s,th);
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

%% Normalization of vector field
veca=vecX./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);
vecb=vecY./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);
vecc=vecZ./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);

vecX=veca;
vecY=vecb;
vecZ=vecc;

ntemp=40;

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


%Writing Vector Field Data.
dlmwrite('Sec4VecField',[U,V,vecU,vecV,RHO,F,RHOS,RHOTH,FS,FTH,S,TH,vecS,vecTH,Xprime,Yprime,vecXprime,vecYprime,X,Y,Z,vecX,vecY,vecZ]);


