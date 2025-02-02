% VecComparisonTest
%function [allPossibleFields]= VecComparisonPooledTest(pooledmap,alpha,beta,polarity,fieldType,VizAngle,allPossibleFields,k,h)

% cellType can be either 'ONDS', 'ONOFFDS', 'OFFDS', or 'retroONDS'.
% fieldType can be either 'ASC', 'PSC', 'LSC', 'cardinalUp', etc.
%VecComparison

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This script compares the vector fields obtained from the elastic model
%   with that from Shai's data. We compare the data in the four releveant
%   coordinate systems for the problem. We write the data to four seperate
%   files called
%`
%   1. angleDifUV
%   2. angleDifSTH
%   3. angleDifGSTH
%   4. angleDifXYZ
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

%% Input Parameters
VizAngle=180;
normalASC=[0.592,0.764,0.247];     %ASC
[alpha,beta]=AB(normalASC);
polarity=1;
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

VizAngle=deg2rad(VizAngle);

n=40; %Number of points in discretization of vector field on the sphere.


%% Creating Vector field in XYZ coordinates.
% Construction of domain in s th coordinates
s=linspace(0,M,n+1);
s=linspace(s(2),M,n);
th=linspace(0,2*pi,n);


[S,TH]=meshgrid(s,th);
%Calculation of the vector field of rotations in the coordinates s and th.
%The variables (t,phi) denote geodesic polar coordinate around the axis of
%rotation. 

%% Creating Vector field in Visual Space coordinates.
SV=-VizAngle*R/(2*M)*S+pi*R;
THV=TH+pi;



t=zeros(n,n);
phi=zeros(n,n);

%This loop constructs the components of the vector field of rotations
%in the the s, th coordinates. The function InverseSTH determines the
%coordinates of (t,phi) given coordiantes s and th and is necessary to
%calculate the components of the vector field. 


%This loop constructs the components of the vector field of rotations
%in the the s, th coordinates. The function InverseSTH determines the
%coordinates of (t,phi) given coordinates s and th and is necessary to
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

temp=(vecSV.^2+vecTHV.^2).^(1/2);

vecSV=vecSV./temp;
vecTHV=vecTHV./temp;

figure;
hold on;
quiver(SV,THV,vecSV,vecTHV,.2);
plot(pi*R-R*beta,alpha+pi,'o');
hold off;

figure;
XV=R*sin(SV/R).*cos(THV);
YV=R*sin(SV/R).*sin(THV);
ZV=-R*cos(SV/R);

vecXV=vecSV.*cos(SV/R).*cos(THV)-R*vecTHV.*sin(SV/R).*sin(THV);
vecYV=vecSV.*cos(SV/R).*sin(THV)+R*vecTHV.*sin(SV/R).*cos(THV);
vecZV=vecSV.*sin(SV/R);

quiver3(XV,YV,ZV,vecXV,vecYV,vecZV);

temp=(vecSV.^2+vecTHV.^2).^(1/2);
vecSV=vecSV./temp;
vecTHV=vecTHV./temp;

vecS=-2*M/(VizAngle*R)*vecSV;

vecTH=vecTHV.*polarity;

figure;
quiver(S,TH,vecS,vecTH,.2);

vecX=vecS.*cos(S/R).*cos(TH)-R*vecTH.*sin(S/R).*sin(TH);
vecY=vecS.*cos(S/R).*sin(TH)+R*vecTH.*sin(S/R).*cos(TH);
vecZ=vecS.*sin(S/R);

% Normalization of vector field
veca=vecX./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);
vecb=vecY./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);
vecc=vecZ./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);

XM=R.*cos(TH).*sin(S./R);
YM=R.*sin(TH).*sin(S./R);
ZM=R-R.*cos(S./R);

vecXM=veca;
vecYM=vecb;
vecZM=vecc;

figure;
quiver3(XM,YM,ZM,vecXM,vecYM,vecZM);

load pooledmap.mat
% model ASC canal
X=pooledmap.XYZvecComp.X';
Y=pooledmap.XYZvecComp.Y'; 
Z=pooledmap.XYZvecComp.Z'; 
VecX=pooledmap.XYZvecComp.VecXcompASC'; 
VecY=pooledmap.XYZvecComp.VecYcompASC'; 
VecZ=pooledmap.XYZvecComp.VecZcompASC'; 

% % model ASC and LSC canal
% X=[pooledmap.XYZvecComp.X';pooledmap.XYZvecComp.X'];
% Y=[pooledmap.XYZvecComp.Y';pooledmap.XYZvecComp.Y']; 
% Z=[pooledmap.XYZvecComp.Z';pooledmap.XYZvecComp.Z']; 
% VecX=[pooledmap.XYZvecComp.VecXcompASC';pooledmap.XYZvecComp.VecXcompLSC']; 
% VecY=[pooledmap.XYZvecComp.VecYcompASC';pooledmap.XYZvecComp.VecYcompLSC']; 
% VecZ=[pooledmap.XYZvecComp.VecZcompASC';pooledmap.XYZvecComp.VecZcompLSC']; 



%% Interpolating vector field at data points

XMshaped=reshape(XM,n^2,1);
YMshaped=reshape(YM,n^2,1);
ZMshaped=reshape(ZM,n^2,1);
VecXMshaped=reshape(vecXM,n^2,1);
VecYMshaped=reshape(vecYM,n^2,1);
VecZMshaped=reshape(vecZM,n^2,1);

VecXMInterp=scatteredInterpolant(XMshaped,YMshaped,ZMshaped,VecXMshaped);
VecYMInterp=scatteredInterpolant(XMshaped,YMshaped,ZMshaped,VecYMshaped);
VecZMInterp=scatteredInterpolant(XMshaped,YMshaped,ZMshaped,VecZMshaped);

figure;
quiver3(XM,YM,ZM,vecXM,vecYM,vecZM);

Npooled=size(X,2);
for i=1:Npooled
            
            %XYZ data
            VecXcomp(i)=VecXMInterp(X(i),Y(i),Z(i));
            VecYcomp(i)=VecYMInterp(X(i),Y(i),Z(i));
            VecZcomp(i)=VecZMInterp(X(i),Y(i),Z(i));
                      
end



figure;
quiver3(X,Y,Z,VecX,VecY,VecZ);

%% Normalizing Vector Fields

temp=(VecXcomp.^2+VecYcomp.^2+VecZcomp.^2).^(1/2);
VecXcomp=VecXcomp./temp;
VecYcomp=VecYcomp./temp;
VecZcomp=VecZcomp./temp;

figure;
quiver3(X,Y,Z,VecX,VecY,VecZ,'r');
hold on
quiver3(X,Y,Z,VecXcomp,VecYcomp,VecZcomp,'b');
