function [allPossibleFields]= VecComparisonPooled2(pooledmap,alpha,beta,polarity,fieldType,VizAngle,allPossibleFields,k,h)

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

n=100; %Number of points in discretization of vector field on the sphere.

%% Creating Vector field in XYZ coordinates.
% Construction of domain in s th coordinates
s=linspace(0,M,n+1);
s=linspace(s(2),M,n);
th=linspace(0,2*pi,n);


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

vecS=-2*M/(VizAngle*R)*vecS;

vecTH=vecTH.*polarity;

vecX=vecS.*cos(S/R).*cos(TH)-R*vecTH.*sin(S/R).*sin(TH);
vecY=vecS.*cos(S/R).*sin(TH)+R*vecTH.*sin(S/R).*cos(TH);
vecZ=vecS.*sin(S/R);

% Normalization of vector field
veca=vecX./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);
vecb=vecY./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);
vecc=vecZ./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);

XM=R*cos(TH)*sin(S/R);
YM=R*sin(TH)*sin(S/R);
ZM=R-R*cos(S/R);

vecXM=veca;
vecYM=vecb;
vecZM=vecc;

%% Input measured data mapped to XYZ coordinates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   We now input measured Data from the retina in the flattened
%   coordinates.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% import XYZ data.
X=pooledmap.XYZvecComp.X';
Y=pooledmap.XYZvecComp.Y'; 
Z=pooledmap.XYZvecComp.Z'; 
VecX=pooledmap.XYZvecComp.VecX'; 
VecY=pooledmap.XYZvecComp.VecY'; 
VecZ=pooledmap.XYZvecComp.VecZ'; 


n=100; %Number of points used in discretization of vector field.

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

Npooled=size(X,2);
for i=1:Npooled
            
            %XYZ data
            VecXcomp(i)=VecXMInterp(X(i),Y(i),Z(i));
            VecYcomp(i)=VecYMInterp(X(i),Y(i),Z(i));
            VecZcomp(i)=VecZMInterp(X(i),Y(i),Z(i));
                      
end

%% Normalizing Vector Fields

temp=(VecXcomp.^2+VecYcomp.^2+VecZcomp.^2).^(1/2);
VecXcomp=VecXcomp./temp;
VecYcomp=VecYcomp./temp;
VecZcomp=VecZcomp./temp;


%% Computing Angles Between Vectors

matchind180=[];
angleXYZ=rad2deg(acos(VecX.*VecXcomp+VecY.*VecYcomp+VecZ.*VecZcomp));
for c=5:5:45
    cutoff=(['c',num2str(c)]);
    matchind180.(cutoff)=100*(numel(find(angleXYZ<c))/numel(angleXYZ));
end
matchind90=[];
angleXYZ90=[angleXYZ(find(angleXYZ<=90)),180-angleXYZ(find(angleXYZ>90))];
for c=5:5:45
    cutoff=(['c',num2str(c)]);
matchind90.(cutoff)=100*(numel(find(angleXYZ90<c))/numel(angleXYZ));
end

%figure;
subplot(1,2,1); histogram(angleXYZ,15); 
ax=gca; xlim([0 180]); ax.XTick=[0:45:180];
subplot(1,2,2); histogram(angleXYZ90,15);
ax=gca; xlim([0 90]); ax.XTick=[0:45:90];
drawnow expose

allPossibleFields.field{k,h}=fieldType;
allPossibleFields.anglediff180{k,h}=angleXYZ';
allPossibleFields.anglediff90{k,h}=angleXYZ90';
allPossibleFields.VecXcomp{k,h}=VecXcomp';
allPossibleFields.VecYcomp{k,h}=VecYcomp';
allPossibleFields.VecZcomp{k,h}=VecZcomp';
allPossibleFields.matchind180(k,h)=matchind180;
allPossibleFields.matchind90(k,h)=matchind90;
if k==1 && h==1
allPossibleFields.VecX=VecX';
allPossibleFields.VecY=VecY';
allPossibleFields.VecZ=VecZ';
end




end

