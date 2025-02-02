function [allPossibleFields]= VecComparisonPooled4TransBoot2(data,a,maxVec,alpha,beta,polarity,fieldID,VizAngle,allPossibleFields,k,h,interpType,histType)

% histType can be either 1 (predicted vector for ASC), 2 (predicted vector
% for PSC), 3 (predicted vector for LSC), 4 2 (predicted vector for all
% three canals), or 5 (real preferred direction of cells).

% cellType can be either 'ONDS', 'ONOFFDS', 'OFFDS', or 'retroONDS'.
% fieldID can be either 'ASC', 'PSC', 'LSC', 'cardinalUp', etc.
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

VizAngle=deg2rad(VizAngle);

%% Input measured data mapped to XYZ coordinates
% import XYZ data.

Xt=data.XYZvecComp.X';
Yt=data.XYZvecComp.Y';
Zt=data.XYZvecComp.Z';
VecXt=data.XYZvecComp.VecX';
VecYt=data.XYZvecComp.VecY';
VecZt=data.XYZvecComp.VecZ';

nsam=size(Xt,2);


% resample with replacement
[X,idx]=datasample(Xt,nsam);
Y=Yt(idx);
Z=Zt(idx);
VecX=VecXt(idx)./(VecXt(idx).^2+VecYt(idx).^2+VecZt(idx).^2).^(1/2);
VecY=VecYt(idx)./(VecXt(idx).^2+VecYt(idx).^2+VecZt(idx).^2).^(1/2);
VecZ=VecZt(idx)./(VecXt(idx).^2+VecYt(idx).^2+VecZt(idx).^2).^(1/2);

% X=Xt;
% Y=Yt;
% Z=Zt;
% VecX=VecXt./(VecXt.^2+VecYt.^2+VecZt.^2).^(1/2);
% VecY=VecYt./(VecXt.^2+VecYt.^2+VecZt.^2).^(1/2);
% VecZ=VecZt./(VecXt.^2+VecYt.^2+VecZt.^2).^(1/2);

%% Creating Vector at data points
S=R*asin((X.^2+Y.^2).^(1/2)/R);
S=real(S);
TH=atan2(Y,X);

SV=-VizAngle*R/(2*M)*S+pi*R;
THV=TH+pi;


Npooled=size(X,2);

t=zeros(1,Npooled);
phi=zeros(1,Npooled);
vecSV=zeros(1,Npooled);
vecTHV=zeros(1,Npooled);

for i=1:Npooled
    [t(i),phi(i)]=InverseSTH(alpha,beta,R,SV(i),THV(i));
    
    vecSV(i)=sin(SV(i)/R)*(cos(beta)*sin(t(i)/R)+...
        cos(t(i)/R)*cos(phi(i))*sin(beta));
    vecSV(i)=vecSV(i)+cos(SV(i)/R)*(-cos(alpha-THV(i))*sin(t(i)/R)*sin(beta)+...
        cos(t(i)/R)*(cos(beta)*cos(alpha-THV(i))*cos(phi(i))-...
        sin(alpha-THV(i))*sin(phi(i))));
    
    vecTHV(i)=-sin(t(i)/R)*sin(beta)*sin(alpha-THV(i))+cos(t(i)/R)*...
        (cos(beta)*cos(phi(i))*sin(alpha-THV(i))+cos(alpha-THV(i))*sin(phi(i)));
    vecTHV(i)=1/R*csc(SV(i)/R)*vecTHV(i);
    
end

temp=(vecSV.^2+vecTHV.^2).^(1/2);

vecSV=vecSV./temp;
vecTHV=vecTHV./temp;

%quiver3(XV,YV,ZV,vecXV,vecYV,vecZV);

temp=(vecSV.^2+vecTHV.^2).^(1/2);
vecSV=vecSV./temp;
vecTHV=vecTHV./temp;

vecS=-2*M/(VizAngle*R)*vecSV;

vecTH=vecTHV.*polarity;

vecX=vecS.*cos(S/R).*cos(TH)-R*vecTH.*sin(S/R).*sin(TH);
vecY=vecS.*cos(S/R).*sin(TH)+R*vecTH.*sin(S/R).*cos(TH);
vecZ=vecS.*sin(S/R);

% Normalization of vector field
veca=vecX./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);
vecb=vecY./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);
vecc=vecZ./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);

VecXcomp=veca;
VecYcomp=vecb;
VecZcomp=vecc;


%% calculate matchind matrices
if ~isempty(VecX)
    matchind180=[];
    matchindMix=[];
    angleXYZ=rad2deg(real(acos(VecX.*VecXcomp+VecY.*VecYcomp+VecZ.*VecZcomp)));
    
    for c=10:10:10
        cutoff=(['c',num2str(c)]);
        matchind180.(cutoff)=100*(numel(find(angleXYZ<c))/numel(angleXYZ));
        matchindMix.(cutoff)=100*(((numel(intersect(find(angleXYZ>=90-c),find(angleXYZ<90+c)))/2)+numel(find(angleXYZ<c)))/numel(angleXYZ));
    end
end

%%  save
allPossibleFields.field{k,h}=fieldID;
allPossibleFields.anglediff180{k,h}=angleXYZ';
allPossibleFields.VecXcomp{k,h}=VecXcomp';
allPossibleFields.VecYcomp{k,h}=VecYcomp';
allPossibleFields.VecZcomp{k,h}=VecZcomp';
allPossibleFields.matchind180(k,h)=matchind180;
allPossibleFields.matchindMix(k,h)=matchindMix;
if k==1 && h==1
    allPossibleFields.X=X';
    allPossibleFields.Y=Y';
    allPossibleFields.Z=Z';
    allPossibleFields.VecX=VecX';
    allPossibleFields.VecY=VecY';
    allPossibleFields.VecZ=VecZ';
end

end


