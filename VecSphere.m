%VecSphere

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This script maps Shai's data back to a sphere and plots the resulting
%   vector field on the sphere.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Physical Parameters.
% This data is imported from data exported from the mapping script. These
% values correspond to measured values on the flattened retina.

clear all;
close all;

figure;
hold on;
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

%% Collect data in the (S,TH) coordiante system.

[filedir,go]=bbdirselector('select data folder',cd);
if ~go disp('no folder selected'); return; end

names=getFileList(filedir,'STHData',0,'anywhere');
        if ~isempty(names)
            STHdata=dlmread(cell2mat(names));
        end

%STHdata=dlmread('STHData');

Sn=STHdata(:,1);
THn=STHdata(:,2);
vecSn=STHdata(:,3);
vecTHn=STHdata(:,4);

Sn=Sn/(2*pi)*M; %This converts the angular representation of S back into a length 

%% Cooredinates in cartesian space cooresponding to s and th
x=R*sin(Sn/R).*cos(THn);
y=R*sin(Sn/R).*sin(THn);
z=-R*cos(Sn/R);

%% Vector components in cartesian space
vecX=vecSn.*cos(Sn/R).*cos(THn)-R*vecTHn.*sin(Sn/R).*sin(THn);
vecY=vecSn.*cos(Sn/R).*sin(THn)+R*vecTHn.*sin(Sn/R).*cos(THn);
vecZ=vecSn.*sin(Sn/R);

%% Normalization of vector field
veca=vecX./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);
vecb=vecY./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);
vecc=vecZ./(vecX.^2+vecY.^2+vecZ.^2).^(1/2);

vecX=veca;
vecY=vecb;
vecZ=vecc;

%% Plotting vector
figure;
hold on;
quiver3(x,y,z,vecX,vecY,vecZ,.3,'linewidth',2);


%% Plotting retina
n=40;

s=linspace(0,M,n);
th=linspace(0,2*pi,n);

[S,TH]=meshgrid(s,th);

surf(R*sin(S/R).*cos(TH),R*sin(S/R).*sin(TH),-R*cos(S/R),S,'EdgeColor','none','FaceAlpha','interp','AlphaDataMapping','scaled',...
    'AlphaData',S);
alphamap('default')

hold off;




