
function []=alphaCorr

%% load meridian_radius and map files, apply the alpha correction, and save as physAlphaCorr

discPoints=40;
filedir=cd;

namesMR=getFileList(filedir,'meridian_radius',0,'anywhere');
if ~isempty(namesMR)
    MR=load(num2str(namesMR{1}));
end

namesMap=getFileList(filedir,'mapAlphaCorr_',0,'anywhere');
for i=1:size(namesMap,2)   
    if strfind(namesMap{i},'.mat')
    load(num2str(namesMap{i}));
    end
end

% meridian and radius of retina
M=MR(1);
R=MR(2); 
R=R/M;

% radial lengths of cuts
m1=floor(map.D1)/M;
m2=floor(map.D2)/M;
m3=floor(map.D3)/M;
m4=floor(map.D4)/M;
nu=1/2;

%Angular coordinates of cuts
a1=map.a1;
a2=map.a2;
a3=map.a3;
a4=map.a4;

M=1;

%%Output physical parameters to be used in later scripts
%Physical parameters
phys(1)=R;
phys(2)=M;

phys(3)=m1;
phys(4)=m2;
phys(5)=m3;
phys(6)=m4;

phys(7)=nu;

phys(8)=a1;
phys(9)=a2;
phys(10)=a3;
phys(11)=a4;

dlmwrite('physAlphaCorr2',phys);

%% load sector 1-4 files
Sec1Data=dlmread(['Sector1_',num2str(discPoints)]);
Sec2Data=dlmread(['Sector2_',num2str(discPoints)]);
Sec3Data=dlmread(['Sector3_',num2str(discPoints)]);
Sec4Data=dlmread(['Sector4_',num2str(discPoints)]);

%% extract rho, f,S, and TH
rho1=Sec1Data(:,1:discPoints);
f1=Sec1Data(:,discPoints+1:2*discPoints);
S1=Sec1Data(:,2*discPoints+1:3*discPoints);
TH1=Sec1Data(:,3*discPoints+1:4*discPoints);

rho2=Sec2Data(:,1:discPoints);
f2=Sec2Data(:,discPoints+1:2*discPoints);
S2=Sec2Data(:,2*discPoints+1:3*discPoints);
TH2=Sec2Data(:,3*discPoints+1:4*discPoints);

rho3=Sec3Data(:,1:discPoints);
f3=Sec3Data(:,discPoints+1:2*discPoints);
S3=Sec3Data(:,2*discPoints+1:3*discPoints);
TH3=Sec3Data(:,3*discPoints+1:4*discPoints);

rho4=Sec4Data(:,1:discPoints);
f4=Sec4Data(:,discPoints+1:2*discPoints);
S4=Sec4Data(:,2*discPoints+1:3*discPoints);
TH4=Sec4Data(:,3*discPoints+1:4*discPoints);

%% rotate sectors 1-4 by corrFactor

fnew1=mod(f1+deg2rad(map.alphaCorr),2*pi);
fnew2=mod(f2+deg2rad(map.alphaCorr),2*pi);
fnew3=mod(f3+deg2rad(map.alphaCorr),2*pi);
fnew4=mod(f4+deg2rad(map.alphaCorr),2*pi);

THnew1=mod(TH1+deg2rad(map.alphaCorr),2*pi);
THnew2=mod(TH2+deg2rad(map.alphaCorr),2*pi);
THnew3=mod(TH3+deg2rad(map.alphaCorr),2*pi);
THnew4=mod(TH4+deg2rad(map.alphaCorr),2*pi);



%% Plot of flatened retina
n=discPoints;
set(0, 'CurrentFigure', f1) %make f1 the active  figure
surf(rho1.*cos(fnew1),rho1.*sin(fnew1),zeros(n,n));
surf(rho2.*cos(fnew2),rho2.*sin(fnew2),zeros(n,n));
surf(rho3.*cos(fnew3),rho3.*sin(fnew3),zeros(n,n));
surf(rho4.*cos(fnew4),rho4.*sin(fnew4),zeros(n,n));
view(0,90);
colormap(jet);
caxis([0,.05]);
colorbar;
shading interp;
hold on;


% U1=rho1.*cos(f1);
% V1=rho1.*sin(f1);
% unew1=U1*cos(corrFactor)-V1*sin(corrFactor);
% vnew1=U1*sin(corrFactor)+V1*cos(corrFactor);
% 
% U2=rho2.*cos(f2);
% V2=rho2.*sin(f2);
% unew2=U2*cos(corrFactor)-V2*sin(corrFactor);
% vnew2=U2*sin(corrFactor)+V2*cos(corrFactor);
% 
% U3=rho3.*cos(f3);
% V3=rho3.*sin(f3);
% unew3=U3*cos(corrFactor)-V3*sin(corrFactor);
% vnew3=U3*sin(corrFactor)+V3*cos(corrFactor);
% 
% U4=rho4.*cos(f4);
% V4=rho4.*sin(f4);
% unew4=U4*cos(corrFactor)-V4*sin(corrFactor);
% vnew4=U4*sin(corrFactor)+V4*cos(corrFactor);

    

end










