function []=ContourPlotAlphaBeta(resolution)

VizAngle=deg2rad(180);
discpoints=40;

phys=dlmread('phys');

R=phys(1);
M=phys(2);

% VizAngle should be in radians


        [FileName,PathName] = uigetfile('*.mat','Select a allPossibleFields or fitAlphaCorr file');

        if strfind(FileName,'allPossibleFields')
            load(FileName);
            
            c=10;
            cutoff=(['c',num2str(c)]);
            k=1;
            h=1;
            for Alpha=0:resolution:360
                for beta=0:resolution:180
                    mat180(k,h)=allPossibleFields.matchind180(k,h).(cutoff);
                    mat90(k,h)=allPossibleFields.matchind90(k,h).(cutoff);
                    k=k+1;
                end
                k=1;
                h=h+1;
            end
            
            data=mat180;
           
            
        elseif strfind(FileName,'fitAlphaCorr')
            load(FileName);
            
            data=CfitRs180;
            %data=mat180ASC;
            %data=mat180PSC;
            %data=mat180LSC;
            
            
            
            
            
        end
    

%% Standard Contour
[a,b]=size(data);

Alpha1=0.6257+pi;
beta1=(pi/2-1.2624)+pi/2;
Alpha2=0.8848;
beta2=3.0135;
Alpha3=-1.0012+pi;
beta3=(pi/2-1.14502)+pi/2;

Alpha=linspace(0,2*pi,b+1);
Alpha=Alpha(2:b+1);

beta=linspace(pi-VizAngle/2,pi,a);


[A,B]=meshgrid(Alpha,beta);

% figure;
% hold on;
% contourf(A,B,data,'EdgeColor','none','LineStyle','none');
% plot(Alpha1,beta1,'ok',Alpha2,beta2,'ok',Alpha3,beta3,'ok');
% hold off;
% shading interp;

%% Contour on Sphere

X=R*cos(A).*sin(B);
Y=R*sin(A).*sin(B);
Z=-R*cos(B);

% shift array by 30 degrees (2 positions)
% NEED TO COMMENT
% X = circshift(X,2,2);
% Y = circshift(Y,2,2);
% Z = circshift(Z,2,2);


Xshift=X;
Yshift=Y;
Zshift=Z;
datashift=data;

Xshift(:,b+1)=X(:,1);
Yshift(:,b+1)=Y(:,1);
Zshift(:,b+1)=Z(:,1);
datashift(:,b+1)=data(:,1);



figure('Color',[1 1 1]);
hold on;
surf(Xshift,Yshift,Zshift,datashift,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
shading interp;
axis equal

% Plot mesh outside of sphere
% thetafine=linspace(0,2*pi,40);
% sfine=linspace(pi/2,pi,10);
% 
% [THfine,Sfine]=meshgrid(thetafine,sfine);
% 
% Xfine=(R+0.01)*cos(THfine).*sin(Sfine);
% Yfine=(R+0.01)*sin(THfine).*sin(Sfine);
% Zfine=-(R+0.01)*cos(Sfine);
% 
% surf(Xfine,Yfine,Zfine,'FaceColor','none');


% Plot mesh outside of sphere
thetafine=linspace(0,2*pi,40);
sfine=linspace(pi/2,pi,10);

[THfine,Sfine]=meshgrid(thetafine,sfine);

Xfine=(R-0.01)*cos(THfine).*sin(Sfine);
Yfine=(R-0.01)*sin(THfine).*sin(Sfine);
Zfine=-(R-0.01)*cos(Sfine);

surf(Xfine,Yfine,Zfine,'FaceColor','none','EdgeAlpha',0.1);
%%
%%%%% plot symbols for canal locations

normalASC=[0.592,0.764,0.247];     %ASC
%normalASC=[0.724,0.615,0.300];      %ASC Prime
[alphaASC,betaASC]=AB(normalASC);
normalPSC=[0.470,-0.764,0.438];       %PSC
%normalPSC=[0.678,-0.638,0.359];       %PSC Prime
[alphaPSC,betaPSC]=AB(normalPSC);
normalLSC=[0.421,-0.056,-0.901];       %LSC
%normalLSC=[0.533,-0.102,-0.829];       %LSC Prime
[alphaLSC,betaLSC]=AB(normalLSC);


% calculate the reciprocal singularities
if alphaASC>180
    alpha2ASC=alphaASC-pi;
else
    alpha2ASC=alphaASC+pi;
end
beta2ASC=pi/2-(betaASC-pi/2);
    
if alphaPSC>180
    alpha2PSC=alphaPSC-pi;
else
    alpha2PSC=alphaPSC+pi;
end
beta2PSC=pi/2-(betaPSC-pi/2);

if alphaLSC>180
    alpha2LSC=alphaLSC-pi;
else
    alpha2LSC=alphaLSC+pi;
end
beta2LSC=pi/2-(betaLSC-pi/2);


% inside sphere ASC
xASCi=(R-0.01)*cos(alphaASC)*sin(betaASC);
yASCi=(R-0.01)*sin(alphaASC)*sin(betaASC);
zASCi=-(R-0.01)*cos(betaASC);
hASCi=plot3(xASCi,yASCi,zASCi,'ok','MarkerSize',10, 'LineWidth',1','MarkerFaceColor','b');

% inside sphere ASC secondary singularity
xASC2i=(R-0.01)*cos(alpha2ASC)*sin(beta2ASC);
yASC2i=(R-0.01)*sin(alpha2ASC)*sin(beta2ASC);
zASC2i=-(R-0.01)*cos(beta2ASC);
hASC2i=plot3(xASC2i,yASC2i,zASC2i,'+','MarkerSize',5, 'LineWidth',5','MarkerEdgeColor','b');

% inside sphere PSC
xPSCi=(R-0.01)*cos(alphaPSC)*sin(betaPSC);
yPSCi=(R-0.01)*sin(alphaPSC)*sin(betaPSC);
zPSCi=-(R-0.01)*cos(betaPSC);
hPSCi=plot3(xPSCi,yPSCi,zPSCi,'ok','MarkerSize',10, 'LineWidth',1','MarkerFaceColor','m');

% inside sphere PSC secondary singularity
xPSC2i=(R-0.01)*cos(alpha2PSC)*sin(beta2PSC);
yPSC2i=(R-0.01)*sin(alpha2PSC)*sin(beta2PSC);
zPSC2i=-(R-0.01)*cos(beta2PSC);
hPSC2i=plot3(xPSC2i,yPSC2i,zPSC2i,'+','MarkerSize',5, 'LineWidth',5,'MarkerEdgeColor','m');

% inside sphere LSC
xLSCi=(R-0.01)*cos(alphaLSC)*sin(betaLSC);
yLSCi=(R-0.01)*sin(alphaLSC)*sin(betaLSC);
zLSCi=-(R-0.01)*cos(betaLSC);
hLSCi=plot3(xLSCi,yLSCi,zLSCi,'ok','MarkerSize',10, 'LineWidth',1','MarkerFaceColor','r');

% inside sphere LSC secondary singularity
xLSC2i=(R-0.01)*cos(alpha2LSC)*sin(beta2LSC);
yLSC2i=(R-0.01)*sin(alpha2LSC)*sin(beta2LSC);
zLSC2i=-(R-0.01)*cos(beta2LSC);
hLSC2i=plot3(xLSC2i,yLSC2i,zLSC2i,'+','MarkerSize',5, 'LineWidth',5','MarkerEdgeColor','r');

% plot a marker for the optic disc
hOpt=plot3(0,0,0.5,'ok','MarkerSize',10, 'LineWidth',1','MarkerFaceColor','k');

%%

az = 0;
el = -90;
view(az, el);
set(gca, 'visible', 'off');
set(gca,'color',[1 1 1]);


%% draw the normal of canals as lines
lb=0.52;
ub=1;


lb=-1;
ub=1;

lambda=linspace(lb,ub,100);

%plot axis of ASC
xl=lambda*cos(alphaASC)*sin(betaASC);
yl=lambda*sin(alphaASC)*sin(betaASC);
zl=-lambda*cos(betaASC);
h=plot3(xl,yl,zl,'b');
set(h,'LineWidth',5);

xl=lambda*cos(alphaPSC)*sin(betaPSC);
yl=lambda*sin(alphaPSC)*sin(betaPSC);
zl=-lambda*cos(betaPSC);
h=plot3(xl,yl,zl,'m');
set(h,'LineWidth',5);

xl=lambda*cos(alphaLSC)*sin(betaLSC);
yl=lambda*sin(alphaLSC)*sin(betaLSC);
zl=-lambda*cos(betaLSC);
h=plot3(xl,yl,zl,'r');
set(h,'LineWidth',5);

% outside sphere ASC
xASCo=(R+0.005)*cos(alphaASC)*sin(betaASC);
yASCo=(R+0.005)*sin(alphaASC)*sin(betaASC);
zASCo=-(R+0.005)*cos(betaASC);
hASCo=plot3(xASCo,yASCo,zASCo,'ok','MarkerSize',7, 'LineWidth',1, 'MarkerFaceColor','b');

% outside sphere PSC
xPSCo=(R+0.005)*cos(alphaPSC)*sin(betaPSC);
yPSCo=(R+0.005)*sin(alphaPSC)*sin(betaPSC);
zPSCo=-(R+0.005)*cos(betaPSC);
hPSCo=plot3(xPSCo,yPSCo,zPSCo,'ok','MarkerSize',7, 'LineWidth',1, 'MarkerFaceColor','m');

% outside sphere LSC
xLSCo=(R+0.005)*cos(alphaLSC)*sin(betaLSC);
yLSCo=(R+0.005)*sin(alphaLSC)*sin(betaLSC);
zLSCo=-(R+0.005)*cos(betaLSC);
hLSCo=plot3(xLSCo,yLSCo,zLSCo,'ok','MarkerSize',7, 'LineWidth',1, 'MarkerFaceColor','r');

hold off;


%% Contour on Standard retina

%Mapping to S-TH coordinates
Sdata=-2*M/VizAngle*acos(-Z/R)+2*M*pi/VizAngle;
THdata=atan2(Y,X)+pi;


StSec1Data=dlmread(['StandardSector1_',num2str(discpoints)]);
StSec2Data=dlmread(['StandardSector2_',num2str(discpoints)]);
StSec3Data=dlmread(['StandardSector3_',num2str(discpoints)]);
StSec4Data=dlmread(['StandardSector4_',num2str(discpoints)]);

n=size(StSec1Data,1);


StS1=StSec1Data(1:n,2*n+1:3*n);
StTH1=StSec1Data(1:n,3*n+1:4*n);
StRHO1=StSec1Data(1:n,1:n);
StF1=StSec1Data(1:n,n+1:2*n);

StU1=StRHO1.*cos(StF1);
StV1=StRHO1.*sin(StF1);


StS2=StSec2Data(1:n,2*n+1:3*n);
StTH2=StSec2Data(1:n,3*n+1:4*n);
StRHO2=StSec2Data(1:n,1:n);
StF2=StSec2Data(1:n,n+1:2*n);

StU2=StRHO2.*cos(StF2);
StV2=StRHO2.*sin(StF2);


StS3=StSec3Data(1:n,2*n+1:3*n);
StTH3=StSec3Data(1:n,3*n+1:4*n);
StRHO3=StSec3Data(1:n,1:n);
StF3=StSec3Data(1:n,n+1:2*n);

StU3=StRHO3.*cos(StF3);
StV3=StRHO3.*sin(StF3);


StS4=StSec4Data(1:n,2*n+1:3*n);
StTH4=StSec4Data(1:n,3*n+1:4*n);
StRHO4=StSec4Data(1:n,1:n);
StF4=StSec4Data(1:n,n+1:2*n);

StU4=StRHO4.*cos(StF4);
StV4=StRHO4.*sin(StF4);

for i=1:a,
    for j=1:b,
        if THdata(i,j)>=0 &&  THdata(i,j)<pi/2,
            StUdata(i,j)=griddata(StS1,StTH1,StU1,Sdata(i,j),THdata(i,j));
            StVdata(i,j)=griddata(StS1,StTH1,StV1,Sdata(i,j),THdata(i,j));
        elseif THdata(i,j)>=pi/2 &&  THdata(i,j)<pi,
          StUdata(i,j)=griddata(StS2,StTH2,StU2,Sdata(i,j),THdata(i,j));
          StVdata(i,j)=griddata(StS2,StTH2,StV2,Sdata(i,j),THdata(i,j));
        elseif THdata(i,j)>=pi &&  THdata(i,j)<3*pi/2,
          StUdata(i,j)=griddata(StS3,StTH3,StU3,Sdata(i,j),THdata(i,j));
          StVdata(i,j)=griddata(StS3,StTH3,StV3,Sdata(i,j),THdata(i,j)); 
        else
            StUdata(i,j)=griddata(StS4,StTH4,StU4,Sdata(i,j),THdata(i,j));
            StVdata(i,j)=griddata(StS4,StTH4,StV4,Sdata(i,j),THdata(i,j));
        end
    end
end

StUdata(isnan(StUdata))=0;
StVdata(isnan(StVdata))=0;



StUdata1=StUdata;
StVdata1=StVdata;
data1=data;

StUdata1(~(THdata>=0 &  THdata<pi/2))=NaN;
StVdata1(~(THdata>=0 &  THdata<pi/2))=NaN;
data1(~(THdata>=0 &  THdata<pi/2))=NaN;

StUdata2=StUdata;
StVdata2=StVdata;
data2=data;

StUdata2(~(THdata>=pi/2 &  THdata<pi))=NaN;
StVdata2(~(THdata>=pi/2 &  THdata<pi))=NaN;
data2(~(THdata>=pi/2 &  THdata<pi))=NaN;

StUdata3=StUdata;
StVdata3=StVdata;
data3=data;

StUdata3(~(THdata>=pi &  THdata<3*pi/2))=NaN;
StVdata3(~(THdata>=pi &  THdata<3*pi/2))=NaN;
data3(~(THdata>=pi &  THdata<3*pi/2))=NaN;

StUdata4=StUdata;
StVdata4=StVdata;
data4=data;

StUdata4(~(THdata>=3*pi/2 &  THdata<2*pi))=NaN;
StVdata4(~(THdata>=3*pi/2 &  THdata<2*pi))=NaN;
data4(~(THdata>=3*pi/2 &  THdata<2*pi))=NaN;

figure;
hold on;
contourf(StUdata1,StVdata1,data1,'EdgeColor','none','LineStyle','none');
contourf(StUdata2,StVdata2,data2,'EdgeColor','none','LineStyle','none');
contourf(StUdata3,StVdata3,data3,'EdgeColor','none','LineStyle','none');
contourf(StUdata4,StVdata4,data4,'EdgeColor','none','LineStyle','none');
%plot(Alpha1,beta1,'ok',Alpha2,beta2,'ok',Alpha3,beta3,'ok');
hold off;

shading interp;
