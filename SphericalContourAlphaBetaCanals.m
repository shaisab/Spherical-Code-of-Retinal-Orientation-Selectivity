function []=SphericalContourAlphaBetaCanals(resolution)

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
    
%% interpolate data
%resolution=5;
alphaold=0:resolution:360;
betaold=0:resolution:180;
[Aold,Bold]=meshgrid(alphaold,betaold);

newRes=1;
alphanew=0:newRes:360;
betanew=0:newRes:180;
[Anew,Bnew]=meshgrid(alphanew,betanew);
data=griddata(Aold,Bold,data,Anew,Bnew);

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

%beta=linspace(pi-VizAngle/2,2*pi,a);
beta=linspace(0,pi,a);


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
colormap jet
xlabel('X');ylabel('Y');zlabel('Z');

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
sfine=linspace(0,pi,20);

[THfine,Sfine]=meshgrid(thetafine,sfine);

Xfine=(R+0.005)*cos(THfine).*sin(Sfine);
Yfine=(R+0.005)*sin(THfine).*sin(Sfine);
Zfine=-(R+0.005)*cos(Sfine);

surf(Xfine,Yfine,Zfine,'FaceColor','none','EdgeAlpha',0.1);

%% draw shaded hemisphere
c=ceil(size(beta,2)/2)-1;
d=size(beta,2);

betaS=beta(c:d);
[A,B]=meshgrid(Alpha,betaS);

X=R*cos(A).*sin(B);
Y=R*sin(A).*sin(B);
Z=-R*cos(B);

Xshift=X;
Yshift=Y;
Zshift=Z;
datashift=data(c:d,:);

Xshift(:,b+1)=X(:,1);
Yshift(:,b+1)=Y(:,1);
Zshift(:,b+1)=Z(:,1);
datashift(:,b+1)=data(c:d,1);

datashift=ones(size(datashift,1),size(datashift,2));

hold on;
surf(Xshift,Yshift,Zshift,datashift,'EdgeColor','none','LineStyle','none','FaceColor',[1,1,1],'FaceAlpha',0.2);


%% plot symbols for canal locations
%normalASC=[0.596,0.766,0.241];     %ASC
normalASC=[0.604,0.758,0.243];     %ASC Symmetric
[alphaASC,betaASC]=AB(normalASC);
%normalPSC=[0.471,-0.766,0.438];       %PSC
normalPSC=[0.447,-0.76,0.469];       %PSC Symmetric
[alphaPSC,betaPSC]=AB(normalPSC);
%normalLSC=[0.421,-0.060,-0.905];       %LSC
normalLSC=[0.449,-0.085,-0.886];       %LSC Symmetric
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
 
 
% % inside sphere ASC
% xASCi=(R-0.01)*cos(alphaASC)*sin(betaASC);
% yASCi=(R-0.01)*sin(alphaASC)*sin(betaASC);
% zASCi=-(R-0.01)*cos(betaASC);
% hASCi=plot3(xASCi,yASCi,zASCi,'ok','MarkerSize',10, 'LineWidth',1','MarkerFaceColor','b');
% 
% % inside sphere ASC secondary singularity
% xASC2i=(R-0.01)*cos(alpha2ASC)*sin(beta2ASC);
% yASC2i=(R-0.01)*sin(alpha2ASC)*sin(beta2ASC);
% zASC2i=-(R-0.01)*cos(beta2ASC);
% hASC2i=plot3(xASC2i,yASC2i,zASC2i,'+','MarkerSize',5, 'LineWidth',5','MarkerEdgeColor','b');
% 
% % inside sphere PSC
% xPSCi=(R-0.01)*cos(alphaPSC)*sin(betaPSC);
% yPSCi=(R-0.01)*sin(alphaPSC)*sin(betaPSC);
% zPSCi=-(R-0.01)*cos(betaPSC);
% hPSCi=plot3(xPSCi,yPSCi,zPSCi,'ok','MarkerSize',10, 'LineWidth',1','MarkerFaceColor','m');
% 
% % inside sphere PSC secondary singularity
% xPSC2i=(R-0.01)*cos(alpha2PSC)*sin(beta2PSC);
% yPSC2i=(R-0.01)*sin(alpha2PSC)*sin(beta2PSC);
% zPSC2i=-(R-0.01)*cos(beta2PSC);
% hPSC2i=plot3(xPSC2i,yPSC2i,zPSC2i,'+','MarkerSize',5, 'LineWidth',5,'MarkerEdgeColor','m');
% 
% % inside sphere LSC
% xLSCi=(R-0.01)*cos(alphaLSC)*sin(betaLSC);
% yLSCi=(R-0.01)*sin(alphaLSC)*sin(betaLSC);
% zLSCi=-(R-0.01)*cos(betaLSC);
% hLSCi=plot3(xLSCi,yLSCi,zLSCi,'ok','MarkerSize',10, 'LineWidth',1','MarkerFaceColor','r');
% 
% % inside sphere LSC secondary singularity
% xLSC2i=(R-0.01)*cos(alpha2LSC)*sin(beta2LSC);
% yLSC2i=(R-0.01)*sin(alpha2LSC)*sin(beta2LSC);
% zLSC2i=-(R-0.01)*cos(beta2LSC);
% hLSC2i=plot3(xLSC2i,yLSC2i,zLSC2i,'+','MarkerSize',5, 'LineWidth',5','MarkerEdgeColor','r');

% plot a marker for the optic disc
hOpt=plot3(0,0,0.5,'ok','MarkerSize',10, 'LineWidth',1','MarkerFaceColor','k');

%%

az = -45;
el = 225;
view(az, el);
set(gca, 'visible', 'off');
set(gca,'color',[1 1 1]);


%% draw the normal of canals as lines
lb=-0.75;
ub=0.75;


% lb=-1;
% ub=1;

lambda=linspace(lb,ub,100);

%plot axis of ASC
xl=lambda*cos(alphaASC)*sin(betaASC);
yl=lambda*sin(alphaASC)*sin(betaASC);
zl=-lambda*cos(betaASC);
hASC=plot3(xl,yl,zl,'b');
set(hASC,'LineWidth',5);
Cylinder([xl(1),yl(1),zl(1)],[xl(end),yl(end),zl(end)],0.01,20,'b',1,0)

%plot axis of PSC
xl=lambda*cos(alphaPSC)*sin(betaPSC);
yl=lambda*sin(alphaPSC)*sin(betaPSC);
zl=-lambda*cos(betaPSC);
hPSC=plot3(xl,yl,zl,'g');
set(hPSC,'LineWidth',5);
Cylinder([xl(1),yl(1),zl(1)],[xl(end),yl(end),zl(end)],0.01,20,'g',1,0)

%plot axis of LSC
xl=lambda*cos(alphaLSC)*sin(betaLSC);
yl=lambda*sin(alphaLSC)*sin(betaLSC);
zl=-lambda*cos(betaLSC);
hLSC=plot3(xl,yl,zl,'r');
set(hLSC,'LineWidth',5);
Cylinder([xl(1),yl(1),zl(1)],[xl(end),yl(end),zl(end)],0.01,20,'r',1,0)

% optic axis
xl=lambda*cos(0)*sin(pi);
yl=lambda*sin(0)*sin(pi);
zl=-lambda*cos(pi);
hOpt=plot3(xl,yl,zl,'k');
set(hOpt,'LineWidth',5);
Cylinder([xl(1),yl(1),zl(1)],[xl(end),yl(end),zl(end)],0.012,20,[0.5,0.5,0.5],1,0)

% outside sphere ASC
xASCo=(R+0.005)*cos(alphaASC)*sin(betaASC);
yASCo=(R+0.005)*sin(alphaASC)*sin(betaASC);
zASCo=-(R+0.005)*cos(betaASC);
hASCo=plot3(xASCo,yASCo,zASCo,'ok','MarkerSize',14, 'LineWidth',1, 'MarkerFaceColor','b');

% outside sphere PSC
xPSCo=(R+0.005)*cos(alphaPSC)*sin(betaPSC);
yPSCo=(R+0.005)*sin(alphaPSC)*sin(betaPSC);
zPSCo=-(R+0.005)*cos(betaPSC);
hPSCo=plot3(xPSCo,yPSCo,zPSCo,'ok','MarkerSize',14, 'LineWidth',1, 'MarkerFaceColor','g');

% outside sphere LSC
xLSCo=(R+0.005)*cos(alphaLSC)*sin(betaLSC);
yLSCo=(R+0.005)*sin(alphaLSC)*sin(betaLSC);
zLSCo=-(R+0.005)*cos(betaLSC);
hLSCo=plot3(xLSCo,yLSCo,zLSCo,'ok','MarkerSize',14, 'LineWidth',1, 'MarkerFaceColor','r');

% outside sphere optic axis
xOpto=(R+0.005)*cos(0)*sin(pi);
yOpto=(R+0.005)*sin(0)*sin(pi);
zOpto=-(R+0.005)*cos(pi);
hOpto=plot3(xOpto,yOpto,zOpto,'ok','MarkerSize',16, 'LineWidth',1, 'MarkerFaceColor',[0.5,0.5,0.5]);


legend([hASC,hPSC,hLSC,hOpt],{' ASC',' PSC',' LSC',' Optic axis'},'Location','northeast','FontSize',14);
legend('boxoff');
zoom(2.2)

%% draw the margin of the retina
% hold on
% 
% betaMargin=repmat(pi/2,size(Alpha,2),1)';
% Xmargin=R*cos(Alpha).*sin(betaMargin);
% Ymargin=R*sin(Alpha).*sin(betaMargin);
% Zmargin=-R*cos(betaMargin);
% 
% plot3(Xmargin,Ymargin,Zmargin,'m--','LineWidth',4)

%hold off;

end
