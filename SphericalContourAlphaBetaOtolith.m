function []=SphericalContourAlphaBetaOtolith(resolution)

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
            
            data=fliplr(data);      % flip data along the x axis of alpha-beta map
            
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
thetafine=linspace(0,2*pi,10);
sfine=linspace(0,pi,5);

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

%% draw intersection lines of saccule
hold on

n=1000;
normal=[0.795000000000000;0.544000000000000;0.254000000000000];       % normal vector of saccule
[alphaCirc,betaCirc]=GreatCircle(normal,n);
Xsaccule=1.01*R*cos(alphaCirc).*sin(betaCirc);
Ysaccule=1.01*R*sin(alphaCirc).*sin(betaCirc);
Zsaccule=-1.01*R*cos(betaCirc);

plot3(Xsaccule,Ysaccule,Zsaccule,'k','LineWidth',3)

hSaccule=plot3(Xsaccule,Ysaccule,Zsaccule,'k','LineWidth',5);

%% draw intersection lines of utricle
n=1000;
normal=[-0.447172870422993;0.229457818946242;0.854651702907152];       % normal vector of utricle
[alphaCirc,betaCirc]=GreatCircle(normal,n);
Xutricle=1.01*R*cos(alphaCirc).*sin(betaCirc);
Yutricle=1.01*R*sin(alphaCirc).*sin(betaCirc);
Zutricle=-1.01*R*cos(betaCirc);

hUtricle=plot3(Xutricle,Yutricle,Zutricle,'r','LineWidth',5);

%% draw the margin of the retina
% hold on
% 
% betaMargin=repmat(pi/2,size(Alpha,2),1)';
% Xmargin=R*cos(Alpha).*sin(betaMargin);
% Ymargin=R*sin(Alpha).*sin(betaMargin);
% Zmargin=-R*cos(betaMargin);
% 
% plot3(Xmargin,Ymargin,Zmargin,'m--','LineWidth',4)

%% plot symbols for optic axis
hOpt=plot3(0,0,0.5,'ok','MarkerSize',10, 'LineWidth',1','MarkerFaceColor',[0.5,0.5,0.5]);

%% set viewing angle
az = -45;
el = 225;
view(az, el);
set(gca, 'visible', 'off');
set(gca,'color',[1 1 1]);

%% draw the normal of optic axis lines
lb=-0.75;
ub=0.75;

% lb=-1;
% ub=1;

lambda=linspace(lb,ub,100);

% optic axis
xl=lambda*cos(0)*sin(pi);
yl=lambda*sin(0)*sin(pi);
zl=-lambda*cos(pi);
% hOpt=plot3(xl,yl,zl,'Color',[0.5,0.5,0.5]);
% set(hOpt,'LineWidth',7);

Cylinder([xl(1),yl(1),zl(1)],[xl(end),yl(end),zl(end)],0.02,20,[0.5,0.5,0.5],1,0)

% outside sphere optic axis
xOpto=(R+0.005)*cos(0)*sin(pi);
yOpto=(R+0.005)*sin(0)*sin(pi);
zOpto=-(R+0.005)*cos(pi);
hOpto=plot3(xOpto,yOpto,zOpto,'ok','MarkerSize',20, 'LineWidth',1, 'MarkerFaceColor',[0.5,0.5,0.5]);

legend([hSaccule,hUtricle,hOpt],{' Saccular macula',' Utricular macula',' Optic axis'},'Location','northeast','FontSize',14);
legend('boxoff');
zoom(2)


end
