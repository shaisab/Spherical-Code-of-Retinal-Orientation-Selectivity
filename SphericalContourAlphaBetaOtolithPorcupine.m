function []=SphericalContourAlphaBetaOtolithPorcupine(resolution)

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
% sfine=linspace(0,pi,20);
% 
% [THfine,Sfine]=meshgrid(thetafine,sfine);
% 
% Xfine=(R+0.01)*cos(THfine).*sin(Sfine);
% Yfine=(R+0.01)*sin(THfine).*sin(Sfine);
% Zfine=-(R+0.01)*cos(Sfine);
% 
% surf(Xfine,Yfine,Zfine,'FaceColor','none','EdgeAlpha',0.1);


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


%% draw the margin of the retina
% hold on
% 
% betaMargin=repmat(pi/2,size(Alpha,2),1)';
% Xmargin=R*cos(Alpha).*sin(betaMargin);
% Ymargin=R*sin(Alpha).*sin(betaMargin);
% Zmargin=-R*cos(betaMargin);
% 
% plot3(Xmargin,Ymargin,Zmargin,'Color',[0.5,0.5,0.5],'LineStyle','--','LineWidth',3)

%% plot symbols for optic axis
hOpt=plot3(0,0,0.5,'ok','MarkerSize',10, 'LineWidth',1','MarkerFaceColor','k');

%% set viewing angle
az = -45;
el = 225;
view(az, el);
set(gca, 'visible', 'off');
set(gca,'color',[1 1 1]);

%% draw optic axis as a line
lb=-0.75;
ub=0.75;

% lb=-1;
% ub=1;

lambda=linspace(lb,ub,100);

% optic axis
xl=lambda*cos(0)*sin(pi);
yl=lambda*sin(0)*sin(pi);
zl=-lambda*cos(pi);
h=plot3(xl,yl,zl,'Color',[0.5,0.5,0.5]);
set(h,'LineWidth',5);
Cylinder([xl(1),yl(1),zl(1)],[xl(end),yl(end),zl(end)],0.02,20,[0.5,0.5,0.5],1,0);

% outside sphere optic axis
xOpto=(R+0.005)*cos(0)*sin(pi);
yOpto=(R+0.005)*sin(0)*sin(pi);
zOpto=-(R+0.005)*cos(pi);
hOpto=plot3(xOpto,yOpto,zOpto,'ok','MarkerSize',8, 'LineWidth',1, 'MarkerFaceColor',[0.5,0.5,0.5]);

%% draw the normal of singularities as lines
%lb=-0.75;
lb=0;
ubNasal=2.4;          %3
ubDorsal=3;       %3      3 for nonCART
ubTemporal=1;     %1.1
ubVentral=1.7;      %1.7    3 for Hb9

% lb=-1;
% ub=1;

lambdaNasal=linspace(lb,ubNasal,100);
lambdaDorsal=linspace(lb,ubDorsal,100);
lambdaTemporal=linspace(lb,ubTemporal,100);
lambdaVentral=linspace(lb,ubVentral,100);

% transform alphasing from retinal coordinates to visual space
allPossibleFields.alphasing=mod(360-allPossibleFields.alphasing,360);


%plot axis of nasal singularity
xl=lambdaNasal*cosd(allPossibleFields.alphasing(1))*sind(allPossibleFields.betasing(1));
yl=lambdaNasal*sind(allPossibleFields.alphasing(1))*sind(allPossibleFields.betasing(1));
zl=-lambdaNasal*cosd(allPossibleFields.betasing(1));
hnasal=plot3(xl,yl,zl,'b');
set(hnasal,'LineWidth',5);
Cylinder([xl(1),yl(1),zl(1)],[xl(end),yl(end),zl(end)],0.01,20,'b',1,0);

%plot axis of dorsal singularity
xl=lambdaDorsal*cosd(allPossibleFields.alphasing(2))*sind(allPossibleFields.betasing(2));
yl=lambdaDorsal*sind(allPossibleFields.alphasing(2))*sind(allPossibleFields.betasing(2));
zl=-lambdaDorsal*cosd(allPossibleFields.betasing(2));
hdorsal=plot3(xl,yl,zl,'g');
set(hdorsal,'LineWidth',5);
Cylinder([xl(1),yl(1),zl(1)],[xl(end),yl(end),zl(end)],0.01,20,'g',1,0);

%plot axis of temporal
xl=lambdaTemporal*cosd(allPossibleFields.alphasing(3))*sind(allPossibleFields.betasing(3));
yl=lambdaTemporal*sind(allPossibleFields.alphasing(3))*sind(allPossibleFields.betasing(3));
zl=-lambdaTemporal*cosd(allPossibleFields.betasing(3));
htemporal=plot3(xl,yl,zl,'r');
set(htemporal,'LineWidth',5);
Cylinder([xl(1),yl(1),zl(1)],[xl(end),yl(end),zl(end)],0.01,20,'r',1,0);

%plot axis of ventral
xl=lambdaVentral*cosd(allPossibleFields.alphasing(4))*sind(allPossibleFields.betasing(4));
yl=lambdaVentral*sind(allPossibleFields.alphasing(4))*sind(allPossibleFields.betasing(4));
zl=-lambdaVentral*cosd(allPossibleFields.betasing(4));
hventral=plot3(xl,yl,zl,'m');
set(hventral,'LineWidth',5);
Cylinder([xl(1),yl(1),zl(1)],[xl(end),yl(end),zl(end)],0.005,20,'m',1,0)

% % optic axis
% xl=lambda*cos(0)*sin(pi);
% yl=lambda*sin(0)*sin(pi);
% zl=-lambda*cos(pi);
% hOpt=plot3(xl,yl,zl,'k');
% set(hOpt,'LineWidth',5);

% outside sphere nasal singularity
xnasalo=(R+0.005)*cosd(allPossibleFields.alphasing(1))*sind(allPossibleFields.betasing(1));
ynasalo=(R+0.005)*sind(allPossibleFields.alphasing(1))*sind(allPossibleFields.betasing(1));
znasalo=-(R+0.005)*cosd(allPossibleFields.betasing(1));
hnasalo=plot3(xnasalo,ynasalo,znasalo,'ok','MarkerSize',12, 'LineWidth',1, 'MarkerFaceColor','b');

% outside sphere dorsal singularity
xdorsalo=(R+0.005)*cosd(allPossibleFields.alphasing(2))*sind(allPossibleFields.betasing(2));
ydorsalo=(R+0.005)*sind(allPossibleFields.alphasing(2))*sind(allPossibleFields.betasing(2));
zdorsalo=-(R+0.005)*cosd(allPossibleFields.betasing(2));
hdorsalo=plot3(xdorsalo,ydorsalo,zdorsalo,'ok','MarkerSize',12, 'LineWidth',1, 'MarkerFaceColor','g');

% outside sphere temporal singularity
xtemporalo=(R+0.005)*cosd(allPossibleFields.alphasing(3))*sind(allPossibleFields.betasing(3));
ytemporalo=(R+0.005)*sind(allPossibleFields.alphasing(3))*sind(allPossibleFields.betasing(3));
ztemporalo=-(R+0.005)*cosd(allPossibleFields.betasing(3));
htemporalo=plot3(xtemporalo,ytemporalo,ztemporalo,'ok','MarkerSize',12, 'LineWidth',1, 'MarkerFaceColor','r');

% outside sphere ventral singularity
xventralo=(R+0.005)*cosd(allPossibleFields.alphasing(4))*sind(allPossibleFields.betasing(4));
yventralo=(R+0.005)*sind(allPossibleFields.alphasing(4))*sind(allPossibleFields.betasing(4));
zventralo=-(R+0.005)*cosd(allPossibleFields.betasing(4));
hventralo=plot3(xventralo,yventralo,zventralo,'ok','MarkerSize',8, 'LineWidth',1, 'MarkerFaceColor','m');



%% plot normal heat vectors
%load('allPossibleFieldsAlphaCorr_Calcium imaging data_real_translation_40_DS_created_Oct_09_2015_13_13_ONOFF.mat');
resolution=floor(370/size(allPossibleFields.matchind180,2));
c=10;
cutoff=(['c',num2str(c)]);
k=1;
h=1;
for alpha=0:resolution:360
for beta=0:resolution:180
    data(k,h)=allPossibleFields.matchind180(k,h).(cutoff);
    k=k+1;
end
k=1;
h=h+1;
end

data=fliplr(data);      % flip data along the x axis of alpha-beta map

% plot surface
%figure;
alpha=0:resolution:360;
beta=0:resolution:180;
%contourf(alpha,beta,data);
%hold on;

% down sampling by a factor of 2 in order to see vectors bettern
downFactor=2;
alphaDown=downsample(alpha,downFactor);
betaDown=downsample(beta,downFactor);
dataDown=downsample(data,downFactor);
dataDown=downsample(dataDown',downFactor);
dataDown=dataDown';
zeroheight=zeros(size(dataDown,1),size(dataDown,2));

% calculating normals
[X,Y]=meshgrid(0:resolution*downFactor:360,0:resolution*downFactor:180);
[xn,yn,zn]=surfnorm(X,Y,dataDown);

% for Hb9 physiology
%dataDown=(dataDown.^2)./40;

% for ONOFF
dataDown=(dataDown.^2)./40;

%colorVector=real2rgb(dataDown,'jet');

tmp=rand(size(dataDown,1),3)
% plot normals
%quiver3(alphaDown,betaDown,zeroheight,xn,yn,dataDown,'AutoScale','off','Color','r','ShowArrowHead','off','LineWidth',1.5);
%h.ShowArrowHead='off';

%figure;
% plot normals
hold on;
[alphaDown,betaDown]=meshgrid(alphaDown,betaDown);
dataDown=dataDown/20;       % /30 for ONOFF and Hb9
Vx=cosd(alphaDown).*sind(betaDown).*dataDown;
Vy=sind(alphaDown).*sind(betaDown).*dataDown;
Vz=-cosd(betaDown).*dataDown;

surf(R*cosd(alphaDown).*sind(betaDown),R*sind(alphaDown).*sind(betaDown),-R*cosd(betaDown),'FaceColor','none','EdgeAlpha',0.2);
for i=1:size(alphaDown,1)
    for j=1:size(alphaDown,2)
%h=quiver3(R*cosd(alphaDown(i,j)).*sind(betaDown(i,j)),R*sind(alphaDown(i,j)).*sind(betaDown(i,j)),-R*cosd(betaDown(i,j)),Vx(i,j),Vy(i,j),Vz(i,j),'AutoScale','off','Color',squeeze(colorVector(i,j,:))','LineWidth',dataDown(i,j).*3);

% for ONOFF
%h=quiver3(R*cosd(alphaDown(i,j)).*sind(betaDown(i,j)),R*sind(alphaDown(i,j)).*sind(betaDown(i,j)),-R*cosd(betaDown(i,j)),Vx(i,j),Vy(i,j),Vz(i,j),'AutoScale','off','Color','k','LineWidth',(dataDown(i,j)+0.001)*2.5);

% for Hb9 physiology
%h=quiver3(R*cosd(alphaDown(i,j)).*sind(betaDown(i,j)),R*sind(alphaDown(i,j)).*sind(betaDown(i,j)),-R*cosd(betaDown(i,j)),Vx(i,j),Vy(i,j),Vz(i,j),'AutoScale','off','Color','k','LineWidth',(dataDown(i,j)+0.001)*0.8);

% for Trhr physiology
h=quiver3(R*cosd(alphaDown(i,j)).*sind(betaDown(i,j)),R*sind(alphaDown(i,j)).*sind(betaDown(i,j)),-R*cosd(betaDown(i,j)),Vx(i,j),Vy(i,j),Vz(i,j),'AutoScale','off','Color','k','LineWidth',(dataDown(i,j)+0.001)*0.8);


h.ShowArrowHead='off';
    end
end


axis equal;
axis off
end
