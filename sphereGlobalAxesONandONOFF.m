function channels = sphereGlobalAxesONandONOFF()


%% draw visual field of both eyes
drawVisualField
f1=gcf;
f2=copyobj(f1,0);  %duplicate figure
set(0,'CurrentFigure',f1); 

%% import data
phys=dlmread('phys');
R=phys(1);

names=getFileList(cd,'singularities of global maps',0,'anywhere');
load(num2str(names{1}));

names=getFileList(cd,'singularities of global maps_ON',0,'anywhere');
ON=load(num2str(names{1}));

%% Standard Contour
% a=361;
% b=181;
% Alpha=linspace(0,2*pi,b+1);
% Alpha=Alpha(2:b+1);
% beta=linspace(0,pi,a);
% [A,B]=meshgrid(Alpha,beta);

%% Contour on Sphere
% X=R*cos(A).*sin(B);
% Y=R*sin(A).*sin(B);
% Z=-R*cos(B);
% 
% Xshift=X;
% Yshift=Y;
% Zshift=Z;
% 
% Xshift(:,b+1)=X(:,1);
% Yshift(:,b+1)=Y(:,1);
% Zshift(:,b+1)=Z(:,1);

%figure('Color',[1 1 1]);
%hold on;
% surf(Xshift,Yshift,Zshift,'EdgeColor','none','LineStyle','none','FaceLighting','phong','FaceAlpha',0.4);
% shading(h,'interp');
% colormap(h,'gray');
% axis equal
% xlabel('X');ylabel('Y');zlabel('Z');

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

%% set viewing angle
% az = -45;
% el = 225;
% view(az, el);
% set(gca, 'visible', 'off');
% set(gca,'color',[1 1 1]);
% 
%% draw the normal of ONOFF translation singularities as cylinders
lb=0;
ub=1;
lambda=linspace(lb,ub,100);

for i=1:size(Right.translation,2)
alphasing(1,i)=Right.translation{1,i}(1);       % I flipped alpha to correspond to vestibular axes
betasing(1,i)=Right.translation{1,i}(2);        
end

% transform alphasing from retinal coordinates to visual space
%allPossibleFields.alphasing=mod(360-allPossibleFields.alphasing,360);

%plot axis of nasal singularity
xT1=lambda*cosd(alphasing(1))*sind(betasing(1));
yT1=lambda*sind(alphasing(1))*sind(betasing(1));
zT1=-lambda*cosd(betasing(1));
Cylinder([xT1(1),yT1(1),zT1(1)],[xT1(end),yT1(end),zT1(end)],0.015,20,'b',1,0);       
channels.xT1=xT1; channels.yT1=yT1; channels.zT1=zT1;

%plot axis of dorsal singularity
xT2=lambda*cosd(alphasing(2))*sind(betasing(2));
yT2=lambda*sind(alphasing(2))*sind(betasing(2));
zT2=-lambda*cosd(betasing(2));
Cylinder([xT2(1),yT2(1),zT2(1)],[xT2(end),yT2(end),zT2(end)],0.015,20,'g',1,0);       
channels.xT2=xT2; channels.yT2=yT2; channels.zT2=zT2;

%plot axis of temporal
xT3=lambda*cosd(alphasing(3))*sind(betasing(3));
yT3=lambda*sind(alphasing(3))*sind(betasing(3));
zT3=-lambda*cosd(betasing(3));
Cylinder([xT3(1),yT3(1),zT3(1)],[xT3(end),yT3(end),zT3(end)],0.015,20,'r',1,0);      
channels.xT3=xT3; channels.yT3=yT3; channels.zT3=zT3;

%plot axis of ventral
xT4=lambda*cosd(alphasing(4))*sind(betasing(4));
yT4=lambda*sind(alphasing(4))*sind(betasing(4));
zT4=-lambda*cosd(betasing(4));
Cylinder([xT4(1),yT4(1),zT4(1)],[xT4(end),yT4(end),zT4(end)],0.015,20,'m',1,0);       
channels.xT4=xT4; channels.yT4=yT4; channels.zT4=zT4;



%% draw the normal of ON translation singularities as cylinders
lb=0;
ub=1;
lambda=linspace(lb,ub,100);

for i=1:size(ON.Right.translation,2)
alphasing(1,i)=ON.Right.translation{1,i}(1);       % I flipped alpha to correspond to vestibular axes
betasing(1,i)=ON.Right.translation{1,i}(2);        
end

% transform alphasing from retinal coordinates to visual space
%allPossibleFields.alphasing=mod(360-allPossibleFields.alphasing,360);

%plot axis of nasal singularity
xT1=lambda*cosd(alphasing(1))*sind(betasing(1));
yT1=lambda*sind(alphasing(1))*sind(betasing(1));
zT1=-lambda*cosd(betasing(1));
Cylinder([xT1(1),yT1(1),zT1(1)],[xT1(end),yT1(end),zT1(end)],0.015,20,'k',1,0);       
channels.xT1=xT1; channels.yT1=yT1; channels.zT1=zT1;

%plot axis of dorsal singularity
xT2=lambda*cosd(alphasing(2))*sind(betasing(2));
yT2=lambda*sind(alphasing(2))*sind(betasing(2));
zT2=-lambda*cosd(betasing(2));
Cylinder([xT2(1),yT2(1),zT2(1)],[xT2(end),yT2(end),zT2(end)],0.015,20,'k',1,0);       
channels.xT2=xT2; channels.yT2=yT2; channels.zT2=zT2;

%plot axis of temporal
xT3=lambda*cosd(alphasing(3))*sind(betasing(3));
yT3=lambda*sind(alphasing(3))*sind(betasing(3));
zT3=-lambda*cosd(betasing(3));
Cylinder([xT3(1),yT3(1),zT3(1)],[xT3(end),yT3(end),zT3(end)],0.015,20,'k',1,0);      
channels.xT3=xT3; channels.yT3=yT3; channels.zT3=zT3;

%plot axis of ventral
xT4=lambda*cosd(alphasing(4))*sind(betasing(4));
yT4=lambda*sind(alphasing(4))*sind(betasing(4));
zT4=-lambda*cosd(betasing(4));
Cylinder([xT4(1),yT4(1),zT4(1)],[xT4(end),yT4(end),zT4(end)],0.015,20,'k',1,0);       
channels.xT4=xT4; channels.yT4=yT4; channels.zT4=zT4;

xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);
axis equal

%% draw the normal of ONOFF rotation singularities as cylinders

set(0, 'CurrentFigure', f2) 

lb=0;
ub=1;
lambda=linspace(lb,ub,100);

for i=1:size(Right.rotation,2)
alphasingRot(1,i)=Right.rotation{1,i}(1);
betasingRot(1,i)=Right.rotation{1,i}(2);
end

% transform alphasing from retinal coordinates to visual space
%allPossibleFields.alphasing=mod(360-allPossibleFields.alphasing,360);

%plot axis of nasal singularity
xR1=lambda*cosd(alphasingRot(1))*sind(betasingRot(1));
yR1=lambda*sind(alphasingRot(1))*sind(betasingRot(1));
zR1=-lambda*cosd(betasingRot(1));
cl=Cylinder([xR1(1),yR1(1),zR1(1)],[xR1(end),yR1(end),zR1(end)],0.015,20,'b',1,0);        
set(cl,'FaceAlpha',0.3);
channels.xR1=xR1; channels.yR1=yR1; channels.zR1=zR1;

%plot axis of dorsal singularity
xR2=lambda*cosd(alphasingRot(2))*sind(betasingRot(2));
yR2=lambda*sind(alphasingRot(2))*sind(betasingRot(2));
zR2=-lambda*cosd(betasingRot(2));
cl=Cylinder([xR2(1),yR2(1),zR2(1)],[xR2(end),yR2(end),zR2(end)],0.015,20,'g',1,0);        
set(cl,'FaceAlpha',0.3);
channels.xR2=xR2; channels.yR2=yR2; channels.zR2=zR2;

%plot axis of temporal
xR3=lambda*cosd(alphasingRot(3))*sind(betasingRot(3));
yR3=lambda*sind(alphasingRot(3))*sind(betasingRot(3));
zR3=-lambda*cosd(betasingRot(3));
cl=Cylinder([xR3(1),yR3(1),zR3(1)],[xR3(end),yR3(end),zR3(end)],0.015,20,'r',1,0);        
set(cl,'FaceAlpha',0.3);
channels.xR3=xR3; channels.yR3=yR3; channels.zR3=zR3;

%plot axis of ventral
xR4=lambda*cosd(alphasingRot(4))*sind(betasingRot(4));
yR4=lambda*sind(alphasingRot(4))*sind(betasingRot(4));
zR4=-lambda*cosd(betasingRot(4));
cl=Cylinder([xR4(1),yR4(1),zR4(1)],[xR4(end),yR4(end),zR4(end)],0.015,20,'m',1,0);     
set(cl,'FaceAlpha',0.3);
channels.xR4=xR4; channels.yR4=yR4; channels.zR4=zR4;


%% draw the normal of ON rotation singularities as cylinders

set(0, 'CurrentFigure', f2) 

lb=0;
ub=1;
lambda=linspace(lb,ub,100);

for i=1:size(ON.Right.rotation,2)
alphasingRot(1,i)=ON.Right.rotation{1,i}(1);
betasingRot(1,i)=ON.Right.rotation{1,i}(2);
end

% transform alphasing from retinal coordinates to visual space
%allPossibleFields.alphasing=mod(360-allPossibleFields.alphasing,360);

%plot axis of nasal singularity
xR1=lambda*cosd(alphasingRot(1))*sind(betasingRot(1));
yR1=lambda*sind(alphasingRot(1))*sind(betasingRot(1));
zR1=-lambda*cosd(betasingRot(1));
cl=Cylinder([xR1(1),yR1(1),zR1(1)],[xR1(end),yR1(end),zR1(end)],0.015,20,'k',1,0);        
set(cl,'FaceAlpha',0.3);
channels.xR1=xR1; channels.yR1=yR1; channels.zR1=zR1;

%plot axis of dorsal singularity
xR2=lambda*cosd(alphasingRot(2))*sind(betasingRot(2));
yR2=lambda*sind(alphasingRot(2))*sind(betasingRot(2));
zR2=-lambda*cosd(betasingRot(2));
cl=Cylinder([xR2(1),yR2(1),zR2(1)],[xR2(end),yR2(end),zR2(end)],0.015,20,'k',1,0);        
set(cl,'FaceAlpha',0.3);
channels.xR2=xR2; channels.yR2=yR2; channels.zR2=zR2;

%plot axis of temporal
xR3=lambda*cosd(alphasingRot(3))*sind(betasingRot(3));
yR3=lambda*sind(alphasingRot(3))*sind(betasingRot(3));
zR3=-lambda*cosd(betasingRot(3));
cl=Cylinder([xR3(1),yR3(1),zR3(1)],[xR3(end),yR3(end),zR3(end)],0.015,20,'k',1,0);        
set(cl,'FaceAlpha',0.3);
channels.xR3=xR3; channels.yR3=yR3; channels.zR3=zR3;

%plot axis of ventral
xR4=lambda*cosd(alphasingRot(4))*sind(betasingRot(4));
yR4=lambda*sind(alphasingRot(4))*sind(betasingRot(4));
zR4=-lambda*cosd(betasingRot(4));
cl=Cylinder([xR4(1),yR4(1),zR4(1)],[xR4(end),yR4(end),zR4(end)],0.015,20,'k',1,0);     
set(cl,'FaceAlpha',0.3);
channels.xR4=xR4; channels.yR4=yR4; channels.zR4=zR4;

xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);
axis equal

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
%hOpt=plot3(0,0,0.5,'ok','MarkerSize',10, 'LineWidth',1','MarkerFaceColor','k');

end
