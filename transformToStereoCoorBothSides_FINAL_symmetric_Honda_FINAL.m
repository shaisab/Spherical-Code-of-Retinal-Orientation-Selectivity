
%clear all;

%% Right sacular macula for 4 mice SMi-SMs
meanV=[-0.119,0.951,-0.286];
RightSM1=-meanV;
hold on
%% Right sacular macula for 4 mice SMf-SMi
meanV=[0.246,0.675,-0.695];
RightSM2=-meanV;
xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);

%% Right sacular macula for 4 mice SMs-SMf
meanV=[-0.381,0.917,0.121];
RightSM3=-meanV;
xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);

%% Left sacular macula for 4 mice SMi-SMs
meanV=[-0.119,0.951,0.286];
LeftSM1=-meanV;

%% Left sacular macula for 4 mice SMf-SMi
meanV=[0.246,0.675,0.695];
LeftSM2=-meanV;
xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);

%% Left sacular macula for 4 mice SMs-SMf
meanV=[-0.381,0.917,-0.121];
LeftSM3=-meanV;
xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);

%% Lambda for 4 mice
meanV=[0.094791309,-0.102330783,9.42255E-18];
Lambda=meanV;
xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);

%% input mean vectors and coordinates from Eye-Canal-Vectors-21-11-14

LB=-[1,0,0];     % Lambda - Bregma mean vector

RightoptN=[-0.237,-0.490,0.839];    % right optic axis vector from Nathan
RightMR=[-0.612,-0.577,-0.540];    % medial rectus mean vector right side (flipped to make the vector to go from origin to insertion)
RightLR=[-0.230,-0.461,-0.857];    % lateral rectus mean vector right side (flipped to make the vector to go from origin to insertion)
RightASC=-[0.604,0.243,0.758];      % mean vector right side
RightPSC=-[0.447,0.469,-0.760];     % mean vector right side
RightLSC=[0.449,-0.886,-0.085];    % mean vector right side
RightSM=[-0.795,-0.254,-0.544];           % mean vector right side

RightSMh=[-0.469,-0.252,-0.825];           % mean vector right side
RightUT=[0.584,-0.798,-0.072];             % mean vector right side

% find the rotation matrix needed to transform the SM from Honda to SM from Nathan, and use it to rotate UT of Honda vectors  
R=vrrotvec2mat(vrrotvec(RightSMh',RightSM'));
RightSMhr=R*RightSMh';
RightUTr=R*RightUT';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LeftoptN=[-0.237,-0.490,-0.839];      % left optic axis vector from Nathan
LeftMR=[-0.612,-0.577,0.540];    % medial rectus mean vector left side (flipped to make the vector to go from origin to insertion)
LeftLR=[-0.230,-0.461,0.857];    % lateral rectus mean vector left side (flipped to make the vector to go from origin to insertion)
LeftASC=[-0.604,-0.243,0.758];     % mean vector left side
LeftPSC=[-0.447,-0.469,-0.760];     % mean vector left side
LeftLSC=-[-0.449,0.886,-0.085];     % mean vector left side
LeftSM=-[0.795,0.254,-0.544];
% LeftSMh=[-0.523,-0.211,0.816];           % mean vector right side
% LeftUT=[0.579,-0.807,0.006];             % mean vector right side

LeftSMh=[RightSMhr(1),RightSMhr(2),-RightSMhr(3)];           % mean vector right side - standardized for our data and then y and z switched and the polarity inverted
LeftUT=[RightUTr(1),RightUTr(2),-RightUTr(3)];           % mean vector right side - standardized for our data and then y and z switched and the polarity inverted

%normal=[0,-0.998571,0.0534409];        % this normal was calculated by the mathematica code
normal=[0,-1,0];                        % minizing the roll was already done by Nathan

%% find the rotation matrix needed to transform the normal caculated in Mathematica to the vertical axis
% [0,0,1], and use it to rotate all vectors  
R=vrrotvec2mat(vrrotvec(normal',[0,0,1]'));

LBr=R*LB';
RightoptNr=R*RightoptN';
RightMRr=R*RightMR';
RightLRr=R*RightLR';
RightASCr=R*RightASC';
RightPSCr=R*RightPSC';
RightLSCr=R*RightLSC';
RightSM1r=R*RightSM1';
RightSM2r=R*RightSM2';
RightSM3r=R*RightSM3';
RightSMr=R*RightSM';
RightSMhr=R*RightSMhr;
RightUTr=R*RightUTr;

LeftoptNr=R*LeftoptN';
LeftMRr=R*LeftMR';
LeftLRr=R*LeftLR';
LeftASCr=R*LeftASC';
LeftPSCr=R*LeftPSC';
LeftLSCr=R*LeftLSC';
LeftSM1r=R*LeftSM1';
LeftSM2r=R*LeftSM2';
LeftSM3r=R*LeftSM3';
LeftSMr=R*LeftSM';
LeftSMhr=R*LeftSMh';
LeftUTr=R*LeftUT';

%% find the rotation matrix needed to transform LB to the midsagital axis
% [1,0,0], and use it to rotate all vectors 
R=vrrotvec2mat(vrrotvec(LBr',[1,0,0]'));

LBr=R*LBr;

RightoptNr=R*RightoptNr;
RightMRr=R*RightMRr;
RightLRr=R*RightLRr;
RightASCr=R*RightASCr;
RightPSCr=R*RightPSCr;
RightLSCr=R*RightLSCr;
RightSM1r=R*RightSM1r;
RightSM2r=R*RightSM2r;
RightSM3r=R*RightSM3r;
RightSMr=R*RightSMr;
RightSMhr=R*RightSMhr;
RightUTr=R*RightUTr;
%RightLSCampPSCintr=R*RightLSCampPSCintr;

LeftoptNr=R*LeftoptNr;
LeftMRr=R*LeftMRr;
LeftLRr=R*LeftLRr;
LeftASCr=R*LeftASCr;
LeftPSCr=R*LeftPSCr;
LeftLSCr=R*LeftLSCr;
LeftSM1r=R*LeftSM1r;
LeftSM2r=R*LeftSM2r;
LeftSM3r=R*LeftSM3r;
LeftSMr=R*LeftSMr;
LeftSMhr=R*LeftSMhr;
LeftUTr=R*LeftUTr;

%% plot all vectors

f1=figure;
quiver3(0,0,0,LBr(1),LBr(2),LBr(3),'m','LineWidth',4);
hold on;

%quiver3(0,0,0,1*RightoptNr(1),1*RightoptNr(2),1*RightoptNr(3),'k--','LineWidth',4);
% quiver3(0,0,0,1*RightMRr(1),1*RightMRr(2),1*RightMRr(3),'c','LineWidth',2);
% quiver3(0,0,0,1*RightLRr(1),1*RightLRr(2),1*RightLRr(3),'c','LineWidth',2);
quiver3(0,0,0,RightASCr(1),RightASCr(2),RightASCr(3),'b','LineWidth',2);
quiver3(0,0,0,RightPSCr(1),RightPSCr(2),RightPSCr(3),'g','LineWidth',2);
quiver3(0,0,0,RightLSCr(1),RightLSCr(2),RightLSCr(3),'r','LineWidth',2);
% quiver3(0,0,0,RightSM1r(1),RightSM1r(2),RightSM1r(3),'c','LineWidth',2);
% quiver3(0,0,0,RightSM2r(1),RightSM2r(2),RightSM2r(3),'c','LineWidth',2);
% quiver3(0,0,0,RightSM3r(1),RightSM3r(2),RightSM3r(3),'c','LineWidth',2);
quiver3(0,0,0,RightSMr(1),RightSMr(2),RightSMr(3),'y','LineWidth',2);
%quiver3(0,0,0,RightSMhr(1),RightSMhr(2),RightSMhr(3),'c','LineWidth',2);
quiver3(0,0,0,RightUTr(1),RightUTr(2),RightUTr(3),'Color',[0 .7 .7],'LineStyle','-','LineWidth',2);
%quiver3(0,0,0,RightLSCampPSCintr(1),RightLSCampPSCintr(2),RightLSCampPSCintr(3),'Color',[0.2,0.4,0.6],'LineWidth',2);

%quiver3(0,0,0,1*LeftoptNr(1),1*LeftoptNr(2),1*LeftoptNr(3),'k--','LineWidth',2);
% quiver3(0,0,0,1*LeftMRr(1),1*LeftMRr(2),1*LeftMRr(3),'c','LineWidth',2);
% quiver3(0,0,0,1*LeftLRr(1),1*LeftLRr(2),1*LeftLRr(3),'c','LineWidth',2);
quiver3(0,0,0,LeftASCr(1),LeftASCr(2),LeftASCr(3),'b--','LineWidth',2);
quiver3(0,0,0,LeftPSCr(1),LeftPSCr(2),LeftPSCr(3),'g--','LineWidth',2);
quiver3(0,0,0,LeftLSCr(1),LeftLSCr(2),LeftLSCr(3),'r--','LineWidth',2);
%quiver3(0,0,0,LeftSM1r(1),LeftSM1r(2),LeftSM1r(3),'c--','LineWidth',2);
%quiver3(0,0,0,LeftSM2r(1),LeftSM2r(2),LeftSM2r(3),'c--','LineWidth',2);
%quiver3(0,0,0,LeftSM3r(1),LeftSM3r(2),LeftSM3r(3),'c--','LineWidth',2);
quiver3(0,0,0,LeftSMr(1),LeftSMr(2),LeftSMr(3),'y--','LineWidth',2);
%quiver3(0,0,0,LeftSMhr(1),LeftSMhr(2),LeftSMhr(3),'c--','LineWidth',2);
quiver3(0,0,0,LeftUTr(1),LeftUTr(2),LeftUTr(3),'Color',[0 .7 .7],'LineStyle','--','LineWidth',2);

set(gcf, 'color', [1 1 1]);
xlabel('X');ylabel('Y');zlabel('Z');
xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);

% Now LB correspond to the x-axis and the horizon.

%% rotate all vectors by 29 deg head down around an axis perpendicualr to lambda-bregma
tmp=rodrigues_rot(LBr,[0;1;0],deg2rad(29));

%% find the rotation matrix needed to transform LB to tmp
% [1,0,0], and use it to rotate all vectors 
R=vrrotvec2mat(vrrotvec(LBr',tmp'));

LBr=R*LBr;

RightoptNr=R*RightoptNr;
RightMRr=R*RightMRr;
RightLRr=R*RightLRr;
RightASCr=R*RightASCr;
RightPSCr=R*RightPSCr;
RightLSCr=R*RightLSCr;
RightSM1r=R*RightSM1r;
RightSM2r=R*RightSM2r;
RightSM3r=R*RightSM3r;
RightSMr=R*RightSMr;
RightSMhr=R*RightSMhr;
RightUTr=R*RightUTr;
%RightLSCampPSCintr=R*RightLSCampPSCintr;

LeftoptNr=R*LeftoptNr;
LeftMRr=R*LeftMRr;
LeftLRr=R*LeftLRr;
LeftASCr=R*LeftASCr;
LeftPSCr=R*LeftPSCr;
LeftLSCr=R*LeftLSCr;
LeftSM1r=R*LeftSM1r;
LeftSM2r=R*LeftSM2r;
LeftSM3r=R*LeftSM3r;
LeftSMr=R*LeftSMr;
LeftSMhr=R*LeftSMhr;
LeftUTr=R*LeftUTr;


%% plot all vectors after applying tilt down by 29 deg

f2=figure;
quiver3(0,0,0,LBr(1),LBr(2),LBr(3),'m','LineWidth',4);
hold on;

%quiver3(0,0,0,1*RightoptNr(1),1*RightoptNr(2),1*RightoptNr(3),'k--','LineWidth',4);
% quiver3(0,0,0,1*RightMRr(1),1*RightMRr(2),1*RightMRr(3),'c','LineWidth',2);
% quiver3(0,0,0,1*RightLRr(1),1*RightLRr(2),1*RightLRr(3),'c','LineWidth',2);
quiver3(0,0,0,RightASCr(1),RightASCr(2),RightASCr(3),'b','LineWidth',2);
quiver3(0,0,0,RightPSCr(1),RightPSCr(2),RightPSCr(3),'g','LineWidth',2);
quiver3(0,0,0,RightLSCr(1),RightLSCr(2),RightLSCr(3),'r','LineWidth',2);
% quiver3(0,0,0,RightSM1r(1),RightSM1r(2),RightSM1r(3),'c','LineWidth',2);
% quiver3(0,0,0,RightSM2r(1),RightSM2r(2),RightSM2r(3),'c','LineWidth',2);
% quiver3(0,0,0,RightSM3r(1),RightSM3r(2),RightSM3r(3),'c','LineWidth',2);
quiver3(0,0,0,RightSMr(1),RightSMr(2),RightSMr(3),'y','LineWidth',2);
%quiver3(0,0,0,RightSMhr(1),RightSMhr(2),RightSMhr(3),'c','LineWidth',2);
quiver3(0,0,0,RightUTr(1),RightUTr(2),RightUTr(3),'Color',[0 .7 .7],'LineStyle','-','LineWidth',2);
%quiver3(0,0,0,RightLSCampPSCintr(1),RightLSCampPSCintr(2),RightLSCampPSCintr(3),'Color',[0.2,0.4,0.6],'LineWidth',2);

%quiver3(0,0,0,1*LeftoptNr(1),1*LeftoptNr(2),1*LeftoptNr(3),'k--','LineWidth',2);
% quiver3(0,0,0,1*LeftMRr(1),1*LeftMRr(2),1*LeftMRr(3),'c','LineWidth',2);
% quiver3(0,0,0,1*LeftLRr(1),1*LeftLRr(2),1*LeftLRr(3),'c','LineWidth',2);
quiver3(0,0,0,LeftASCr(1),LeftASCr(2),LeftASCr(3),'b--','LineWidth',2);
quiver3(0,0,0,LeftPSCr(1),LeftPSCr(2),LeftPSCr(3),'g--','LineWidth',2);
quiver3(0,0,0,LeftLSCr(1),LeftLSCr(2),LeftLSCr(3),'r--','LineWidth',2);
%quiver3(0,0,0,LeftSM1r(1),LeftSM1r(2),LeftSM1r(3),'c--','LineWidth',2);
%quiver3(0,0,0,LeftSM2r(1),LeftSM2r(2),LeftSM2r(3),'c--','LineWidth',2);
%quiver3(0,0,0,LeftSM3r(1),LeftSM3r(2),LeftSM3r(3),'c--','LineWidth',2);
quiver3(0,0,0,LeftSMr(1),LeftSMr(2),LeftSMr(3),'y--','LineWidth',2);
%quiver3(0,0,0,LeftSMhr(1),LeftSMhr(2),LeftSMhr(3),'c--','LineWidth',2);
quiver3(0,0,0,LeftUTr(1),LeftUTr(2),LeftUTr(3),'Color',[0 .7 .7],'LineStyle','--','LineWidth',2);

set(gcf, 'color', [1 1 1]);
xlabel('X');ylabel('Y');zlabel('Z');
xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);

% Now LB is tilted 29 down from the horizon. Now all vectors are in
% relative to the horizon and can be combined with the axes of the visual
% system.

%% caculate and plot the vector of the right optic disc based on the litrature (This optic axis is relative to the horizon and stright ahead)
RightoptL=[sind(90-22)*cosd(-64),sind(90-22)*sind(-64),cosd(90-22)];      %from Oommen and Stahl, 2008 - "during ambulation, each eye
                                                                          %is deviated 64° lateral to the mid-sagittal plane and 22° superior to the
                                                                          %horizon" and "during ambulation, mice hold the head with
                                                                          %the lambda–bregma (L–B) axis inclined 29° down"

quiver3(0,0,0,1*RightoptL(1),1*RightoptL(2),1*RightoptL(3),'k','LineWidth',3);

%% caculate and plot the vector of the left optic disc based on the litrature (This optic axis is relative to the horizon and stright ahead)
LeftoptL=[sind(90-22)*cosd(64),sind(90-22)*sind(64),cosd(90-22)];      %from Oommen and Stahl, 2008 - "during ambulation, each eye
                                                                          %is deviated 64° lateral to the mid-sagittal plane and 22° superior to the
                                                                          %horizon" and "during ambulation, mice hold the head with
                                                                          %the lambda–bregma (L–B) axis inclined 29° down"

quiver3(0,0,0,1*LeftoptL(1),1*LeftoptL(2),1*LeftoptL(3),'k--','LineWidth',3);


%% draw horizontal and midsagittal planes
[X,Y]=meshgrid(linspace(-1,1,50),linspace(-1,1,50));
Z1=-(0*X+0*Y)./1;
surf(X,Y,Z1,'FaceColor','k','EdgeColor','none','FaceAlpha',0.8);  % horizontal plane
Z2=-(0*X+(-1)*Y)./1E-10;
surf(X,Y,Z2,'FaceColor','k','EdgeColor','none','FaceAlpha',0.8);    % plane of saccula
plot3([-1;1],[0;0],[0;0],'k','LineWidth',1);
xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);

%% draw planes of utricle and saccula
[X,Y]=meshgrid(linspace(-1,1,50),linspace(-1,1,50));
Z1=-(RightUTr(1)*X+RightUTr(2)*Y)./RightUTr(3);
surf(X,Y,Z1,'FaceColor',[0 .7 .7],'EdgeColor','none','FaceAlpha',0.8);  % plane of utricle
Z2=-(RightSMr(1)*X+RightSMr(2)*Y)./RightSMr(3);
surf(X,Y,Z2,'FaceColor','m','EdgeColor','none','FaceAlpha',0.8);    % plane of saccula

%% plot unit sphere
% [x,y,z] = sphere(100);
% x=x.*0.92;
% y=y.*0.92;
% z=z.*0.92;
% hsurf=surf(x,y,z,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
% axis equal
% colormap cool
% shading interp;
% alpha(0.3)
% lightangle(0,180)
% hsurf.FaceLighting = 'gouraud';
% hsurf.AmbientStrength = 0.7;
% hsurf.DiffuseStrength = 0.1;
% hsurf.SpecularStrength = 0.5;
% hsurf.SpecularExponent = 25;
% hsurf.BackFaceLighting = 'lit';
% 
% %% plot equator and prime meridian and x,y,z axes
% aziEquator=1:1:360;
% elevEquator=zeros(1,360);
% [x,y,z]=sph2cart(aziEquator,elevEquator,0.9);
% scatter3(x,y,z,10,'k.');
% 
% aziPM=zeros(1,360);
% elevPM=1:1:360;
% [xPM,yPM,zPM]=sph2cart(aziPM,elevPM,0.9);
% scatter3(xPM,yPM,zPM,10,'k.');
% 
% hold on
% factor2=0.9;
% plot3(factor2.*[-1,1],factor2.*[0,0],factor2.*[0,0],'k:','LineWidth',1);
% plot3(factor2.*[0,0],factor2.*[-1,1],factor2.*[0,0],'k:','LineWidth',1);
% plot3(factor2.*[0,0],factor2.*[0,0],factor2.*[-1,1],'k:','LineWidth',1);



%% plot intersection of all vectors on unit sphere in global space
% hold on
% factor=0.9;
% scatter3(factor*RightoptL(1),factor*RightoptL(2),factor*RightoptL(3),150,'k','LineWidth',2);
% scatter3(factor*RightASCr(1),factor*RightASCr(2),factor*RightASCr(3),150,'b','LineWidth',2);
% scatter3(factor*RightPSCr(1),factor*RightPSCr(2),factor*RightPSCr(3),150,'g','LineWidth',2);
% scatter3(factor*RightLSCr(1),factor*RightLSCr(2),factor*RightLSCr(3),150,'r','LineWidth',2);
% scatter3(factor*LeftoptL(1),factor*LeftoptL(2),factor*LeftoptL(3),150,'k');
% scatter3(factor*LeftASCr(1),factor*LeftASCr(2),factor*LeftASCr(3),150,'b');
% scatter3(factor*LeftPSCr(1),factor*LeftPSCr(2),factor*LeftPSCr(3),150,'g');
% scatter3(factor*LeftLSCr(1),factor*LeftLSCr(2),factor*LeftLSCr(3),150,'r');
% 
%% calculate angle between LB and MR-LR
% dif=(MRr-LRr);
% v1=[dif(1),dif(2),0];
% 
% angle=rad2deg(acos(dif'*v1'/(norm(v1)*norm(dif))));
% quiver3(LRr(1),LRr(2),LRr(3),1*dif(1),1*dif(2),1*dif(3),'c--','LineWidth',2);

%% calculate angle between LB and ASC normal
LB_RightASCr=acosd(dot(LBr,RightASCr)/(norm(LBr)*norm(RightASCr)))
LB_LeftASCr=acosd(dot(LBr,LeftASCr)/(norm(LBr)*norm(LeftASCr)))

LB_RightoptL=acosd(dot(LBr,RightoptL)/(norm(LBr)*norm(RightoptL)))

RightSMr_RightoptL=acosd(dot(RightSMr,RightoptL)/(norm(RightSMr)*norm(RightoptL)))

RightoptNr_RightoptL=acosd(dot(RightoptNr,RightoptL)/(norm(RightoptNr)*norm(RightoptL)))

RightASCr_RightPSCr=acosd(dot(RightASCr,RightPSCr)/(norm(RightASCr)*norm(RightPSCr)))
RightASCr_RightLSCr=acosd(dot(RightASCr,RightLSCr)/(norm(RightASCr)*norm(RightLSCr)))
RightPSCr_RightLSCr=acosd(dot(RightPSCr,RightLSCr)/(norm(RightPSCr)*norm(RightLSCr)))

LeftASCr_LeftPSCr=acosd(dot(LeftASCr,LeftPSCr)/(norm(LeftASCr)*norm(LeftPSCr)))
LeftASCr_LeftLSCr=acosd(dot(LeftASCr,LeftLSCr)/(norm(LeftASCr)*norm(LeftLSCr)))
LeftPSCr_LeftLSCr=acosd(dot(LeftPSCr,LeftLSCr)/(norm(LeftPSCr)*norm(LeftLSCr)))

RightASCr_LeftASCr=acosd(dot(RightASCr,LeftASCr)/(norm(RightASCr)*norm(LeftASCr)))

RightSMr_RightSM1r=acosd(dot(RightSMr,RightSM1r)/(norm(RightSMr)*norm(RightSM1r)))
RightSMr_RightSM2r=acosd(dot(RightSMr,RightSM2r)/(norm(RightSMr)*norm(RightSM2r)))
RightSMr_RightSM3r=acosd(dot(RightSMr,RightSM3r)/(norm(RightSMr)*norm(RightSM3r)))

RightSMr_RightoptL=acosd(dot(RightSMr,RightoptL)/(norm(RightSMr)*norm(RightoptL)))
LeftSMr_LeftoptL=acosd(dot(LeftSMr,LeftoptL)/(norm(LeftSMr)*norm(LeftoptL)))

RightUTr_RightoptL=acosd(dot(RightUTr,RightoptL)/(norm(RightUTr)*norm(RightoptL)))
LeftUTr_LeftoptL=acosd(dot(LeftUTr,LeftoptL)/(norm(LeftUTr)*norm(LeftoptL)))

RightSMr_RightSMhr=acosd(dot(RightSMr,RightSMhr)/(norm(RightSMr)*norm(RightSMhr)))

RightUTr_RightSMhr=acosd(dot(RightUTr,RightSMhr)/(norm(RightUTr)*norm(RightSMhr)))

RightUTr_RightSMr=acosd(dot(RightUTr,RightSMr)/(norm(RightUTr)*norm(RightSMr)))

LeftUTr_LeftSMr=acosd(dot(LeftUTr,LeftSMr)/(norm(LeftUTr)*norm(LeftSMr)))


LeftoptL_RightoptL=acosd(dot(LeftoptL,RightoptL)/(norm(LeftoptL)*norm(RightoptL)))

LeftASCr_RightPSCr=acosd(dot(LeftASCr,RightPSCr)/(norm(LeftASCr)*norm(RightPSCr)))
%% calcuate azimuth abd elevation of vector in global space
[aziLB,elevLB,r]=cart2sph(LBr(1),LBr(2),LBr(3));
[aziRightoptL,elevRightoptL,r]=cart2sph(RightoptL(1),RightoptL(2),RightoptL(3));
[aziRightASC,elevRightASC,r]=cart2sph(RightASCr(1),RightASCr(2),RightASCr(3));
[aziRightPSC,elevRightPSC,r]=cart2sph(RightPSCr(1),RightPSCr(2),RightPSCr(3));
[aziRightLSC,elevRightLSC,r]=cart2sph(RightLSCr(1),RightLSCr(2),RightLSCr(3));
[aziRightSM,elevRightSM,r]=cart2sph(RightSMr(1),RightSMr(2),RightSMr(3));
[aziRightSMh,elevRightSMh,r]=cart2sph(RightSMhr(1),RightSMhr(2),RightSMhr(3));
[aziRightUT,elevRightUT,r]=cart2sph(RightUTr(1),RightUTr(2),RightUTr(3));
[aziLeftoptL,elevLeftoptL,r]=cart2sph(LeftoptL(1),LeftoptL(2),LeftoptL(3));
[aziLeftASC,elevLeftASC,r]=cart2sph(LeftASCr(1),LeftASCr(2),LeftASCr(3));
[aziLeftPSC,elevLeftPSC,r]=cart2sph(LeftPSCr(1),LeftPSCr(2),LeftPSCr(3));
[aziLeftLSC,elevLeftLSC,r]=cart2sph(LeftLSCr(1),LeftLSCr(2),LeftLSCr(3));
[aziLeftSM,elevLeftSM,r]=cart2sph(LeftSMr(1),LeftSMr(2),LeftSMr(3));
[aziLeftSMh,elevLeftSMh,r]=cart2sph(LeftSMhr(1),LeftSMhr(2),LeftSMhr(3));
[aziLeftUT,elevLeftUT,r]=cart2sph(LeftUTr(1),LeftUTr(2),LeftUTr(3));


%% plot vectors
f3=figure;
set(0,'CurrentFigure',f3);
hold on
hLB=scatter(rad2deg(aziLB),rad2deg(elevLB),100,'m','filled'); 
hopt=scatter(rad2deg(aziRightoptL),rad2deg(elevRightoptL),100,'k','LineWidth',2); 
hASC=scatter(rad2deg(aziRightASC),rad2deg(elevRightASC),100,'b','LineWidth',2); 
hPSC=scatter(rad2deg(aziRightPSC),rad2deg(elevRightPSC),100,'g','LineWidth',2); 
hLSC=scatter(rad2deg(aziRightLSC),rad2deg(elevRightLSC),100,'r','LineWidth',2);
hSM=scatter(rad2deg(aziRightSM),rad2deg(elevRightSM),100,'y','LineWidth',2);
%hSMh=scatter(rad2deg(aziRightSMh),rad2deg(elevRightSMh),100,'c','LineWidth',2);
hUT=scatter(rad2deg(aziRightUT),rad2deg(elevRightUT),100,'MarkerEdgeColor',[0 .7 .7],'MarkerFaceColor',[1 1 1],'LineWidth',2);

scatter(rad2deg(aziLeftoptL),rad2deg(elevLeftoptL),100,'k','filled'); 
scatter(rad2deg(aziLeftASC),rad2deg(elevLeftASC),100,'b','filled'); 
scatter(rad2deg(aziLeftPSC),rad2deg(elevLeftPSC),100,'g','filled'); 
scatter(rad2deg(aziLeftLSC),rad2deg(elevLeftLSC),100,'r','filled'); 
scatter(rad2deg(aziLeftSM),rad2deg(elevLeftSM),100,'y','filled'); 
%scatter(rad2deg(aziLeftSMh),rad2deg(elevLeftSMh),100,'c','filled');
scatter(rad2deg(aziLeftUT),rad2deg(elevLeftUT),100,'MarkerEdgeColor',[0 .7 .7],'MarkerFaceColor',[0 .7 .7]);

plot([-180 180],[0 0],'m','LineWidth',1);
plot([0 0],[-90 90],'m','LineWidth',1);

xlim([-180 180]); ylim([-90 90]);
ax = gca;
axis equal
set(gcf, 'color', [1 1 1]);
legend([hLB,hASC,hPSC,hLSC,hSM,hUT,hopt],{'LB',' ASC',' PSC',' LSC',' SM',' UM',' Optic axis'},'Location','northoutside','Orientation','vertical','FontSize',14);
legend('boxoff');
ax.XTick = -180:45:180;
ax.YTick = -90:45:90;
ax.XTickLabel = {'behind','-135','left','-45','ahead','45','right','135','behind'};
ax.FontSize=12;
hold on
ax.YTickLabel = {'nadir','-45', 'horizon','45','zenith'};
xlabel('Azimuth (deg.)','FontSize',16);
%ax.YTickLabel = {'','','','',''};
hold on
ylabh=ylabel('Elevation (deg.)','FontSize',16);
%set(ylabh, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
axis equal
hold off
%% save workspace

save('summary anatomy');

%% plot translation and rotation axes along with visual field
set(0,'CurrentFigure',f2);
pause(1);
channels=sphereGlobalAxes();

%% calculate rotation axes based on optic axis and translation axes
% pT1=[-channels.yT1(end)/sqrt(channels.xT1(end).^2+channels.yT1(end).^2),channels.xT1(end)/sqrt(channels.xT1(end).^2+channels.yT1(end).^2),0];     % a point on the plane perpendicular to translation axis 1
% pT2=[-channels.yT2(end)/sqrt(channels.xT2(end).^2+channels.yT2(end).^2),channels.xT2(end)/sqrt(channels.xT2(end).^2+channels.yT2(end).^2),0];     % a point on the plane perpendicular to translation axis 2
% pT3=[-channels.yT3(end)/sqrt(channels.xT3(end).^2+channels.yT3(end).^2),channels.xT3(end)/sqrt(channels.xT3(end).^2+channels.yT3(end).^2),0];     % a point on the plane perpendicular to translation axis 3
% pT4=[-channels.yT4(end)/sqrt(channels.xT4(end).^2+channels.yT4(end).^2),channels.xT4(end)/sqrt(channels.xT4(end).^2+channels.yT4(end).^2),0];     % a point on the plane perpendicular to translation axis 4
% 
% pOpt=[-RightoptL(2)/sqrt(RightoptL(1).^2+RightoptL(2).^2),RightoptL(1)/sqrt(RightoptL(1).^2+RightoptL(2).^2),0];     % a point on the plane perpendicular to the optic axis
% 
% [P,N,check]=plane_intersect(RightoptL,pOpt,[channels.xT1(end),channels.yT1(end),channels.zT1(end)],pT1)





