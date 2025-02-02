function [] = VectorMovie(Normal,VizAngle)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Inputs
%   1. SCC normal
%   2. OpAxis
%   3. VisAngle NOTE: For this code VizAngle will be assumed to be in
%   radians.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RV=2;
%Computing alphabar, beta given the normal and the optical axis.
[alphabar,beta]=AB(Normal);

temp=linspace(0,-beta+3*pi/2,12);
CL=temp(2)-temp(1);


n=1000; %Number of discretization points for the curve.
t=linspace(0,2*pi,n); %parametrization coordinates for the curve
dt=t(2)-t(1);

%% Creation of Curves in Visual Space

%This part of the code creates the initial curves in the visual space about
%which rotations are taken.

Normal=[cos(alphabar)*sin(beta),sin(alphabar)*sin(beta),-cos(beta)]';

%The curve given by the rotation if (alphabar,beta)=(0,0)
C1=[RV*cos(t)*sin(CL);RV*sin(t)*sin(CL);-RV*cos(CL)*ones(1,n)];
C2=[RV*cos(t)*sin(2*CL);RV*sin(t)*sin(2*CL);-RV*cos(2*CL)*ones(1,n)];
C3=[RV*cos(t)*sin(3*CL);RV*sin(t)*sin(3*CL);-RV*cos(3*CL)*ones(1,n)];
C4=[RV*cos(t)*sin(4*CL);RV*sin(t)*sin(4*CL);-RV*cos(4*CL)*ones(1,n)];
C5=[RV*cos(t)*sin(5*CL);RV*sin(t)*sin(5*CL);-RV*cos(5*CL)*ones(1,n)];
C6=[RV*cos(t)*sin(6*CL);RV*sin(t)*sin(6*CL);-RV*cos(6*CL)*ones(1,n)];
C7=[RV*cos(t)*sin(7*CL);RV*sin(t)*sin(7*CL);-RV*cos(7*CL)*ones(1,n)];
C8=[RV*cos(t)*sin(8*CL);RV*sin(t)*sin(8*CL);-RV*cos(8*CL)*ones(1,n)];
C9=[RV*cos(t)*sin(9*CL);RV*sin(t)*sin(9*CL);-RV*cos(9*CL)*ones(1,n)];
C10=[RV*cos(t)*sin(10*CL);RV*sin(t)*sin(10*CL);-RV*cos(10*CL)*ones(1,n)];

%Rotation matrices for rotating the curve to lie around the optical
%singularity.

Rot1=[cos(beta), 0, -sin(beta); 0, 1, 0; sin(beta), 0, cos(beta)];
Rot2=[cos(alphabar), -sin(alphabar), 0; sin(alphabar), cos(alphabar), 0; 0, 0, 1];

%Curve about the optical singularity.

C1=Rot2*Rot1*C1;
C2=Rot2*Rot1*C2;
C3=Rot2*Rot1*C3;
C4=Rot2*Rot1*C4;
C5=Rot2*Rot1*C5;
C6=Rot2*Rot1*C6;
C7=Rot2*Rot1*C7;
C8=Rot2*Rot1*C8;
C9=Rot2*Rot1*C9;
C10=Rot2*Rot1*C10;

%Removal of data not in visual space.
sv1=RV*acos(-C1(3,:)/RV);
sv2=RV*acos(-C2(3,:)/RV);
sv3=RV*acos(-C3(3,:)/RV);
sv4=RV*acos(-C4(3,:)/RV);
sv5=RV*acos(-C5(3,:)/RV);
sv6=RV*acos(-C6(3,:)/RV);
sv7=RV*acos(-C7(3,:)/RV);
sv8=RV*acos(-C8(3,:)/RV);
sv9=RV*acos(-C9(3,:)/RV);
sv10=RV*acos(-C10(3,:)/RV);

InRange1=sv1>=RV*(pi-VizAngle/2);
InRange2=sv2>=RV*(pi-VizAngle/2);
InRange3=sv3>=RV*(pi-VizAngle/2);
InRange4=sv4>=RV*(pi-VizAngle/2);
InRange5=sv5>=RV*(pi-VizAngle/2);
InRange6=sv6>=RV*(pi-VizAngle/2);
InRange7=sv7>=RV*(pi-VizAngle/2);
InRange8=sv8>=RV*(pi-VizAngle/2);
InRange9=sv9>=RV*(pi-VizAngle/2);
InRange10=sv10>=RV*(pi-VizAngle/2);

C1(1,~InRange1)=NaN;
C2(1,~InRange2)=NaN;
C3(1,~InRange3)=NaN;
C4(1,~InRange4)=NaN;
C5(1,~InRange5)=NaN;
C6(1,~InRange6)=NaN;
C7(1,~InRange7)=NaN;
C8(1,~InRange8)=NaN;
C9(1,~InRange9)=NaN;
C10(1,~InRange10)=NaN;

%% Calculation of tangent vectors in visual space
Dt=DiffX(n)/dt;

C1X=(Dt*C1(1,:)')';
C1Y=(Dt*C1(2,:)')';
C1Z=(Dt*C1(3,:)')';
C2X=(Dt*C2(1,:)')';
C2Y=(Dt*C2(2,:)')';
C2Z=(Dt*C2(3,:)')';
C3X=(Dt*C3(1,:)')';
C3Y=(Dt*C3(2,:)')';
C3Z=(Dt*C3(3,:)')';
C4X=(Dt*C2(1,:)')';
C4Y=(Dt*C2(2,:)')';
C4Z=(Dt*C2(3,:)')';
C5X=(Dt*C5(1,:)')';
C5Y=(Dt*C5(2,:)')';
C5Z=(Dt*C5(3,:)')';
C6X=(Dt*C6(1,:)')';
C6Y=(Dt*C6(2,:)')';
C6Z=(Dt*C6(3,:)')';
C7X=(Dt*C7(1,:)')';
C7Y=(Dt*C7(2,:)')';
C7Z=(Dt*C7(3,:)')';
C8X=(Dt*C8(1,:)')';
C8Y=(Dt*C8(2,:)')';
C8Z=(Dt*C8(3,:)')';
C9X=(Dt*C9(1,:)')';
C9Y=(Dt*C9(2,:)')';
C9Z=(Dt*C9(3,:)')';
C10X=(Dt*C10(1,:)')';
C10Y=(Dt*C10(2,:)')';
C10Z=(Dt*C10(3,:)')';

scale=1/5;

normC1=(C1X.^2+C1Y.^2+C1Z.^2).^(1/2)/scale;
normC2=(C2X.^2+C2Y.^2+C2Z.^2).^(1/2)/scale;
normC3=(C3X.^2+C3Y.^2+C3Z.^2).^(1/2)/scale;
normC4=(C4X.^2+C4Y.^2+C4Z.^2).^(1/2)/scale;
normC5=(C5X.^2+C5Y.^2+C5Z.^2).^(1/2)/scale;
normC6=(C6X.^2+C6Y.^2+C6Z.^2).^(1/2)/scale;
normC7=(C7X.^2+C7Y.^2+C7Z.^2).^(1/2)/scale;
normC8=(C8X.^2+C8Y.^2+C8Z.^2).^(1/2)/scale;
normC9=(C9X.^2+C9Y.^2+C9Z.^2).^(1/2)/scale;
normC10=(C10X.^2+C10Y.^2+C10Z.^2).^(1/2)/scale;

C1X=C1X./normC1;
C1Y=C1Y./normC1;
C1Z=C1Z./normC1;
C2X=C2X./normC2;
C2Y=C2Y./normC2;
C2Z=C2Z./normC2;
C3X=C3X./normC3;
C3Y=C3Y./normC3;
C3Z=C3Z./normC3;
C4X=C4X./normC4;
C4Y=C4Y./normC4;
C4Z=C4Z./normC4;
C5X=C5X./normC5;
C5Y=C5Y./normC5;
C5Z=C5Z./normC5;
C6X=C6X./normC6;
C6Y=C6Y./normC6;
C6Z=C6Z./normC6;
C7X=C7X./normC7;
C7Y=C7Y./normC7;
C7Z=C7Z./normC7;
C8X=C8X./normC8;
C8Y=C8Y./normC8;
C8Z=C8Z./normC8;
C9X=C9X./normC9;
C9Y=C9Y./normC9;
C9Z=C9Z./normC9;
C10X=C10X./normC10;
C10Y=C10Y./normC10;
C10Z=C10Z./normC10;


%% Plot of curves in visual space with normal
figure('Position',get(0,'ScreenSize'));
hold on;

s=linspace((pi-VizAngle/2),pi,100);
th=linspace(0,2*pi,100);

[S,TH]=meshgrid(s,th);
surf(RV*cos(TH).*sin(S),RV*sin(TH).*sin(S),-RV*cos(S),'EdgeColor','none');
alpha(.3);

%Plotting Optical singularity
plot3(RV*cos(alphabar)*sin(beta),RV*sin(alphabar)*sin(beta),-RV*cos(beta),'b.','MarkerSize',40);

%Plotting curve in visual space.
plot3(C1(1,:),C1(2,:),C1(3,:),'k','LineWidth',2);
plot3(C2(1,:),C2(2,:),C2(3,:),'k','LineWidth',2);
plot3(C3(1,:),C3(2,:),C3(3,:),'k','LineWidth',2);
plot3(C4(1,:),C4(2,:),C4(3,:),'k','LineWidth',2);
plot3(C5(1,:),C5(2,:),C5(3,:),'k','LineWidth',2);
plot3(C6(1,:),C6(2,:),C6(3,:),'k','LineWidth',2);
plot3(C7(1,:),C7(2,:),C7(3,:),'k','LineWidth',2);
plot3(C8(1,:),C8(2,:),C8(3,:),'k','LineWidth',2);
plot3(C9(1,:),C9(2,:),C9(3,:),'k','LineWidth',2);
plot3(C10(1,:),C10(2,:),C10(3,:),'k','LineWidth',2);

% %Plotting tangents to curves
% scale=0;
% step=40;
% quiver3(C1(1,1:step:n),C1(2,1:step:n),C1(3,1:step:n),C1X(1:step:n),C1Y(1:step:n),C1Z(1:step:n),scale,'b','LineWidth',2);
% quiver3(C2(1,1:step:n),C2(2,1:step:n),C2(3,1:step:n),C2X(1:step:n),C2Y(1:step:n),C2Z(1:step:n),scale,'b','LineWidth',2);
% quiver3(C3(1,1:step:n),C3(2,1:step:n),C3(3,1:step:n),C3X(1:step:n),C3Y(1:step:n),C3Z(1:step:n),scale,'b','LineWidth',2);
% quiver3(C4(1,1:step:n),C4(2,1:step:n),C4(3,1:step:n),C4X(1:step:n),C4Y(1:step:n),C4Z(1:step:n),scale,'b','LineWidth',2);
% quiver3(C5(1,1:step:n),C5(2,1:step:n),C5(3,1:step:n),C5X(1:step:n),C5Y(1:step:n),C5Z(1:step:n),scale,'b','LineWidth',2);
% quiver3(C6(1,1:step:n),C6(2,1:step:n),C6(3,1:step:n),C6X(1:step:n),C6Y(1:step:n),C6Z(1:step:n),scale,'b','LineWidth',2);
% quiver3(C7(1,1:step:n),C7(2,1:step:n),C7(3,1:step:n),C7X(1:step:n),C7Y(1:step:n),C7Z(1:step:n),scale,'b','LineWidth',2);
% quiver3(C8(1,1:step:n),C8(2,1:step:n),C8(3,1:step:n),C8X(1:step:n),C8Y(1:step:n),C8Z(1:step:n),scale,'b','LineWidth',2);
% quiver3(C9(1,1:step:n),C9(2,1:step:n),C9(3,1:step:n),C9X(1:step:n),C9Y(1:step:n),C9Z(1:step:n),scale,'b','LineWidth',2);
% quiver3(C10(1,1:step:n),C10(2,1:step:n),C10(3,1:step:n),C10X(1:step:n),C10Y(1:step:n),C10Z(1:step:n),scale,'b','LineWidth',2);

%Plotting Representation of mouse eye
s=linspace(0,pi,100);
th=linspace(0,2*pi,100);
[S,TH]=meshgrid(s,th);
surf(RV/10*cos(TH).*sin(S),RV/10*sin(TH).*sin(S),-RV/10*cos(S)-RV/10,'FaceColor','green','EdgeColor','none');
%axis([-RV,RV,-RV,RV,-2*RV/10,RV]);
alpha(.3);

%Plotting Representation of mouse lens
s=linspace(0,pi,100);
th=linspace(0,2*pi,100);
[S,TH]=meshgrid(s,th);
surf(RV/28*cos(TH).*sin(S),RV/28*sin(TH).*sin(S),-RV/28*cos(S)-RV/28,'FaceColor','blue','EdgeColor','none');
alpha(.3);

%Plotting Representation of retina
s=linspace(0,deg2rad(105),100);
th=linspace(0,2*pi,100);
[S,TH]=meshgrid(s,th);
surf((RV/12)*cos(TH).*sin(S),(RV/12)*sin(TH).*sin(S),-(RV/12)*cos(S)-RV/10,'FaceColor','red','EdgeColor','none');
alpha(.5);

%Ploting normal axis
t=linspace(0,3/2*RV,100);
plot3(t*Normal(1),t*Normal(2),t*Normal(3),'k','LineWidth',4);

%Plotting line connecting optical singularity on retina to lens
t=linspace(0,1,100)';
betashift=-2*deg2rad(105)/pi*beta+2*deg2rad(105);
alphashift=alphabar+pi;

L(:,1)=t*RV/12*cos(alphashift)*sin(betashift);
L(:,2)=t*RV/12*sin(alphashift)*sin(betashift);
L(:,3)=t*(-RV/12*cos(betashift)-RV/10);
plot3(L(:,1),L(:,2),L(:,3),'k','LineWidth',4);

%Plotting Vector Field on retina
%quiver3(C1(1,:),C1(2,:),C1(3,:),vecXC1,vecYC1,vecZC1);

axis equal;
view(90,18);
axis off;
set(gcf, 'color', [1 1 1])
hold off;

% %% Movie in Visual Space
%
% rate=10;
% maxvalue=n/rate; %frame rate
%
% AnglePar=linspace(0,2*pi,n);
%
% for i=1:maxvalue,
%     close all
%     j=i*rate;
%
% figure('Position',get(0,'ScreenSize'));
% hold on;
%
% s=linspace((pi-VizAngle/2),pi,100);
% th=linspace(0,2*pi,100);
%
% [S,TH]=meshgrid(s,th);
% surf(RV*cos(TH).*sin(S),RV*sin(TH).*sin(S),-RV*cos(S),'EdgeColor','none');
% alpha(.3);
%
% %Plotting Optical singularity
% plot3(RV*cos(alphabar)*sin(beta),RV*sin(alphabar)*sin(beta),-RV*cos(beta),'b.','MarkerSize',40);
%
% %Plotting curve in visual space.
% plot3(C1(1,:),C1(2,:),C1(3,:),'k','LineWidth',2);
% plot3(C2(1,:),C2(2,:),C2(3,:),'k','LineWidth',2);
% plot3(C3(1,:),C3(2,:),C3(3,:),'k','LineWidth',2);
% plot3(C4(1,:),C4(2,:),C4(3,:),'k','LineWidth',2);
% plot3(C5(1,:),C5(2,:),C5(3,:),'k','LineWidth',2);
% plot3(C6(1,:),C6(2,:),C6(3,:),'k','LineWidth',2);
% plot3(C7(1,:),C7(2,:),C7(3,:),'k','LineWidth',2);
% plot3(C8(1,:),C8(2,:),C8(3,:),'k','LineWidth',2);
% plot3(C9(1,:),C9(2,:),C9(3,:),'k','LineWidth',2);
% plot3(C10(1,:),C10(2,:),C10(3,:),'k','LineWidth',2);
%
% % %Plotting tangents to curves
% % scale=0;
% % step=40;
% % quiver3(C1(1,1:step:n),C1(2,1:step:n),C1(3,1:step:n),C1X(1:step:n),C1Y(1:step:n),C1Z(1:step:n),scale,'b','LineWidth',2);
% % quiver3(C2(1,1:step:n),C2(2,1:step:n),C2(3,1:step:n),C2X(1:step:n),C2Y(1:step:n),C2Z(1:step:n),scale,'b','LineWidth',2);
% % quiver3(C3(1,1:step:n),C3(2,1:step:n),C3(3,1:step:n),C3X(1:step:n),C3Y(1:step:n),C3Z(1:step:n),scale,'b','LineWidth',2);
% % quiver3(C4(1,1:step:n),C4(2,1:step:n),C4(3,1:step:n),C4X(1:step:n),C4Y(1:step:n),C4Z(1:step:n),scale,'b','LineWidth',2);
% % quiver3(C5(1,1:step:n),C5(2,1:step:n),C5(3,1:step:n),C5X(1:step:n),C5Y(1:step:n),C5Z(1:step:n),scale,'b','LineWidth',2);
% % quiver3(C6(1,1:step:n),C6(2,1:step:n),C6(3,1:step:n),C6X(1:step:n),C6Y(1:step:n),C6Z(1:step:n),scale,'b','LineWidth',2);
% % quiver3(C7(1,1:step:n),C7(2,1:step:n),C7(3,1:step:n),C7X(1:step:n),C7Y(1:step:n),C7Z(1:step:n),scale,'b','LineWidth',2);
% % quiver3(C8(1,1:step:n),C8(2,1:step:n),C8(3,1:step:n),C8X(1:step:n),C8Y(1:step:n),C8Z(1:step:n),scale,'b','LineWidth',2);
% % quiver3(C9(1,1:step:n),C9(2,1:step:n),C9(3,1:step:n),C9X(1:step:n),C9Y(1:step:n),C9Z(1:step:n),scale,'b','LineWidth',2);
% % quiver3(C10(1,1:step:n),C10(2,1:step:n),C10(3,1:step:n),C10X(1:step:n),C10Y(1:step:n),C10Z(1:step:n),scale,'b','LineWidth',2);
%
% %Plotting Representation of mouse eye
% s=linspace(0,pi,100);
% th=linspace(0,2*pi,100);
% [S,TH]=meshgrid(s,th);
% surf(RV/10*cos(TH).*sin(S),RV/10*sin(TH).*sin(S),-RV/10*cos(S)-RV/10,'FaceColor','green','EdgeColor','none');
% %axis([-RV,RV,-RV,RV,-2*RV/10,RV]);
% alpha(.3);
%
% %Plotting Representation of mouse lens
% s=linspace(0,pi,100);
% th=linspace(0,2*pi,100);
% [S,TH]=meshgrid(s,th);
% surf(RV/28*cos(TH).*sin(S),RV/28*sin(TH).*sin(S),-RV/28*cos(S)-RV/28,'FaceColor','blue','EdgeColor','none');
% alpha(.3);
%
% %Plotting Representation of retina
% s=linspace(0,deg2rad(105),100);
% th=linspace(0,2*pi,100);
% [S,TH]=meshgrid(s,th);
% surf((RV/12)*cos(TH).*sin(S),(RV/12)*sin(TH).*sin(S),-(RV/12)*cos(S)-RV/10,'FaceColor','red','EdgeColor','none');
% alpha(.5);
%
% %Ploting normal axis
% t=linspace(0,3/2*RV,100);
% plot3(t*Normal(1),t*Normal(2),t*Normal(3),'k','LineWidth',4);
%
% %Plotting line connecting optical singularity on retina to lens
% t=linspace(0,1,100)';
% betashift=-2*deg2rad(105)/pi*beta+2*deg2rad(105);
% alphashift=alphabar+pi;
%
% L(:,1)=t*RV/12*cos(alphashift)*sin(betashift);
% L(:,2)=t*RV/12*sin(alphashift)*sin(betashift);
% L(:,3)=t*(-RV/12*cos(betashift)-RV/10);
% plot3(L(:,1),L(:,2),L(:,3),'k','LineWidth',4);
%
% %Plotting Vector Field on retina
% %quiver3(C1(1,:),C1(2,:),C1(3,:),vecXC1,vecYC1,vecZC1);
%
% %Plotting Point on curve in visual space.
% plot3(C1(1,j),C1(2,j),C1(3,j),'b.','MarkerSize',40);
% plot3(C2(1,j),C2(2,j),C2(3,j),'b.','MarkerSize',40);
% plot3(C3(1,j),C3(2,j),C3(3,j),'b.','MarkerSize',40);
% plot3(C4(1,j),C4(2,j),C4(3,j),'b.','MarkerSize',40);
% plot3(C5(1,j),C5(2,j),C5(3,j),'b.','MarkerSize',40);
% plot3(C6(1,j),C6(2,j),C6(3,j),'b.','MarkerSize',40);
% plot3(C7(1,j),C7(2,j),C7(3,j),'b.','MarkerSize',40);
% plot3(C8(1,j),C8(2,j),C8(3,j),'b.','MarkerSize',40);
% plot3(C9(1,j),C9(2,j),C9(3,j),'b.','MarkerSize',40);
% plot3(C10(1,j),C10(2,j),C10(3,j),'b.','MarkerSize',40);
%
% axis equal;
% view(90,18);
% axis off;
% set(gcf, 'color', [1 1 1])
% hold off;
%
% F(i)=getframe(gca);
%     F(i)=getframe(gca);
%     hold off;
%     movie2avi(F,'RotationVisualSpace.avi');
%
% end

%% Reading in data from retina
%This part of the code computes all of the data needed for the flattening
%map.

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

ntemp=40;
%% Creation of curves on retina

% S, TH coordinates for curves on the retina
sret1=-2*M/VizAngle*acos(-C1(3,:)/RV)+2*M*pi/VizAngle;
thret1=atan2(C1(2,:),C1(1,:))+pi;
sret2=-2*M/VizAngle*acos(-C2(3,:)/RV)+2*M*pi/VizAngle;
thret2=atan2(C2(2,:),C2(1,:))+pi;
sret3=-2*M/VizAngle*acos(-C3(3,:)/RV)+2*M*pi/VizAngle;
thret3=atan2(C3(2,:),C3(1,:))+pi;
sret4=-2*M/VizAngle*acos(-C4(3,:)/RV)+2*M*pi/VizAngle;
thret4=atan2(C4(2,:),C4(1,:))+pi;
sret5=-2*M/VizAngle*acos(-C5(3,:)/RV)+2*M*pi/VizAngle;
thret5=atan2(C5(2,:),C5(1,:))+pi;
sret6=-2*M/VizAngle*acos(-C6(3,:)/RV)+2*M*pi/VizAngle;
thret6=atan2(C6(2,:),C6(1,:))+pi;
sret7=-2*M/VizAngle*acos(-C7(3,:)/RV)+2*M*pi/VizAngle;
thret7=atan2(C7(2,:),C7(1,:))+pi;
sret8=-2*M/VizAngle*acos(-C8(3,:)/RV)+2*M*pi/VizAngle;
thret8=atan2(C8(2,:),C8(1,:))+pi;
sret9=-2*M/VizAngle*acos(-C9(3,:)/RV)+2*M*pi/VizAngle;
thret9=atan2(C9(2,:),C9(1,:))+pi;
sret10=-2*M/VizAngle*acos(-C10(3,:)/RV)+2*M*pi/VizAngle;
thret10=atan2(C10(2,:),C10(1,:))+pi;

% 3-dimensional coordinate of curves on the retina.
C1ret(1,:)=R*cos(thret1).*sin(sret1/R);
C1ret(2,:)=R*sin(thret1).*sin(sret1/R);
C1ret(3,:)=-R*cos(sret1/R);
C2ret(1,:)=R*cos(thret2).*sin(sret2/R);
C2ret(2,:)=R*sin(thret2).*sin(sret2/R);
C2ret(3,:)=-R*cos(sret2/R);
C3ret(1,:)=R*cos(thret3).*sin(sret3/R);
C3ret(2,:)=R*sin(thret3).*sin(sret3/R);
C3ret(3,:)=-R*cos(sret3/R);
C4ret(1,:)=R*cos(thret4).*sin(sret4/R);
C4ret(2,:)=R*sin(thret4).*sin(sret4/R);
C4ret(3,:)=-R*cos(sret4/R);
C5ret(1,:)=R*cos(thret5).*sin(sret5/R);
C5ret(2,:)=R*sin(thret5).*sin(sret5/R);
C5ret(3,:)=-R*cos(sret5/R);
C6ret(1,:)=R*cos(thret6).*sin(sret6/R);
C6ret(2,:)=R*sin(thret6).*sin(sret6/R);
C6ret(3,:)=-R*cos(sret6/R);
C7ret(1,:)=R*cos(thret7).*sin(sret7/R);
C7ret(2,:)=R*sin(thret7).*sin(sret7/R);
C7ret(3,:)=-R*cos(sret7/R);
C8ret(1,:)=R*cos(thret8).*sin(sret8/R);
C8ret(2,:)=R*sin(thret8).*sin(sret8/R);
C8ret(3,:)=-R*cos(sret8/R);
C9ret(1,:)=R*cos(thret9).*sin(sret9/R);
C9ret(2,:)=R*sin(thret9).*sin(sret9/R);
C9ret(3,:)=-R*cos(sret9/R);
C10ret(1,:)=R*cos(thret10).*sin(sret10/R);
C10ret(2,:)=R*sin(thret10).*sin(sret10/R);
C10ret(3,:)=-R*cos(sret10/R);

%% Creation of tangent vectors on the retina
C1retX=(Dt*C1ret(1,:)')';
C1retY=(Dt*C1ret(2,:)')';
C1retZ=(Dt*C1ret(3,:)')';
C2retX=(Dt*C2ret(1,:)')';
C2retY=(Dt*C2ret(2,:)')';
C2retZ=(Dt*C2ret(3,:)')';
C3retX=(Dt*C3ret(1,:)')';
C3retY=(Dt*C3ret(2,:)')';
C3retZ=(Dt*C3ret(3,:)')';
C4retX=(Dt*C2ret(1,:)')';
C4retY=(Dt*C2ret(2,:)')';
C4retZ=(Dt*C2ret(3,:)')';
C5retX=(Dt*C5ret(1,:)')';
C5retY=(Dt*C5ret(2,:)')';
C5retZ=(Dt*C5ret(3,:)')';
C6retX=(Dt*C6ret(1,:)')';
C6retY=(Dt*C6ret(2,:)')';
C6retZ=(Dt*C6ret(3,:)')';
C7retX=(Dt*C7ret(1,:)')';
C7retY=(Dt*C7ret(2,:)')';
C7retZ=(Dt*C7ret(3,:)')';
C8retX=(Dt*C8ret(1,:)')';
C8retY=(Dt*C8ret(2,:)')';
C8retZ=(Dt*C8ret(3,:)')';
C9retX=(Dt*C9ret(1,:)')';
C9retY=(Dt*C9ret(2,:)')';
C9retZ=(Dt*C9ret(3,:)')';
C10retX=(Dt*C10ret(1,:)')';
C10retY=(Dt*C10ret(2,:)')';
C10retZ=(Dt*C10ret(3,:)')';

scale=1/15;

normC1ret=(C1retX.^2+C1retY.^2+C1retZ.^2).^(1/2)/scale;
normC2ret=(C2retX.^2+C2retY.^2+C2retZ.^2).^(1/2)/scale;
normC3ret=(C3retX.^2+C3retY.^2+C3retZ.^2).^(1/2)/scale;
normC4ret=(C4retX.^2+C4retY.^2+C4retZ.^2).^(1/2)/scale;
normC5ret=(C5retX.^2+C5retY.^2+C5retZ.^2).^(1/2)/scale;
normC6ret=(C6retX.^2+C6retY.^2+C6retZ.^2).^(1/2)/scale;
normC7ret=(C7retX.^2+C7retY.^2+C7retZ.^2).^(1/2)/scale;
normC8ret=(C8retX.^2+C8retY.^2+C8retZ.^2).^(1/2)/scale;
normC9ret=(C9retX.^2+C9retY.^2+C9retZ.^2).^(1/2)/scale;
normC10ret=(C10retX.^2+C10retY.^2+C10retZ.^2).^(1/2)/scale;

C1retX=C1retX./normC1ret;
C1retY=C1retY./normC1ret;
C1retZ=C1retZ./normC1ret;
C2retX=C2retX./normC2ret;
C2retY=C2retY./normC2ret;
C2retZ=C2retZ./normC2ret;
C3retX=C3retX./normC3ret;
C3retY=C3retY./normC3ret;
C3retZ=C3retZ./normC3ret;
C4retX=C4retX./normC4ret;
C4retY=C4retY./normC4ret;
C4retZ=C4retZ./normC4ret;
C5retX=C5retX./normC5ret;
C5retY=C5retY./normC5ret;
C5retZ=C5retZ./normC5ret;
C6retX=C6retX./normC6ret;
C6retY=C6retY./normC6ret;
C6retZ=C6retZ./normC6ret;
C7retX=C7retX./normC7ret;
C7retY=C7retY./normC7ret;
C7retZ=C7retZ./normC7ret;
C8retX=C8retX./normC8ret;
C8retY=C8retY./normC8ret;
C8retZ=C8retZ./normC8ret;
C9retX=C9retX./normC9ret;
C9retY=C9retY./normC9ret;
C9retZ=C9retZ./normC9ret;
C10retX=C10retX./normC10ret;
C10retY=C10retY./normC10ret;
C10retZ=C10retZ./normC10ret;

%% Plotting curves on retina

figure('Position',get(0,'ScreenSize'));
hold on;

%Plotting Surface of retina
s=linspace(0,M,100);
th=linspace(0,2*pi,100);
[S,TH]=meshgrid(s,th);
surf(R*cos(TH).*sin(S/R),R*sin(TH).*sin(S/R),-R*cos(S/R),'FaceColor','red','EdgeColor','none');
alpha(.5);
camlight left;
lighting phong;

%Plotting Optical singularity on retina
betashift=-2*M/(R*pi)*beta+2*M/R;
alphashift=alphabar+pi;
plot3(R*cos(alphashift)*sin(betashift),R*sin(alphashift)*sin(betashift),-R*cos(betashift),'b.','MarkerSize',40);


%Plotting curves on retina.

plot3(C1ret(1,:),C1ret(2,:),C1ret(3,:),'k','LineWidth',2);
plot3(C2ret(1,:),C2ret(2,:),C2ret(3,:),'k','LineWidth',2);
plot3(C3ret(1,:),C3ret(2,:),C3ret(3,:),'k','LineWidth',2);
plot3(C4ret(1,:),C4ret(2,:),C4ret(3,:),'k','LineWidth',2);
plot3(C5ret(1,:),C5ret(2,:),C5ret(3,:),'k','LineWidth',2);
plot3(C6ret(1,:),C6ret(2,:),C6ret(3,:),'k','LineWidth',2);
plot3(C7ret(1,:),C7ret(2,:),C7ret(3,:),'k','LineWidth',2);
plot3(C8ret(1,:),C8ret(2,:),C8ret(3,:),'k','LineWidth',2);
plot3(C9ret(1,:),C9ret(2,:),C9ret(3,:),'k','LineWidth',2);
plot3(C10ret(1,:),C10ret(2,:),C10ret(3,:),'k','LineWidth',2);


%Plotting tangents to curves
scale=0;
step=40;
quiver3(C1ret(1,1:step:n),C1ret(2,1:step:n),C1ret(3,1:step:n),C1retX(1:step:n),C1retY(1:step:n),C1retZ(1:step:n),scale,'b','LineWidth',2);
quiver3(C2ret(1,1:step:n),C2ret(2,1:step:n),C2ret(3,1:step:n),C2retX(1:step:n),C2retY(1:step:n),C2retZ(1:step:n),scale,'b','LineWidth',2);
quiver3(C3ret(1,1:step:n),C3ret(2,1:step:n),C3ret(3,1:step:n),C3retX(1:step:n),C3retY(1:step:n),C3retZ(1:step:n),scale,'b','LineWidth',2);
quiver3(C4ret(1,1:step:n),C4ret(2,1:step:n),C4ret(3,1:step:n),C4retX(1:step:n),C4retY(1:step:n),C4retZ(1:step:n),scale,'b','LineWidth',2);
quiver3(C5ret(1,1:step:n),C5ret(2,1:step:n),C5ret(3,1:step:n),C5retX(1:step:n),C5retY(1:step:n),C5retZ(1:step:n),scale,'b','LineWidth',2);
quiver3(C6ret(1,1:step:n),C6ret(2,1:step:n),C6ret(3,1:step:n),C6retX(1:step:n),C6retY(1:step:n),C6retZ(1:step:n),scale,'b','LineWidth',2);
quiver3(C7ret(1,1:step:n),C7ret(2,1:step:n),C7ret(3,1:step:n),C7retX(1:step:n),C7retY(1:step:n),C7retZ(1:step:n),scale,'b','LineWidth',2);
quiver3(C8ret(1,1:step:n),C8ret(2,1:step:n),C8ret(3,1:step:n),C8retX(1:step:n),C8retY(1:step:n),C8retZ(1:step:n),scale,'b','LineWidth',2);
quiver3(C9ret(1,1:step:n),C9ret(2,1:step:n),C9ret(3,1:step:n),C9retX(1:step:n),C9retY(1:step:n),C9retZ(1:step:n),scale,'b','LineWidth',2);
quiver3(C10ret(1,1:step:n),C10ret(2,1:step:n),C10ret(3,1:step:n),C10retX(1:step:n),C10retY(1:step:n),C10retZ(1:step:n),scale,'b','LineWidth',2);

axis off;
axis equal;
view(90,18);
set(gcf, 'color', [1 1 1]);

hold off

% %% Movie On retina
% rate=10;
% maxvalue=n/rate; %frame rate
%
% AnglePar=linspace(0,2*pi,n);
%
% for i=1:maxvalue,
%     j=i*rate;
%
%     close all;
%     figure('Position',get(0,'ScreenSize'));
%
%
%     hold on;
%
%     %Plotting Surface of retina
% s=linspace(0,M,100);
% th=linspace(0,2*pi,100);
% [S,TH]=meshgrid(s,th);
% surf(R*cos(TH).*sin(S/R),R*sin(TH).*sin(S/R),-R*cos(S/R),'FaceColor','red','EdgeColor','none');
% alpha(.5);
% camlight left;
% lighting phong;
%
% %Plotting Optical singularity on retina
% betashift=-2*M/(R*pi)*beta+2*M/R;
% alphashift=alphabar+pi;
% plot3(R*cos(alphashift)*sin(betashift),R*sin(alphashift)*sin(betashift),-R*cos(betashift),'b.','MarkerSize',40);
%
%
% %Plotting curves on retina.
%
% plot3(C1ret(1,:),C1ret(2,:),C1ret(3,:),'k','LineWidth',2);
% plot3(C2ret(1,:),C2ret(2,:),C2ret(3,:),'k','LineWidth',2);
% plot3(C3ret(1,:),C3ret(2,:),C3ret(3,:),'k','LineWidth',2);
% plot3(C4ret(1,:),C4ret(2,:),C4ret(3,:),'k','LineWidth',2);
% plot3(C5ret(1,:),C5ret(2,:),C5ret(3,:),'k','LineWidth',2);
% plot3(C6ret(1,:),C6ret(2,:),C6ret(3,:),'k','LineWidth',2);
% plot3(C7ret(1,:),C7ret(2,:),C7ret(3,:),'k','LineWidth',2);
% plot3(C8ret(1,:),C8ret(2,:),C8ret(3,:),'k','LineWidth',2);
% plot3(C9ret(1,:),C9ret(2,:),C9ret(3,:),'k','LineWidth',2);
% plot3(C10ret(1,:),C10ret(2,:),C10ret(3,:),'k','LineWidth',2);
%
% %Plotting Point on curve in visual space.
% plot3(C1ret(1,j),C1ret(2,j),C1ret(3,j),'b.','MarkerSize',40);
% plot3(C2ret(1,j),C2ret(2,j),C2ret(3,j),'b.','MarkerSize',40);
% plot3(C3ret(1,j),C3ret(2,j),C3ret(3,j),'b.','MarkerSize',40);
% plot3(C4ret(1,j),C4ret(2,j),C4ret(3,j),'b.','MarkerSize',40);
% plot3(C5ret(1,j),C5ret(2,j),C5ret(3,j),'b.','MarkerSize',40);
% plot3(C6ret(1,j),C6ret(2,j),C6ret(3,j),'b.','MarkerSize',40);
% plot3(C7ret(1,j),C7ret(2,j),C7ret(3,j),'b.','MarkerSize',40);
% plot3(C8ret(1,j),C8ret(2,j),C8ret(3,j),'b.','MarkerSize',40);
% plot3(C9ret(1,j),C9ret(2,j),C9ret(3,j),'b.','MarkerSize',40);
% plot3(C10ret(1,j),C10ret(2,j),C10ret(3,j),'b.','MarkerSize',40);
%
% axis square;
%     axis equal;
%     view(90,18);
%     axis off;
%     set(gcf, 'color', [1 1 1])
%
%     F(i)=getframe(gca);
%     F(i)=getframe(gca);
%     hold off;
%     movie2avi(F,'RotationRetina.avi');
%
% %Use the below compression on a window machine.
% %movie2avi(F,'RotationSpace.avi','compression','indeo5');
% end


%% Sector 1 data
% In this code we use the data from the flattened retina to create induce
% the vector field on the flattened retina.

s=linspace(0,M,ntemp+1);
s=linspace(s(2),M,ntemp);
th=linspace(a1,a2,ntemp);
[S,TH]=meshgrid(s,th);

%Construction of domain in visual optic space. This corresponds to the
%the coordinates on the hemisphere about the optical axis which the eye can
%see.

SV=-VizAngle*RV/(2*M)*S+pi*RV;
THV=TH+pi;

%Calculation of the vector field of rotations in the coordinates svis and
%thvis. The variables (t,phi) denote geodesic polar coordinate around the axis of
%rotation.

t=zeros(n,n);
phi=zeros(n,n);
vecS=zeros(n,n);
vecTH=zeros(n,n);


%This loop constructs the components of the vector field of rotations
%in the visual space of the eye.  The function InverseSTH determines the
%coordinates of (t,phi) given coordiantes svis and thvis and is necessary to
%calculate the components of the vector field.

for i=1:ntemp
    for j=1:ntemp,
        [t(i,j),phi(i,j)]=InverseSTH(alphabar,beta,RV,SV(i,j),THV(i,j));
        
        vecSV(i,j)=-cos(SV(i,j)/RV)*cos(beta)*sin(phi(i,j))*cos(THV(i,j)-alphabar)+...
            cos(SV(i,j)/RV)*cos(phi(i,j))*(sin(THV(i,j)-alphabar))-...
            sin(SV(i,j)/RV)*sin(phi(i,j))*sin(beta);
        
        vecSV(i,j)=RV*sin(t(i,j)/RV)*vecSV(i,j);
        
        vecTHV(i,j)=(cos(beta)*sin(phi(i,j))*sin(THV(i,j)-alphabar)+...
            cos(phi(i,j))*cos(alphabar-THV(i,j)));
        
        vecTHV(i,j)=RV*sin(t(i,j)/RV)*vecTHV(i,j);
        vecTHV(i,j)=1/(RV*sin(SV(i,j)/RV))*vecTHV(i,j);
        
    end
end

%Conversion of the vector fields in visual space to vector fields on the
%retina.

vecS=-2*M/(VizAngle*RV)*vecSV;
vecTH=vecTHV;

%We now push forward this vector onto the flattened retina.

StSec1Data=dlmread('StandardSector1_40');

StRHO=StSec1Data(1:ntemp,1:ntemp);
StF=StSec1Data(1:ntemp,ntemp+1:2*ntemp);
StS=StSec1Data(1:ntemp,2*ntemp+1:3*ntemp);
StTH=StSec1Data(1:ntemp,3*ntemp+1:4*ntemp);

StU=StRHO.*cos(StF);
StV=StRHO.*sin(StF);

%Calcution of derivatives
ds=StS(1,2)-StS(1,1);
DS=DiffX(ntemp)/ds;
StRHOS=(DS*StRHO')';
StFS=(DS*StF')';

dth=StTH(2,1)-StTH(1,1);
DTH=DiffX(ntemp)/dth;
StRHOTH=(DTH*StRHO);
StFTH=(DTH*StF);

StVecU=vecS.*(StRHOS.*cos(StF)-StRHO.*StFS.*sin(StF))+vecTH.*(StRHOTH.*cos(StF)-StRHO.*StFTH.*sin(StF));
StVecV=vecS.*(StRHOS.*sin(StF)+StRHO.*StFS.*cos(StF))+vecTH.*(StRHOTH.*sin(StF)+StRHO.*StFTH.*cos(StF));

temp=(StVecU.^2+StVecV.^2).^(1/2);

StVecU=StVecU./temp;
StVecV=StVecV./temp;

S1=S;
TH1=TH;
StU1=StU;
StV1=StV;
StVecU1=StVecU;
StVecV1=StVecV;

%% Sector 2 data
% In this code we use the data from the flattened retina to create induce
% the vector field on the flattened retina.

s=linspace(0,M,ntemp+1);
s=linspace(s(2),M,ntemp);
th=linspace(a2,a3,ntemp);
[S,TH]=meshgrid(s,th);

%Construction of domain in visual optic space. This corresponds to the
%the coordinates on the hemisphere about the optical axis which the eye can
%see.

SV=-VizAngle*RV/(2*M)*S+pi*RV;
THV=TH+pi;

%Calculation of the vector field of rotations in the coordinates svis and
%thvis. The variables (t,phi) denote geodesic polar coordinate around the axis of
%rotation.

t=zeros(n,n);
phi=zeros(n,n);
vecS=zeros(n,n);
vecTH=zeros(n,n);


%This loop constructs the components of the vector field of rotations
%in the visual space of the eye.  The function InverseSTH determines the
%coordinates of (t,phi) given coordiantes svis and thvis and is necessary to
%calculate the components of the vector field.

for i=1:ntemp
    for j=1:ntemp,
        [t(i,j),phi(i,j)]=InverseSTH(alphabar,beta,RV,SV(i,j),THV(i,j));
        
        vecSV(i,j)=-cos(SV(i,j)/RV)*cos(beta)*sin(phi(i,j))*cos(THV(i,j)-alphabar)+...
            cos(SV(i,j)/RV)*cos(phi(i,j))*(sin(THV(i,j)-alphabar))-...
            sin(SV(i,j)/RV)*sin(phi(i,j))*sin(beta);
        
        vecSV(i,j)=RV*sin(t(i,j)/RV)*vecSV(i,j);
        
        vecTHV(i,j)=(cos(beta)*sin(phi(i,j))*sin(THV(i,j)-alphabar)+...
            cos(phi(i,j))*cos(alphabar-THV(i,j)));
        
        vecTHV(i,j)=RV*sin(t(i,j)/RV)*vecTHV(i,j);
        vecTHV(i,j)=1/(RV*sin(SV(i,j)/RV))*vecTHV(i,j);
        
    end
end

%Conversion of the vector fields in visual space to vector fields on the
%retina.

vecS=-2*M/(VizAngle*RV)*vecSV;
vecTH=vecTHV;

%We now push forward this vector onto the flattened retina.

StSec2Data=dlmread('StandardSector2_40');

StRHO=StSec2Data(1:ntemp,1:ntemp);
StF=StSec2Data(1:ntemp,ntemp+1:2*ntemp);
StS=StSec2Data(1:ntemp,2*ntemp+1:3*ntemp);
StTH=StSec2Data(1:ntemp,3*ntemp+1:4*ntemp);

StU=StRHO.*cos(StF);
StV=StRHO.*sin(StF);

%Calcution of derivatives
ds=StS(1,2)-StS(1,1);
DS=DiffX(ntemp)/ds;
StRHOS=(DS*StRHO')';
StFS=(DS*StF')';

dth=StTH(2,1)-StTH(1,1);
DTH=DiffX(ntemp)/dth;
StRHOTH=(DTH*StRHO);
StFTH=(DTH*StF);

StVecU=vecS.*(StRHOS.*cos(StF)-StRHO.*StFS.*sin(StF))+vecTH.*(StRHOTH.*cos(StF)-StRHO.*StFTH.*sin(StF));
StVecV=vecS.*(StRHOS.*sin(StF)+StRHO.*StFS.*cos(StF))+vecTH.*(StRHOTH.*sin(StF)+StRHO.*StFTH.*cos(StF));

temp=(StVecU.^2+StVecV.^2).^(1/2);

StVecU=StVecU./temp;
StVecV=StVecV./temp;

S2=S;
TH2=TH;
StU2=StU;
StV2=StV;
StVecU2=StVecU;
StVecV2=StVecV;
%% Sector 3 data
% In this code we use the data from the flattened retina to create induce
% the vector field on the flattened retina.

s=linspace(0,M,ntemp+1);
s=linspace(s(2),M,ntemp);
th=linspace(a3,a4,ntemp);
[S,TH]=meshgrid(s,th);

%Construction of domain in visual optic space. This corresponds to the
%the coordinates on the hemisphere about the optical axis which the eye can
%see.

SV=-VizAngle*RV/(2*M)*S+pi*RV;
THV=TH+pi;

%Calculation of the vector field of rotations in the coordinates svis and
%thvis. The variables (t,phi) denote geodesic polar coordinate around the axis of
%rotation.

t=zeros(n,n);
phi=zeros(n,n);
vecS=zeros(n,n);
vecTH=zeros(n,n);


%This loop constructs the components of the vector field of rotations
%in the visual space of the eye.  The function InverseSTH determines the
%coordinates of (t,phi) given coordiantes svis and thvis and is necessary to
%calculate the components of the vector field.

for i=1:ntemp
    for j=1:ntemp,
        [t(i,j),phi(i,j)]=InverseSTH(alphabar,beta,RV,SV(i,j),THV(i,j));
        
        vecSV(i,j)=-cos(SV(i,j)/RV)*cos(beta)*sin(phi(i,j))*cos(THV(i,j)-alphabar)+...
            cos(SV(i,j)/RV)*cos(phi(i,j))*(sin(THV(i,j)-alphabar))-...
            sin(SV(i,j)/RV)*sin(phi(i,j))*sin(beta);
        
        vecSV(i,j)=RV*sin(t(i,j)/RV)*vecSV(i,j);
        
        vecTHV(i,j)=(cos(beta)*sin(phi(i,j))*sin(THV(i,j)-alphabar)+...
            cos(phi(i,j))*cos(alphabar-THV(i,j)));
        
        vecTHV(i,j)=RV*sin(t(i,j)/RV)*vecTHV(i,j);
        vecTHV(i,j)=1/(RV*sin(SV(i,j)/RV))*vecTHV(i,j);
        
    end
end

%Conversion of the vector fields in visual space to vector fields on the
%retina.

vecS=-2*M/(VizAngle*RV)*vecSV;
vecTH=vecTHV;

%We now push forward this vector onto the flattened retina.

StSec3Data=dlmread('StandardSector3_40');


StRHO=StSec3Data(1:ntemp,1:ntemp);
StF=StSec3Data(1:ntemp,ntemp+1:2*ntemp);
StS=StSec3Data(1:ntemp,2*ntemp+1:3*ntemp);
StTH=StSec3Data(1:ntemp,3*ntemp+1:4*ntemp);

StU=StRHO.*cos(StF);
StV=StRHO.*sin(StF);

%Calcution of derivatives
ds=StS(1,2)-StS(1,1);
DS=DiffX(ntemp)/ds;
StRHOS=(DS*StRHO')';
StFS=(DS*StF')';

dth=StTH(2,1)-StTH(1,1);
DTH=DiffX(ntemp)/dth;
StRHOTH=(DTH*StRHO);
StFTH=(DTH*StF);

StVecU=vecS.*(StRHOS.*cos(StF)-StRHO.*StFS.*sin(StF))+vecTH.*(StRHOTH.*cos(StF)-StRHO.*StFTH.*sin(StF));
StVecV=vecS.*(StRHOS.*sin(StF)+StRHO.*StFS.*cos(StF))+vecTH.*(StRHOTH.*sin(StF)+StRHO.*StFTH.*cos(StF));

temp=(StVecU.^2+StVecV.^2).^(1/2);

StVecU=StVecU./temp;
StVecV=StVecV./temp;

S3=S;
TH3=TH;
StU3=StU;
StV3=StV;
StVecU3=StVecU;
StVecV3=StVecV;

%% Sector 4 data
% In this code we use the data from the flattened retina to create induce
% the vector field on the flattened retina.

s=linspace(0,M,ntemp+1);
s=linspace(s(2),M,ntemp);
th=linspace(a4,a1+2*pi,ntemp);
[S,TH]=meshgrid(s,th);

%Construction of domain in visual optic space. This corresponds to the
%the coordinates on the hemisphere about the optical axis which the eye can
%see.

SV=-VizAngle*RV/(2*M)*S+pi*RV;
THV=TH+pi;

%Calculation of the vector field of rotations in the coordinates svis and
%thvis. The variables (t,phi) denote geodesic polar coordinate around the axis of
%rotation.

t=zeros(n,n);
phi=zeros(n,n);
vecS=zeros(n,n);
vecTH=zeros(n,n);


%This loop constructs the components of the vector field of rotations
%in the visual space of the eye.  The function InverseSTH determines the
%coordinates of (t,phi) given coordiantes svis and thvis and is necessary to
%calculate the components of the vector field.

for i=1:ntemp
    for j=1:ntemp,
        [t(i,j),phi(i,j)]=InverseSTH(alphabar,beta,RV,SV(i,j),THV(i,j));
        
        vecSV(i,j)=-cos(SV(i,j)/RV)*cos(beta)*sin(phi(i,j))*cos(THV(i,j)-alphabar)+...
            cos(SV(i,j)/RV)*cos(phi(i,j))*(sin(THV(i,j)-alphabar))-...
            sin(SV(i,j)/RV)*sin(phi(i,j))*sin(beta);
        
        vecSV(i,j)=RV*sin(t(i,j)/RV)*vecSV(i,j);
        
        vecTHV(i,j)=(cos(beta)*sin(phi(i,j))*sin(THV(i,j)-alphabar)+...
            cos(phi(i,j))*cos(alphabar-THV(i,j)));
        
        vecTHV(i,j)=RV*sin(t(i,j)/RV)*vecTHV(i,j);
        vecTHV(i,j)=1/(RV*sin(SV(i,j)/RV))*vecTHV(i,j);
        
    end
end

%Conversion of the vector fields in visual space to vector fields on the
%retina.

vecS=-2*M/(VizAngle*RV)*vecSV;
vecTH=vecTHV;

%We now push forward this vector onto the flattened retina.

StSec4Data=dlmread('StandardSector4_40');


StRHO=StSec4Data(1:ntemp,1:ntemp);
StF=StSec4Data(1:ntemp,ntemp+1:2*ntemp);
StS=StSec4Data(1:ntemp,2*ntemp+1:3*ntemp);
StTH=StSec4Data(1:ntemp,3*ntemp+1:4*ntemp);

StU=StRHO.*cos(StF);
StV=StRHO.*sin(StF);

%Calcution of derivatives
ds=StS(1,2)-StS(1,1);
DS=DiffX(ntemp)/ds;
StRHOS=(DS*StRHO')';
StFS=(DS*StF')';

dth=StTH(2,1)-StTH(1,1);
DTH=DiffX(ntemp)/dth;
StRHOTH=(DTH*StRHO);
StFTH=(DTH*StF);

StVecU=vecS.*(StRHOS.*cos(StF)-StRHO.*StFS.*sin(StF))+vecTH.*(StRHOTH.*cos(StF)-StRHO.*StFTH.*sin(StF));
StVecV=vecS.*(StRHOS.*sin(StF)+StRHO.*StFS.*cos(StF))+vecTH.*(StRHOTH.*sin(StF)+StRHO.*StFTH.*cos(StF));

temp=(StVecU.^2+StVecV.^2).^(1/2);

StVecU=StVecU./temp;
StVecV=StVecV./temp;

S4=S;
TH4=TH;
StU4=StU;
StV4=StV;
StVecU4=StVecU;
StVecV4=StVecV;

%% Construction of curves on flattened retina
% This part of the code creates the curves on the flattened retina. The
% idea is that we simply find the closet point on the discrete mesh to the
% points on the curve.

CF1U=zeros(1,n);
CF1V=zeros(1,n);
CF2U=zeros(1,n);
CF2V=zeros(1,n);
CF3U=zeros(1,n);
CF3V=zeros(1,n);
CF4U=zeros(1,n);
CF4V=zeros(1,n);
CF5U=zeros(1,n);
CF5V=zeros(1,n);
CF6U=zeros(1,n);
CF6V=zeros(1,n);
CF7U=zeros(1,n);
CF7V=zeros(1,n);
CF8U=zeros(1,n);
CF8V=zeros(1,n);

for i=1:n
    
    if isnan(sret1(i)),
        CF1U(i)=NaN;
        CF1V(i)=NaN;
        
        CF1TU(i)=NaN;
        CF1TV(i)=NaN;
        
    elseif thret1(i)>=a1 & thret1(i)<a2,
        %Calculation of closest grid point.
        CF1U(i)=interp2(M*S1,TH1,StU1,sret1(i),thret1(i));
        CF1V(i)=interp2(M*S1,TH1,StV1,sret1(i),thret1(i));
        
        CF1TU(i)=interp2(M*S1,TH1,StVecU1,sret1(i),thret1(i));
        CF1TV(i)=interp2(M*S1,TH1,StVecV1,sret1(i),thret1(i));
        
        
    elseif thret1(i)>=a2 & thret1(i)<a3,
        %Calculation of closest grid point.
        CF1U(i)=interp2(M*S2,TH2,StU2,sret1(i),thret1(i));
        CF1V(i)=interp2(M*S2,TH2,StV2,sret1(i),thret1(i));
        
        CF1TU(i)=interp2(M*S2,TH2,StVecU2,sret1(i),thret1(i));
        CF1TV(i)=interp2(M*S2,TH2,StVecV2,sret1(i),thret1(i));
        
    elseif thret1(i)>=a3 & thret1(i)<a4,
        %Calculation of closest grid point.
        CF1U(i)=interp2(M*S3,TH3,StU3,sret1(i),thret1(i));
        CF1V(i)=interp2(M*S3,TH3,StV3,sret1(i),thret1(i));
        
        CF1TU(i)=interp2(M*S3,TH3,StVecU3,sret1(i),thret1(i));
        CF1TV(i)=interp2(M*S3,TH3,StVecV3,sret1(i),thret1(i));
        
    else
        if thret1(i)<pi/2,
            thret1(i)=thret1(i)+2*pi;
        end
        %Calculation of closest grid point.
        CF1U(i)=interp2(M*S4,TH4,StU4,sret1(i),thret1(i));
        CF1V(i)=interp2(M*S4,TH4,StV4,sret1(i),thret1(i));
        
        CF1TU(i)=interp2(M*S4,TH4,StVecU4,sret1(i),thret1(i));
        CF1TV(i)=interp2(M*S4,TH4,StVecV4,sret1(i),thret1(i));
    end
    
    if isnan(sret2(i)),
        CF2U(i)=NaN;
        CF2V(i)=NaN;
        
        CF2TU(i)=NaN;
        CF2TV(i)=NaN;
        
    elseif thret2(i)>=a1 & thret2(i)<a2,
        %Calculation of closest grid point.
        CF2U(i)=interp2(M*S1,TH1,StU1,sret2(i),thret2(i));
        CF2V(i)=interp2(M*S1,TH1,StV1,sret2(i),thret2(i));
        
        CF2TU(i)=interp2(M*S3,TH3,StVecU3,sret2(i),thret2(i));
        CF2TV(i)=interp2(M*S3,TH3,StVecV3,sret2(i),thret2(i));
        
    elseif thret2(i)>=a2 & thret2(i)<a3,
        %Calculation of closest grid point.
        CF2U(i)=interp2(M*S2,TH2,StU2,sret2(i),thret2(i));
        CF2V(i)=interp2(M*S2,TH2,StV2,sret2(i),thret2(i));
        
        CF2TU(i)=interp2(M*S2,TH2,StVecU2,sret2(i),thret2(i));
        CF2TV(i)=interp2(M*S2,TH2,StVecV2,sret2(i),thret2(i));
        
    elseif thret2(i)>=a3 & thret2(i)<a4,
        %Calculation of closest grid point.
        CF2U(i)=interp2(M*S3,TH3,StU3,sret2(i),thret2(i));
        CF2V(i)=interp2(M*S3,TH3,StV3,sret2(i),thret2(i));
        
        CF2TU(i)=interp2(M*S3,TH3,StVecU3,sret2(i),thret2(i));
        CF2TV(i)=interp2(M*S3,TH3,StVecV3,sret2(i),thret2(i));
        
    else
        if thret2(i)<pi/2,
            thret2(i)=thret2(i)+2*pi;
        end
        %Calculation of closest grid point.
        CF2U(i)=interp2(M*S4,TH4,StU4,sret2(i),thret2(i));
        CF2V(i)=interp2(M*S4,TH4,StV4,sret2(i),thret2(i));
        
        CF2TU(i)=interp2(M*S4,TH4,StVecU4,sret2(i),thret2(i));
        CF2TV(i)=interp2(M*S4,TH4,StVecV4,sret2(i),thret2(i));
    end
    
    if isnan(sret3(i)),
        CF3U(i)=NaN;
        CF3V(i)=NaN;
        
        CF3TU(i)=NaN;
        CF3TV(i)=NaN;
    elseif thret3(i)>=a1 & thret3(i)<a2,
        %Calculation of closest grid point.
        CF3U(i)=interp2(M*S1,TH1,StU1,sret3(i),thret3(i));
        CF3V(i)=interp2(M*S1,TH1,StV1,sret3(i),thret3(i));
        
        CF3TU(i)=interp2(M*S1,TH1,StVecU1,sret3(i),thret3(i));
        CF3TV(i)=interp2(M*S1,TH1,StVecV1,sret3(i),thret3(i));
    elseif thret3(i)>=a2 & thret3(i)<a3,
        %Calculation of closest grid point.
        CF3U(i)=interp2(M*S2,TH2,StU2,sret3(i),thret3(i));
        CF3V(i)=interp2(M*S2,TH2,StV2,sret3(i),thret3(i));
        
        CF3TU(i)=interp2(M*S2,TH2,StVecU2,sret3(i),thret3(i));
        CF3TV(i)=interp2(M*S2,TH2,StVecV2,sret3(i),thret3(i));
        
    elseif thret3(i)>=a3 & thret3(i)<a4,
        %Calculation of closest grid point.
        CF3U(i)=interp2(M*S3,TH3,StU3,sret3(i),thret3(i));
        CF3V(i)=interp2(M*S3,TH3,StV3,sret3(i),thret3(i));
        
        CF3TU(i)=interp2(M*S3,TH3,StVecU3,sret3(i),thret3(i));
        CF3TV(i)=interp2(M*S3,TH3,StVecV3,sret3(i),thret3(i));
        
    else
        if thret3(i)<pi/2,
            thret3(i)=thret3(i)+2*pi;
        end
        %Calculation of closest grid point.
        CF3U(i)=interp2(M*S4,TH4,StU4,sret3(i),thret3(i));
        CF3V(i)=interp2(M*S4,TH4,StV4,sret3(i),thret3(i));
        
        CF3TU(i)=interp2(M*S4,TH4,StVecU4,sret3(i),thret3(i));
        CF3TV(i)=interp2(M*S4,TH4,StVecV4,sret3(i),thret3(i));
    end
    
    if isnan(sret4(i)),
        CF4U(i)=NaN;
        CF4V(i)=NaN;
        
        CF4TU(i)=NaN;
        CF4TV(i)=NaN;
    elseif thret4(i)>=a1 & thret4(i)<a2,
        %Calculation of closest grid point.
        CF4U(i)=interp2(M*S1,TH1,StU1,sret4(i),thret4(i));
        CF4V(i)=interp2(M*S1,TH1,StV1,sret4(i),thret4(i));
        
        CF4TU(i)=interp2(M*S1,TH1,StVecU1,sret4(i),thret4(i));
        CF4TV(i)=interp2(M*S1,TH1,StVecV1,sret4(i),thret4(i));
        
    elseif thret4(i)>=a2 & thret4(i)<a3,
        %Calculation of closest grid point.
        CF4U(i)=interp2(M*S2,TH2,StU2,sret4(i),thret4(i));
        CF4V(i)=interp2(M*S2,TH2,StV2,sret4(i),thret4(i));
        
        CF4TU(i)=interp2(M*S2,TH2,StVecU2,sret4(i),thret4(i));
        CF4TV(i)=interp2(M*S2,TH2,StVecV2,sret4(i),thret4(i));
        
    elseif thret4(i)>=a3 & thret4(i)<a4,
        %Calculation of closest grid point.
        CF4U(i)=interp2(M*S3,TH3,StU3,sret4(i),thret4(i));
        CF4V(i)=interp2(M*S3,TH3,StV3,sret4(i),thret4(i));
        
        CF4TU(i)=interp2(M*S3,TH3,StVecU3,sret4(i),thret4(i));
        CF4TV(i)=interp2(M*S3,TH3,StVecV3,sret4(i),thret4(i));
        
    else
        if thret4(i)<pi/2,
            thret4(i)=thret4(i)+2*pi;
        end
        %Calculation of closest grid point.
        CF4U(i)=interp2(M*S4,TH4,StU4,sret4(i),thret4(i));
        CF4V(i)=interp2(M*S4,TH4,StV4,sret4(i),thret4(i));
        
        CF4TU(i)=interp2(M*S4,TH4,StVecU4,sret4(i),thret4(i));
        CF4TV(i)=interp2(M*S4,TH4,StVecV4,sret4(i),thret4(i));
    end
    
    if isnan(sret5(i)),
        CF5U(i)=NaN;
        CF5V(i)=NaN;
        
        CF5TU(i)=NaN;
        CF5TV(i)=NaN;
        
    elseif thret5(i)>=a1 & thret5(i)<a2,
        %Calculation of closest grid point.
        CF5U(i)=interp2(M*S1,TH1,StU1,sret5(i),thret5(i));
        CF5V(i)=interp2(M*S1,TH1,StV1,sret5(i),thret5(i));
        
        CF5TU(i)=interp2(M*S1,TH1,StVecU1,sret5(i),thret5(i));
        CF5TV(i)=interp2(M*S1,TH1,StVecV1,sret5(i),thret5(i));
    elseif thret5(i)>=a2 & thret5(i)<a3,
        %Calculation of closest grid point.
        CF5U(i)=interp2(M*S2,TH2,StU2,sret5(i),thret5(i));
        CF5V(i)=interp2(M*S2,TH2,StV2,sret5(i),thret5(i));
        
        CF5TU(i)=interp2(M*S2,TH2,StVecU2,sret5(i),thret5(i));
        CF5TV(i)=interp2(M*S2,TH2,StVecV2,sret5(i),thret5(i));
        
    elseif thret5(i)>=a3 & thret5(i)<a4,
        %Calculation of closest grid point.
        CF5U(i)=interp2(M*S3,TH3,StU3,sret5(i),thret5(i));
        CF5V(i)=interp2(M*S3,TH3,StV3,sret5(i),thret5(i));
        
        CF5TU(i)=interp2(M*S3,TH3,StVecU3,sret5(i),thret5(i));
        CF5TV(i)=interp2(M*S3,TH3,StVecV3,sret5(i),thret5(i));
        
    else
        if thret5(i)<pi/2,
            thret5(i)=thret5(i)+2*pi;
        end
        %Calculation of closest grid point.
        CF5U(i)=interp2(M*S4,TH4,StU4,sret5(i),thret5(i));
        CF5V(i)=interp2(M*S4,TH4,StV4,sret5(i),thret5(i));
        
        CF5TU(i)=interp2(M*S4,TH4,StVecU4,sret5(i),thret5(i));
        CF5TV(i)=interp2(M*S4,TH4,StVecV4,sret5(i),thret5(i));
    end
    
    if isnan(sret6(i)),
        CF6U(i)=NaN;
        CF6V(i)=NaN;
        
        CF6TU=NaN;
        CF6TV=NaN;
        
    elseif thret6(i)>=a1 & thret6(i)<a2,
        %Calculation of closest grid point.
        CF6U(i)=interp2(M*S1,TH1,StU1,sret6(i),thret6(i));
        CF6V(i)=interp2(M*S1,TH1,StV1,sret6(i),thret6(i));
        
        CF6TU(i)=interp2(M*S1,TH1,StVecU1,sret6(i),thret6(i));
        CF6TV(i)=interp2(M*S1,TH1,StVecV1,sret6(i),thret6(i));
    elseif thret6(i)>=a2 & thret6(i)<a3,
        %Calculation of closest grid point.
        CF6U(i)=interp2(M*S2,TH2,StU2,sret6(i),thret6(i));
        CF6V(i)=interp2(M*S2,TH2,StV2,sret6(i),thret6(i));
        
        CF6TU(i)=interp2(M*S2,TH2,StVecU2,sret6(i),thret6(i));
        CF6TV(i)=interp2(M*S2,TH2,StVecV2,sret6(i),thret6(i));
        
    elseif thret6(i)>=a3 & thret6(i)<a4,
        %Calculation of closest grid point.
        CF6U(i)=interp2(M*S3,TH3,StU3,sret6(i),thret6(i));
        CF6V(i)=interp2(M*S3,TH3,StV3,sret6(i),thret6(i));
        
        CF6TU(i)=interp2(M*S2,TH2,StVecU3,sret6(i),thret6(i));
        CF6TV(i)=interp2(M*S2,TH2,StVecV3,sret6(i),thret6(i));
        
    else
        if thret6(i)<pi/2,
            thret6(i)=thret6(i)+2*pi;
        end
        %Calculation of closest grid point.
        CF6U(i)=interp2(M*S4,TH4,StU4,sret6(i),thret6(i));
        CF6V(i)=interp2(M*S4,TH4,StV4,sret6(i),thret6(i));
        
        CF6TU(i)=interp2(M*S4,TH4,StVecU4,sret6(i),thret6(i));
        CF6TV(i)=interp2(M*S4,TH4,StVecV4,sret6(i),thret6(i));
    end
    
    if isnan(sret7(i)),
        CF7U(i)=NaN;
        CF7V(i)=NaN;
        
        CF7TU(i)=NaN;
        CF7TV(i)=NaN;
        
    elseif thret7(i)>=a1 & thret7(i)<a2,
        %Calculation of closest grid point.
        CF7U(i)=interp2(M*S1,TH1,StU1,sret7(i),thret7(i));
        CF7V(i)=interp2(M*S1,TH1,StV1,sret7(i),thret7(i));
        
        CF7TU(i)=interp2(M*S1,TH1,StVecU1,sret7(i),thret7(i));
        CF7TV(i)=interp2(M*S1,TH1,StVecV1,sret7(i),thret7(i));
        
    elseif thret7(i)>=a2 & thret7(i)<a3,
        %Calculation of closest grid point.
        CF7U(i)=interp2(M*S2,TH2,StU2,sret7(i),thret7(i));
        CF7V(i)=interp2(M*S2,TH2,StV2,sret7(i),thret7(i));
        
        CF7TU(i)=interp2(M*S2,TH2,StVecU2,sret7(i),thret7(i));
        CF7TV(i)=interp2(M*S2,TH2,StVecV2,sret7(i),thret7(i));
        
    elseif thret7(i)>=a3 & thret7(i)<a4,
        %Calculation of closest grid point.
        CF7U(i)=interp2(M*S3,TH3,StU3,sret7(i),thret7(i));
        CF7V(i)=interp2(M*S3,TH3,StV3,sret7(i),thret7(i));
        
        CF7TU(i)=interp2(M*S2,TH2,StVecU3,sret7(i),thret7(i));
        CF7TV(i)=interp2(M*S2,TH2,StVecV3,sret7(i),thret7(i));
        
    else
        if thret7(i)<pi/2,
            thret7(i)=thret7(i)+2*pi;
        end
        %Calculation of closest grid point.
        CF7U(i)=interp2(M*S4,TH4,StU4,sret7(i),thret7(i));
        CF7V(i)=interp2(M*S4,TH4,StV4,sret7(i),thret7(i));
        
        CF7TU(i)=interp2(M*S4,TH4,StVecU4,sret7(i),thret7(i));
        CF7TV(i)=interp2(M*S4,TH4,StVecV4,sret7(i),thret7(i));
    end
    
    if isnan(sret8(i)),
        CF8U(i)=NaN;
        CF8V(i)=NaN;
        
        CF8TU(i)=NaN;
        CF8TV(i)=NaN;
        
    elseif thret8(i)>=a1 & thret8(i)<a2,
        %Calculation of closest grid point.
        CF8U(i)=interp2(M*S1,TH1,StU1,sret8(i),thret8(i));
        CF8V(i)=interp2(M*S1,TH1,StV1,sret8(i),thret8(i));
        
        CF8TU(i)=interp2(M*S1,TH1,StVecU1,sret8(i),thret8(i));
        CF8TV(i)=interp2(M*S1,TH1,StVecV1,sret8(i),thret8(i));
    elseif thret8(i)>=a2 & thret8(i)<a3,
        %Calculation of closest grid point.
        CF8U(i)=interp2(M*S2,TH2,StU2,sret8(i),thret8(i));
        CF8V(i)=interp2(M*S2,TH2,StV2,sret8(i),thret8(i));
        
        CF8TU(i)=interp2(M*S2,TH2,StVecU2,sret8(i),thret8(i));
        CF8TV(i)=interp2(M*S2,TH2,StVecV2,sret8(i),thret8(i));
        
    elseif thret8(i)>=a3 & thret8(i)<a4,
        %Calculation of closest grid point.
        CF8U(i)=interp2(M*S3,TH3,StU3,sret8(i),thret8(i));
        CF8V(i)=interp2(M*S3,TH3,StV3,sret8(i),thret8(i));
        
        CF8TU(i)=interp2(M*S3,TH3,StVecU3,sret8(i),thret8(i));
        CF8TV(i)=interp2(M*S3,TH3,StVecV3,sret8(i),thret8(i));
    else
        if thret8(i)<pi/2,
            thret8(i)=thret8(i)+2*pi;
        end
        %Calculation of closest grid point.
        CF8U(i)=interp2(M*S4,TH4,StU4,sret8(i),thret8(i));
        CF8V(i)=interp2(M*S4,TH4,StV4,sret8(i),thret8(i));
        
        CF8TU(i)=interp2(M*S4,TH4,StVecU4,sret8(i),thret8(i));
        CF8TV(i)=interp2(M*S4,TH4,StVecV4,sret8(i),thret8(i));
    end
    
    if isnan(sret9(i)),
        CF9U(i)=NaN;
        CF9V(i)=NaN;
        
        CF9TU(i)=NaN;
        CF9TV(i)=NaN;
        
    elseif thret9(i)>=a1 & thret9(i)<a2,
        %Calculation of closest grid point.
        CF9U(i)=interp2(M*S1,TH1,StU1,sret9(i),thret9(i));
        CF9V(i)=interp2(M*S1,TH1,StV1,sret9(i),thret9(i));
        
        CF9TU(i)=interp2(M*S1,TH1,StVecU1,sret9(i),thret9(i));
        CF9TV(i)=interp2(M*S1,TH1,StVecV1,sret9(i),thret9(i));
    elseif thret9(i)>=a2 & thret9(i)<a3,
        %Calculation of closest grid point.
        CF9U(i)=interp2(M*S2,TH2,StU2,sret9(i),thret9(i));
        CF9V(i)=interp2(M*S2,TH2,StV2,sret9(i),thret9(i));
        
        CF9TU(i)=interp2(M*S2,TH2,StVecU2,sret9(i),thret9(i));
        CF9TV(i)=interp2(M*S2,TH2,StVecV2,sret9(i),thret9(i));
        
    elseif thret9(i)>=a3 & thret9(i)<a4,
        %Calculation of closest grid point.
        CF9U(i)=interp2(M*S3,TH3,StU3,sret9(i),thret9(i));
        CF9V(i)=interp2(M*S3,TH3,StV3,sret9(i),thret9(i));
        
        CF9TU(i)=interp2(M*S3,TH3,StVecU3,sret9(i),thret9(i));
        CF9TV(i)=interp2(M*S3,TH3,StVecV3,sret9(i),thret9(i));
        
    else
        if thret9(i)<pi/2,
            thret9(i)=thret9(i)+2*pi;
        end
        %Calculation of closest grid point.
        CF9U(i)=interp2(M*S4,TH4,StU4,sret9(i),thret9(i));
        CF9V(i)=interp2(M*S4,TH4,StV4,sret9(i),thret9(i));
        
        CF9TU(i)=interp2(M*S4,TH4,StVecU4,sret9(i),thret9(i));
        CF9TV(i)=interp2(M*S4,TH4,StVecV4,sret9(i),thret9(i));
    end
    
    if isnan(sret10(i)),
        CF10U(i)=NaN;
        CF10V(i)=NaN;
        
        CF10TU(i)=NaN;
        CF10TV(i)=NaN;
        
    elseif thret10(i)>=a1 & thret10(i)<a2,
        %Calculation of closest grid point.
        CF10U(i)=interp2(M*S1,TH1,StU1,sret10(i),thret10(i));
        CF10V(i)=interp2(M*S1,TH1,StV1,sret10(i),thret10(i));
        
        CF10TU(i)=interp2(M*S1,TH1,StVecU1,sret10(i),thret10(i));
        CF10TV(i)=interp2(M*S1,TH1,StVecV1,sret10(i),thret10(i));
    elseif thret10(i)>=a2 & thret10(i)<a3,
        %Calculation of closest grid point.
        CF10U(i)=interp2(M*S2,TH2,StU2,sret10(i),thret10(i));
        CF10V(i)=interp2(M*S2,TH2,StV2,sret10(i),thret10(i));
        
        CF10TU(i)=interp2(M*S2,TH2,StVecU2,sret10(i),thret10(i));
        CF10TV(i)=interp2(M*S2,TH2,StVecV2,sret10(i),thret10(i));
        
    elseif thret10(i)>=a3 & thret10(i)<a4,
        %Calculation of closest grid point.
        CF10U(i)=interp2(M*S3,TH3,StU3,sret10(i),thret10(i));
        CF10V(i)=interp2(M*S3,TH3,StV3,sret10(i),thret10(i));
        
        CF10TU(i)=interp2(M*S3,TH3,StVecU3,sret10(i),thret10(i));
        CF10TV(i)=interp2(M*S3,TH3,StVecV3,sret10(i),thret10(i));
        
    else
        if thret10(i)<pi/2,
            thret10(i)=thret10(i)+2*pi;
        end
        %Calculation of closest grid point.
        CF10U(i)=interp2(M*S4,TH4,StU4,sret10(i),thret10(i));
        CF10V(i)=interp2(M*S4,TH4,StV4,sret10(i),thret10(i));
        
        CF10TU(i)=interp2(M*S4,TH4,StVecU4,sret10(i),thret10(i));
        CF10TV(i)=interp2(M*S4,TH4,StVecV4,sret10(i),thret10(i));
    end
    
    
end

%% Plotting Curves on flattened retina

figure('Position',get(0,'ScreenSize'));
hold on;

%Plotting Flattened Retina
surf(StU1,StV1,zeros(ntemp,ntemp),'FaceColor','red','EdgeColor','none');
surf(StU2,StV2,zeros(ntemp,ntemp),'FaceColor','red','EdgeColor','none');
surf(StU3,StV3,zeros(ntemp,ntemp),'FaceColor','red','EdgeColor','none');
surf(StU4,StV4,zeros(ntemp,ntemp),'FaceColor','red','EdgeColor','none');
alpha(.25);


%Plotting Tangents on flattened retina
step=25;
scale=1;

quiver(CF1U(1:step:n),CF1V(1:step:n),CF1TU(1:step:n),CF1TV(1:step:n),scale,'b','LineWidth',2);

scale=.5;
quiver(CF2U(1:step:n),CF2V(1:step:n),CF2TU(1:step:n),CF2TV(1:step:n),scale,'b','LineWidth',2);
quiver(CF3U(1:step:n),CF3V(1:step:n),CF3TU(1:step:n),CF3TV(1:step:n),scale,'b','LineWidth',2);
quiver(CF4U(1:step:n),CF4V(1:step:n),CF4TU(1:step:n),CF4TV(1:step:n),scale,'b','LineWidth',2);
quiver(CF5U(1:step:n),CF5V(1:step:n),CF5TU(1:step:n),CF5TV(1:step:n),scale,'b','LineWidth',2);
quiver(CF6U(1:step:n),CF6V(1:step:n),CF6TU(1:step:n),CF6TV(1:step:n),scale,'b','LineWidth',2);
quiver(CF7U(1:step:n),CF7V(1:step:n),CF7TU(1:step:n),CF7TV(1:step:n),scale,'b','LineWidth',2);
quiver(CF8U(1:step:n),CF8V(1:step:n),CF8TU(1:step:n),CF8TV(1:step:n),scale,'b','LineWidth',2);
quiver(CF9U(1:step:n),CF9V(1:step:n),CF9TU(1:step:n),CF9TV(1:step:n),scale,'b','LineWidth',2);

scale=1;
quiver(CF10U(1:step:n),CF10V(1:step:n),CF10TU(1:step:n),CF10TV(1:step:n),scale,'b','LineWidth',2);

%Plotting Curves on Retina
plot(CF1U,CF1V,'k','LineWidth',2);
plot(CF2U,CF2V,'k','LineWidth',2);
plot(CF3U,CF3V,'k','LineWidth',2);
plot(CF4U,CF4V,'k','LineWidth',2);
plot(CF5U,CF5V,'k','LineWidth',2);
plot(CF6U,CF6V,'k','LineWidth',2);
plot(CF7U,CF7V,'k','LineWidth',2);
plot(CF8U,CF8V,'k','LineWidth',2);
plot(CF9U,CF9V,'k','LineWidth',2);
plot(CF10U,CF10V,'k','LineWidth',2);


axis off;
axis equal;
set(gcf, 'color', [1 1 1]);

hold off;

% %% Movie on flattened retina
% 
% rate=10;
% maxvalue=n/rate; %frame rate
% 
% AnglePar=linspace(0,2*pi,n);
% 
% for i=1:maxvalue,
%     close all
%     j=i*rate;
%     
%     figure('Position',get(0,'ScreenSize'));
%     hold on;
%     
%     %Plotting Flattened Retina
%     surf(StU1,StV1,zeros(ntemp,ntemp),'FaceColor','red','EdgeColor','none');
%     surf(StU2,StV2,zeros(ntemp,ntemp),'FaceColor','red','EdgeColor','none');
%     surf(StU3,StV3,zeros(ntemp,ntemp),'FaceColor','red','EdgeColor','none');
%     surf(StU4,StV4,zeros(ntemp,ntemp),'FaceColor','red','EdgeColor','none');
%     alpha(.25);
%     
%     
%     %     %Plotting Tangents on flattened retina
%     %     step=25;
%     %     scale=1;
%     %
%     %     quiver(CF1U(1:step:n),CF1V(1:step:n),CF1TU(1:step:n),CF1TV(1:step:n),scale,'b','LineWidth',2);
%     %
%     %     scale=.5;
%     %     quiver(CF2U(1:step:n),CF2V(1:step:n),CF2TU(1:step:n),CF2TV(1:step:n),scale,'b','LineWidth',2);
%     %     quiver(CF3U(1:step:n),CF3V(1:step:n),CF3TU(1:step:n),CF3TV(1:step:n),scale,'b','LineWidth',2);
%     %     quiver(CF4U(1:step:n),CF4V(1:step:n),CF4TU(1:step:n),CF4TV(1:step:n),scale,'b','LineWidth',2);
%     %     quiver(CF5U(1:step:n),CF5V(1:step:n),CF5TU(1:step:n),CF5TV(1:step:n),scale,'b','LineWidth',2);
%     %     quiver(CF6U(1:step:n),CF6V(1:step:n),CF6TU(1:step:n),CF6TV(1:step:n),scale,'b','LineWidth',2);
%     %     quiver(CF7U(1:step:n),CF7V(1:step:n),CF7TU(1:step:n),CF7TV(1:step:n),scale,'b','LineWidth',2);
%     %     quiver(CF8U(1:step:n),CF8V(1:step:n),CF8TU(1:step:n),CF8TV(1:step:n),scale,'b','LineWidth',2);
%     %     quiver(CF9U(1:step:n),CF9V(1:step:n),CF9TU(1:step:n),CF9TV(1:step:n),scale,'b','LineWidth',2);
%     
%     %     scale=1;
%     %     quiver(CF10U(1:step:n),CF10V(1:step:n),CF10TU(1:step:n),CF10TV(1:step:n),scale,'b','LineWidth',2);
%     
%     %Plotting Curves on Retina
%     plot(CF1U,CF1V,'k','LineWidth',2);
%     plot(CF2U,CF2V,'k','LineWidth',2);
%     plot(CF3U,CF3V,'k','LineWidth',2);
%     plot(CF4U,CF4V,'k','LineWidth',2);
%     plot(CF5U,CF5V,'k','LineWidth',2);
%     plot(CF6U,CF6V,'k','LineWidth',2);
%     plot(CF7U,CF7V,'k','LineWidth',2);
%     plot(CF8U,CF8V,'k','LineWidth',2);
%     plot(CF9U,CF9V,'k','LineWidth',2);
%     plot(CF10U,CF10V,'k','LineWidth',2);
%     
%     %Plotting Points on Retina
%     plot(CF1U(j),CF1V(j),'b.','MarkerSize',40);
%     plot(CF2U(j),CF2V(j),'b.','MarkerSize',40);
%     plot(CF3U(j),CF3V(j),'b.','MarkerSize',40);
%     plot(CF4U(j),CF4V(j),'b.','MarkerSize',40);
%     plot(CF5U(j),CF5V(j),'b.','MarkerSize',40);
%     plot(CF6U(j),CF6V(j),'b.','MarkerSize',40);
%     plot(CF7U(j),CF7V(j),'b.','MarkerSize',40);
%     plot(CF8U(j),CF8V(j),'b.','MarkerSize',40);
%     plot(CF9U(j),CF9V(j),'b.','MarkerSize',40);
%     plot(CF10U(j),CF10V(j),'b.','MarkerSize',40);
%     
%     
%     axis off;
%     axis equal;
%     set(gcf, 'color', [1 1 1]);
%     
%     F(i)=getframe(gca);
%     F(i)=getframe(gca);
%     hold off;
%     movie2avi(F,'RotationFlat.avi');
%     
%     hold off;
%     
% end
%% Plotting Vector field on flattened retina
figure('Position',get(0,'ScreenSize'));
hold on;

%Plotting Flattened Retina
surf(StU1,StV1,zeros(ntemp,ntemp),'FaceColor','red','EdgeColor','none');
surf(StU2,StV2,zeros(ntemp,ntemp),'FaceColor','red','EdgeColor','none');
surf(StU3,StV3,zeros(ntemp,ntemp),'FaceColor','red','EdgeColor','none');
surf(StU4,StV4,zeros(ntemp,ntemp),'FaceColor','red','EdgeColor','none');
alpha(.25);


%Plotting vector field
size(StU1)
scale=.7;
step=3;
quiver(StU1(1:step:ntemp^2),StV1(1:step:ntemp^2),StVecU1(1:step:ntemp^2),StVecV1(1:step:ntemp^2),scale,'b');
quiver(StU2(1:step:ntemp^2),StV2(1:step:ntemp^2),StVecU2(1:step:ntemp^2),StVecV2(1:step:ntemp^2),scale,'b');
quiver(StU3(1:step:ntemp^2),StV3(1:step:ntemp^2),StVecU3(1:step:ntemp^2),StVecV3(1:step:ntemp^2),scale,'b');
quiver(StU4(1:step:ntemp^2),StV4(1:step:ntemp^2),StVecU4(1:step:ntemp^2),StVecV4(1:step:ntemp^2),scale,'b');

axis off;
axis equal;
set(gcf, 'color', [1 1 1]);

%% Vector Field on flatened retina with curves

figure('Position',get(0,'ScreenSize'));
hold on;

%Plotting Flattened Retina
surf(StU1,StV1,zeros(ntemp,ntemp),'FaceColor','red','EdgeColor','none');
surf(StU2,StV2,zeros(ntemp,ntemp),'FaceColor','red','EdgeColor','none');
surf(StU3,StV3,zeros(ntemp,ntemp),'FaceColor','red','EdgeColor','none');
surf(StU4,StV4,zeros(ntemp,ntemp),'FaceColor','red','EdgeColor','none');
alpha(.25);


%Plotting vector field
scale=.7;
step=3;
quiver(StU1(1:step:ntemp^2),StV1(1:step:ntemp^2),StVecU1(1:step:ntemp^2),StVecV1(1:step:ntemp^2),scale,'b');
quiver(StU2(1:step:ntemp^2),StV2(1:step:ntemp^2),StVecU2(1:step:ntemp^2),StVecV2(1:step:ntemp^2),scale,'b');
quiver(StU3(1:step:ntemp^2),StV3(1:step:ntemp^2),StVecU3(1:step:ntemp^2),StVecV3(1:step:ntemp^2),scale,'b');
quiver(StU4(1:step:ntemp^2),StV4(1:step:ntemp^2),StVecU4(1:step:ntemp^2),StVecV4(1:step:ntemp^2),scale,'b');

%Plotting Curves on Retina
plot(CF1U,CF1V,'k','LineWidth',2);
plot(CF2U,CF2V,'k','LineWidth',2);
plot(CF3U,CF3V,'k','LineWidth',2);
plot(CF4U,CF4V,'k','LineWidth',2);
plot(CF5U,CF5V,'k','LineWidth',2);
plot(CF6U,CF6V,'k','LineWidth',2);
plot(CF7U,CF7V,'k','LineWidth',2);
plot(CF8U,CF8V,'k','LineWidth',2);
plot(CF9U,CF9V,'k','LineWidth',2);
plot(CF10U,CF10V,'k','LineWidth',2);

axis off;
axis equal;
set(gcf, 'color', [1 1 1]);

hold off;



%% Combined Movie
rate=10;
MaxValue=n/rate; %frame rate

AnglePar=linspace(0,2*pi,n);

M=1;

for i=1:MaxValue,
    j=i*rate;
    
    close all;
    figure('Position',get(0,'ScreenSize'));
    
    ha=tight_subplot(2,2);
    
    axes(ha(1));
    hold on;
    
    % Movie in Visual Space
    
    s=linspace((pi-VizAngle/2),pi,100);
    th=linspace(0,2*pi,100);
    
    [S,TH]=meshgrid(s,th);
    surf(RV*cos(TH).*sin(S),RV*sin(TH).*sin(S),-RV*cos(S),'EdgeColor','none');
    alpha(.3);
    
    %Plotting Optical singularity
    plot3(RV*cos(alphabar)*sin(beta),RV*sin(alphabar)*sin(beta),-RV*cos(beta),'b.','MarkerSize',40);
    
    %Plotting curve in visual space.
    plot3(C1(1,:),C1(2,:),C1(3,:),'k','LineWidth',2);
    plot3(C2(1,:),C2(2,:),C2(3,:),'k','LineWidth',2);
    plot3(C3(1,:),C3(2,:),C3(3,:),'k','LineWidth',2);
    plot3(C4(1,:),C4(2,:),C4(3,:),'k','LineWidth',2);
    plot3(C5(1,:),C5(2,:),C5(3,:),'k','LineWidth',2);
    plot3(C6(1,:),C6(2,:),C6(3,:),'k','LineWidth',2);
    plot3(C7(1,:),C7(2,:),C7(3,:),'k','LineWidth',2);
    plot3(C8(1,:),C8(2,:),C8(3,:),'k','LineWidth',2);
    plot3(C9(1,:),C9(2,:),C9(3,:),'k','LineWidth',2);
    plot3(C10(1,:),C10(2,:),C10(3,:),'k','LineWidth',2);
    
    
    %Plotting Representation of mouse eye
    s=linspace(0,pi,100);
    th=linspace(0,2*pi,100);
    [S,TH]=meshgrid(s,th);
    surf(RV/10*cos(TH).*sin(S),RV/10*sin(TH).*sin(S),-RV/10*cos(S)-RV/10,'FaceColor','green','EdgeColor','none');
    %axis([-RV,RV,-RV,RV,-2*RV/10,RV]);
    alpha(.3);
    
    %Plotting Representation of mouse lens
    s=linspace(0,pi,100);
    th=linspace(0,2*pi,100);
    [S,TH]=meshgrid(s,th);
    surf(RV/28*cos(TH).*sin(S),RV/28*sin(TH).*sin(S),-RV/28*cos(S)-RV/28,'FaceColor','blue','EdgeColor','none');
    alpha(.3);
    
    %Plotting Representation of retina
    s=linspace(0,deg2rad(105),100);
    th=linspace(0,2*pi,100);
    [S,TH]=meshgrid(s,th);
    surf((RV/12)*cos(TH).*sin(S),(RV/12)*sin(TH).*sin(S),-(RV/12)*cos(S)-RV/10,'FaceColor','red','EdgeColor','none');
    alpha(.5);
    
    %Ploting normal axis
    t=linspace(0,3/2*RV,100);
    plot3(t*Normal(1),t*Normal(2),t*Normal(3),'k','LineWidth',4);
    
    %Plotting line connecting optical singularity on retina to lens
    t=linspace(0,1,100)';
    betashift=-2*deg2rad(105)/pi*beta+2*deg2rad(105);
    alphashift=alphabar+pi;
    
    L(:,1)=t*RV/12*cos(alphashift)*sin(betashift);
    L(:,2)=t*RV/12*sin(alphashift)*sin(betashift);
    L(:,3)=t*(-RV/12*cos(betashift)-RV/10);
    plot3(L(:,1),L(:,2),L(:,3),'k','LineWidth',4);
    
    %Plotting Vector Field on retina
    %quiver3(C1(1,:),C1(2,:),C1(3,:),vecXC1,vecYC1,vecZC1);
    
    %Plotting Point on curve in visual space.
    plot3(C1(1,j),C1(2,j),C1(3,j),'b.','MarkerSize',40);
    plot3(C2(1,j),C2(2,j),C2(3,j),'b.','MarkerSize',40);
    plot3(C3(1,j),C3(2,j),C3(3,j),'b.','MarkerSize',40);
    plot3(C4(1,j),C4(2,j),C4(3,j),'b.','MarkerSize',40);
    plot3(C5(1,j),C5(2,j),C5(3,j),'b.','MarkerSize',40);
    plot3(C6(1,j),C6(2,j),C6(3,j),'b.','MarkerSize',40);
    plot3(C7(1,j),C7(2,j),C7(3,j),'b.','MarkerSize',40);
    plot3(C8(1,j),C8(2,j),C8(3,j),'b.','MarkerSize',40);
    plot3(C9(1,j),C9(2,j),C9(3,j),'b.','MarkerSize',40);
    plot3(C10(1,j),C10(2,j),C10(3,j),'b.','MarkerSize',40);
    
    axis equal;
    axis([-RV,RV,-RV,RV,-2*RV/10,RV]);
    view(90,18);
    axis off;
    set(gcf, 'color', [1 1 1])
    hold off;
    
    %Movie on Retina
    axes(ha(2));
    hold on;
    
    %Plotting Surface of retina
    s=linspace(0,M,100);
    th=linspace(0,2*pi,100);
    [S,TH]=meshgrid(s,th);
    surf(R*cos(TH).*sin(S/R),R*sin(TH).*sin(S/R),-R*cos(S/R),'FaceColor','red','EdgeColor','none');
    alpha(.5);
    camlight left;
    lighting phong;
    
    %Plotting Optical singularity on retina
    betashift=-2*M/(R*pi)*beta+2*M/R;
    alphashift=alphabar+pi;
    plot3(R*cos(alphashift)*sin(betashift),R*sin(alphashift)*sin(betashift),-R*cos(betashift),'b.','MarkerSize',40);
    
    
    %Plotting curves on retina.
    
    plot3(C1ret(1,:),C1ret(2,:),C1ret(3,:),'k','LineWidth',2);
    plot3(C2ret(1,:),C2ret(2,:),C2ret(3,:),'k','LineWidth',2);
    plot3(C3ret(1,:),C3ret(2,:),C3ret(3,:),'k','LineWidth',2);
    plot3(C4ret(1,:),C4ret(2,:),C4ret(3,:),'k','LineWidth',2);
    plot3(C5ret(1,:),C5ret(2,:),C5ret(3,:),'k','LineWidth',2);
    plot3(C6ret(1,:),C6ret(2,:),C6ret(3,:),'k','LineWidth',2);
    plot3(C7ret(1,:),C7ret(2,:),C7ret(3,:),'k','LineWidth',2);
    plot3(C8ret(1,:),C8ret(2,:),C8ret(3,:),'k','LineWidth',2);
    plot3(C9ret(1,:),C9ret(2,:),C9ret(3,:),'k','LineWidth',2);
    plot3(C10ret(1,:),C10ret(2,:),C10ret(3,:),'k','LineWidth',2);
    
    %Plotting Point on curve in visual space.
    plot3(C1ret(1,j),C1ret(2,j),C1ret(3,j),'b.','MarkerSize',40);
    plot3(C2ret(1,j),C2ret(2,j),C2ret(3,j),'b.','MarkerSize',40);
    plot3(C3ret(1,j),C3ret(2,j),C3ret(3,j),'b.','MarkerSize',40);
    plot3(C4ret(1,j),C4ret(2,j),C4ret(3,j),'b.','MarkerSize',40);
    plot3(C5ret(1,j),C5ret(2,j),C5ret(3,j),'b.','MarkerSize',40);
    plot3(C6ret(1,j),C6ret(2,j),C6ret(3,j),'b.','MarkerSize',40);
    plot3(C7ret(1,j),C7ret(2,j),C7ret(3,j),'b.','MarkerSize',40);
    plot3(C8ret(1,j),C8ret(2,j),C8ret(3,j),'b.','MarkerSize',40);
    plot3(C9ret(1,j),C9ret(2,j),C9ret(3,j),'b.','MarkerSize',40);
    plot3(C10ret(1,j),C10ret(2,j),C10ret(3,j),'b.','MarkerSize',40);
    
    axis equal;
    axis([-R,R,-R,R,-R,R/2]);
    view(90,18);
    axis off;
    set(gcf, 'color', [1 1 1])
    
    %Plotting Curves on flattened retina
    axes(ha(3));
    
    hold on;
    
    %Plotting Flattened Retina
    surf(StU1,StV1,zeros(ntemp,ntemp),'FaceColor','red','EdgeColor','none');
    surf(StU2,StV2,zeros(ntemp,ntemp),'FaceColor','red','EdgeColor','none');
    surf(StU3,StV3,zeros(ntemp,ntemp),'FaceColor','red','EdgeColor','none');
    surf(StU4,StV4,zeros(ntemp,ntemp),'FaceColor','red','EdgeColor','none');
    alpha(.25);
    
    
    %     %Plotting Tangents on flattened retina
    %     step=25;
    %     scale=1;
    %
    %     quiver(CF1U(1:step:n),CF1V(1:step:n),CF1TU(1:step:n),CF1TV(1:step:n),scale,'b','LineWidth',2);
    %
    %     scale=.5;
    %     quiver(CF2U(1:step:n),CF2V(1:step:n),CF2TU(1:step:n),CF2TV(1:step:n),scale,'b','LineWidth',2);
    %     quiver(CF3U(1:step:n),CF3V(1:step:n),CF3TU(1:step:n),CF3TV(1:step:n),scale,'b','LineWidth',2);
    %     quiver(CF4U(1:step:n),CF4V(1:step:n),CF4TU(1:step:n),CF4TV(1:step:n),scale,'b','LineWidth',2);
    %     quiver(CF5U(1:step:n),CF5V(1:step:n),CF5TU(1:step:n),CF5TV(1:step:n),scale,'b','LineWidth',2);
    %     quiver(CF6U(1:step:n),CF6V(1:step:n),CF6TU(1:step:n),CF6TV(1:step:n),scale,'b','LineWidth',2);
    %     quiver(CF7U(1:step:n),CF7V(1:step:n),CF7TU(1:step:n),CF7TV(1:step:n),scale,'b','LineWidth',2);
    %     quiver(CF8U(1:step:n),CF8V(1:step:n),CF8TU(1:step:n),CF8TV(1:step:n),scale,'b','LineWidth',2);
    %     quiver(CF9U(1:step:n),CF9V(1:step:n),CF9TU(1:step:n),CF9TV(1:step:n),scale,'b','LineWidth',2);
    
    %     scale=1;
    %     quiver(CF10U(1:step:n),CF10V(1:step:n),CF10TU(1:step:n),CF10TV(1:step:n),scale,'b','LineWidth',2);
    
    %Plotting Curves on Retina
    plot(CF1U,CF1V,'k','LineWidth',2);
    plot(CF2U,CF2V,'k','LineWidth',2);
    plot(CF3U,CF3V,'k','LineWidth',2);
    plot(CF4U,CF4V,'k','LineWidth',2);
    plot(CF5U,CF5V,'k','LineWidth',2);
    plot(CF6U,CF6V,'k','LineWidth',2);
    plot(CF7U,CF7V,'k','LineWidth',2);
    plot(CF8U,CF8V,'k','LineWidth',2);
    plot(CF9U,CF9V,'k','LineWidth',2);
    plot(CF10U,CF10V,'k','LineWidth',2);
    
    %Plotting Points on Retina
    plot(CF1U(j),CF1V(j),'b.','MarkerSize',40);
    plot(CF2U(j),CF2V(j),'b.','MarkerSize',40);
    plot(CF3U(j),CF3V(j),'b.','MarkerSize',40);
    plot(CF4U(j),CF4V(j),'b.','MarkerSize',40);
    plot(CF5U(j),CF5V(j),'b.','MarkerSize',40);
    plot(CF6U(j),CF6V(j),'b.','MarkerSize',40);
    plot(CF7U(j),CF7V(j),'b.','MarkerSize',40);
    plot(CF8U(j),CF8V(j),'b.','MarkerSize',40);
    plot(CF9U(j),CF9V(j),'b.','MarkerSize',40);
    plot(CF10U(j),CF10V(j),'b.','MarkerSize',40);
    
    hold off;
    axis off;
    set(gcf, 'color', [1 1 1]);
    

%Plotting Flattened Retina

    axes(ha(4));
    hold on;
    
surf(StU1,StV1,zeros(ntemp,ntemp),'FaceColor','red','EdgeColor','none');
surf(StU2,StV2,zeros(ntemp,ntemp),'FaceColor','red','EdgeColor','none');
surf(StU3,StV3,zeros(ntemp,ntemp),'FaceColor','red','EdgeColor','none');
surf(StU4,StV4,zeros(ntemp,ntemp),'FaceColor','red','EdgeColor','none');
alpha(.25);


%Plotting vector field
scale=.7;
step=3;
quiver(StU1(1:step:ntemp^2),StV1(1:step:ntemp^2),StVecU1(1:step:ntemp^2),StVecV1(1:step:ntemp^2),scale,'b');
quiver(StU2(1:step:ntemp^2),StV2(1:step:ntemp^2),StVecU2(1:step:ntemp^2),StVecV2(1:step:ntemp^2),scale,'b');
quiver(StU3(1:step:ntemp^2),StV3(1:step:ntemp^2),StVecU3(1:step:ntemp^2),StVecV3(1:step:ntemp^2),scale,'b');
quiver(StU4(1:step:ntemp^2),StV4(1:step:ntemp^2),StVecU4(1:step:ntemp^2),StVecV4(1:step:ntemp^2),scale,'b');

%Plotting Curves on Retina
plot(CF1U,CF1V,'k','LineWidth',2);
plot(CF2U,CF2V,'k','LineWidth',2);
plot(CF3U,CF3V,'k','LineWidth',2);
plot(CF4U,CF4V,'k','LineWidth',2);
plot(CF5U,CF5V,'k','LineWidth',2);
plot(CF6U,CF6V,'k','LineWidth',2);
plot(CF7U,CF7V,'k','LineWidth',2);
plot(CF8U,CF8V,'k','LineWidth',2);
plot(CF9U,CF9V,'k','LineWidth',2);
plot(CF10U,CF10V,'k','LineWidth',2);

axis off;
set(gcf, 'color', [1 1 1]);

hold off;
    
   F(i)=getframe(gcf);

    
    
end
    movie2avi(F,'RotationCombined.avi','quality',100);
    movie2avi(F,'RotationCombined_lowRes.avi','quality',75);
    
end