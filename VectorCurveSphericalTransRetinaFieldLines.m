function [] = VectorCurveSphericalTransRetinaFieldLines(alphabar,beta,discPoints,VizAngle,lc)

VizAngle=deg2rad(VizAngle);
alphabar=deg2rad(alphabar);
beta=deg2rad(beta);

RV=2;
%Computing alphabar, beta given the normal and the optical axis.
%[alphabar,beta]=AB(Normal);

% translation
% temp=linspace(0,pi,10);
% CL=temp(2)-temp(1);
% 
% n=1000; %Number of discretization points for the curve.
% t=linspace(0,pi,n); %parametrization coordinates for the curve
% dt=t(2)-t(1);

% rotation
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
C1=[RV*cos(0*CL)*sin(t);RV*sin(0*CL)*sin(t);-RV*cos(t).*ones(1,n)];
C2=[RV*cos(1*CL)*sin(t);RV*sin(1*CL)*sin(t);-RV*cos(t).*ones(1,n)];
C3=[RV*cos(2*CL)*sin(t);RV*sin(2*CL)*sin(t);-RV*cos(t).*ones(1,n)];
C4=[RV*cos(3*CL)*sin(t);RV*sin(3*CL)*sin(t);-RV*cos(t).*ones(1,n)];
C5=[RV*cos(4*CL)*sin(t);RV*sin(4*CL)*sin(t);-RV*cos(t).*ones(1,n)];
C6=[RV*cos(5*CL)*sin(t);RV*sin(5*CL)*sin(t);-RV*cos(t).*ones(1,n)];
C7=[RV*cos(6*CL)*sin(t);RV*sin(6*CL)*sin(t);-RV*cos(t).*ones(1,n)];
C8=[RV*cos(7*CL)*sin(t);RV*sin(7*CL)*sin(t);-RV*cos(t).*ones(1,n)];
C9=[RV*cos(8*CL)*sin(t);RV*sin(8*CL)*sin(t);-RV*cos(t).*ones(1,n)];
C10=[RV*cos(9*CL)*sin(t);RV*sin(9*CL)*sin(t);-RV*cos(t).*ones(1,n)];
C11=[RV*cos(10*CL)*sin(t);RV*sin(10*CL)*sin(t);-RV*cos(t).*ones(1,n)];
C12=[RV*cos(11*CL)*sin(t);RV*sin(11*CL)*sin(t);-RV*cos(t).*ones(1,n)];
C13=[RV*cos(12*CL)*sin(t);RV*sin(12*CL)*sin(t);-RV*cos(t).*ones(1,n)];
C14=[RV*cos(13*CL)*sin(t);RV*sin(13*CL)*sin(t);-RV*cos(t).*ones(1,n)];
C15=[RV*cos(14*CL)*sin(t);RV*sin(14*CL)*sin(t);-RV*cos(t).*ones(1,n)];
C16=[RV*cos(15*CL)*sin(t);RV*sin(15*CL)*sin(t);-RV*cos(t).*ones(1,n)];
C17=[RV*cos(16*CL)*sin(t);RV*sin(16*CL)*sin(t);-RV*cos(t).*ones(1,n)];
C18=[RV*cos(17*CL)*sin(t);RV*sin(17*CL)*sin(t);-RV*cos(t).*ones(1,n)];
C19=[RV*cos(18*CL)*sin(t);RV*sin(18*CL)*sin(t);-RV*cos(t).*ones(1,n)];
C20=[RV*cos(19*CL)*sin(t);RV*sin(19*CL)*sin(t);-RV*cos(t).*ones(1,n)];

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
C11=Rot2*Rot1*C11;
C12=Rot2*Rot1*C12;
C13=Rot2*Rot1*C13;
C14=Rot2*Rot1*C14;
C15=Rot2*Rot1*C15;
C16=Rot2*Rot1*C16;
C17=Rot2*Rot1*C17;
C18=Rot2*Rot1*C18;
C19=Rot2*Rot1*C19;
C20=Rot2*Rot1*C20;

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
sv11=RV*acos(-C11(3,:)/RV);
sv12=RV*acos(-C12(3,:)/RV);
sv13=RV*acos(-C13(3,:)/RV);
sv14=RV*acos(-C14(3,:)/RV);
sv15=RV*acos(-C15(3,:)/RV);
sv16=RV*acos(-C16(3,:)/RV);
sv17=RV*acos(-C17(3,:)/RV);
sv18=RV*acos(-C18(3,:)/RV);
sv19=RV*acos(-C19(3,:)/RV);
sv20=RV*acos(-C20(3,:)/RV);

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
InRange11=sv11>=RV*(pi-VizAngle/2);
InRange12=sv12>=RV*(pi-VizAngle/2);
InRange13=sv13>=RV*(pi-VizAngle/2);
InRange14=sv14>=RV*(pi-VizAngle/2);
InRange15=sv15>=RV*(pi-VizAngle/2);
InRange16=sv16>=RV*(pi-VizAngle/2);
InRange17=sv17>=RV*(pi-VizAngle/2);
InRange18=sv18>=RV*(pi-VizAngle/2);
InRange19=sv19>=RV*(pi-VizAngle/2);
InRange20=sv20>=RV*(pi-VizAngle/2);

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
C11(1,~InRange11)=NaN;
C12(1,~InRange12)=NaN;
C13(1,~InRange13)=NaN;
C14(1,~InRange14)=NaN;
C15(1,~InRange15)=NaN;
C16(1,~InRange16)=NaN;
C17(1,~InRange17)=NaN;
C18(1,~InRange18)=NaN;
C19(1,~InRange19)=NaN;
C20(1,~InRange20)=NaN;

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
C4X=(Dt*C4(1,:)')';
C4Y=(Dt*C4(2,:)')';
C4Z=(Dt*C4(3,:)')';
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
C11X=(Dt*C10(1,:)')';
C11Y=(Dt*C10(2,:)')';
C11Z=(Dt*C10(3,:)')';
C12X=(Dt*C12(1,:)')';
C12Y=(Dt*C12(2,:)')';
C12Z=(Dt*C12(3,:)')';
C13X=(Dt*C13(1,:)')';
C13Y=(Dt*C13(2,:)')';
C13Z=(Dt*C13(3,:)')';
C14X=(Dt*C14(1,:)')';
C14Y=(Dt*C14(2,:)')';
C14Z=(Dt*C14(3,:)')';
C15X=(Dt*C15(1,:)')';
C15Y=(Dt*C15(2,:)')';
C15Z=(Dt*C15(3,:)')';
C16X=(Dt*C16(1,:)')';
C16Y=(Dt*C16(2,:)')';
C16Z=(Dt*C16(3,:)')';
C17X=(Dt*C17(1,:)')';
C17Y=(Dt*C17(2,:)')';
C17Z=(Dt*C17(3,:)')';
C18X=(Dt*C18(1,:)')';
C18Y=(Dt*C18(2,:)')';
C18Z=(Dt*C18(3,:)')';
C19X=(Dt*C19(1,:)')';
C19Y=(Dt*C19(2,:)')';
C19Z=(Dt*C19(3,:)')';
C20X=(Dt*C20(1,:)')';
C20Y=(Dt*C20(2,:)')';
C20Z=(Dt*C20(3,:)')';

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
normC11=(C11X.^2+C11Y.^2+C11Z.^2).^(1/2)/scale;
normC12=(C12X.^2+C12Y.^2+C12Z.^2).^(1/2)/scale;
normC13=(C13X.^2+C13Y.^2+C13Z.^2).^(1/2)/scale;
normC14=(C14X.^2+C14Y.^2+C14Z.^2).^(1/2)/scale;
normC15=(C15X.^2+C15Y.^2+C15Z.^2).^(1/2)/scale;
normC16=(C16X.^2+C16Y.^2+C16Z.^2).^(1/2)/scale;
normC17=(C17X.^2+C17Y.^2+C17Z.^2).^(1/2)/scale;
normC18=(C18X.^2+C18Y.^2+C18Z.^2).^(1/2)/scale;
normC19=(C19X.^2+C19Y.^2+C19Z.^2).^(1/2)/scale;
normC20=(C20X.^2+C20Y.^2+C20Z.^2).^(1/2)/scale;

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
C11X=C11X./normC11;
C11Y=C11Y./normC11;
C11Z=C11Z./normC11;
C12X=C12X./normC12;
C12Y=C12Y./normC12;
C12Z=C12Z./normC12;
C13X=C13X./normC13;
C13Y=C13Y./normC13;
C13Z=C13Z./normC13;
C14X=C14X./normC14;
C14Y=C14Y./normC14;
C14Z=C14Z./normC14;
C15X=C15X./normC15;
C15Y=C15Y./normC15;
C15Z=C15Z./normC15;
C16X=C16X./normC16;
C16Y=C16Y./normC16;
C16Z=C16Z./normC16;
C17X=C17X./normC17;
C17Y=C17Y./normC17;
C17Z=C17Z./normC17;
C18X=C18X./normC18;
C18Y=C18Y./normC18;
C18Z=C18Z./normC18;
C19X=C19X./normC19;
C19Y=C19Y./normC19;
C19Z=C19Z./normC19;
C20X=C20X./normC20;
C20Y=C20Y./normC20;
C20Z=C20Z./normC20;

%% Reading in data from retina
%This part of the code computes all of the data needed for the flattening
%map.

phys=dlmread('phys');

R=phys(1)*1.1;
M=phys(2)*1.1; %This should be normalized to 1.

m1=phys(3);
m2=phys(4);
m3=phys(5);
m4=phys(6);

nu=phys(7);

a1=phys(8);
a2=phys(9);
%a3=phys(10);
a3=pi;
%a4=phys(11);
a4=3*pi/2;

ntemp=discPoints;
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
sret11=-2*M/VizAngle*acos(-C11(3,:)/RV)+2*M*pi/VizAngle;
thret11=atan2(C11(2,:),C11(1,:))+pi;
sret12=-2*M/VizAngle*acos(-C12(3,:)/RV)+2*M*pi/VizAngle;
thret12=atan2(C12(2,:),C12(1,:))+pi;
sret13=-2*M/VizAngle*acos(-C13(3,:)/RV)+2*M*pi/VizAngle;
thret13=atan2(C13(2,:),C13(1,:))+pi;
sret14=-2*M/VizAngle*acos(-C14(3,:)/RV)+2*M*pi/VizAngle;
thret14=atan2(C14(2,:),C14(1,:))+pi;
sret15=-2*M/VizAngle*acos(-C15(3,:)/RV)+2*M*pi/VizAngle;
thret15=atan2(C15(2,:),C15(1,:))+pi;
sret16=-2*M/VizAngle*acos(-C16(3,:)/RV)+2*M*pi/VizAngle;
thret16=atan2(C16(2,:),C16(1,:))+pi;
sret17=-2*M/VizAngle*acos(-C17(3,:)/RV)+2*M*pi/VizAngle;
thret17=atan2(C17(2,:),C17(1,:))+pi;
sret18=-2*M/VizAngle*acos(-C18(3,:)/RV)+2*M*pi/VizAngle;
thret18=atan2(C18(2,:),C18(1,:))+pi;
sret19=-2*M/VizAngle*acos(-C19(3,:)/RV)+2*M*pi/VizAngle;
thret19=atan2(C19(2,:),C19(1,:))+pi;
sret20=-2*M/VizAngle*acos(-C20(3,:)/RV)+2*M*pi/VizAngle;
thret20=atan2(C20(2,:),C20(1,:))+pi;

% 3-dimensional coordinate of curves on the retina.
C1ret(1,:)=R*cos(thret1).*sin(sret1/R);
C1ret(2,:)=R*sin(thret1).*sin(sret1/R);
C1ret(3,:)=R-R*cos(sret1/R);
C2ret(1,:)=R*cos(thret2).*sin(sret2/R);
C2ret(2,:)=R*sin(thret2).*sin(sret2/R);
C2ret(3,:)=R-R*cos(sret2/R);
C3ret(1,:)=R*cos(thret3).*sin(sret3/R);
C3ret(2,:)=R*sin(thret3).*sin(sret3/R);
C3ret(3,:)=R-R*cos(sret3/R);
C4ret(1,:)=R*cos(thret4).*sin(sret4/R);
C4ret(2,:)=R*sin(thret4).*sin(sret4/R);
C4ret(3,:)=R-R*cos(sret4/R);
C5ret(1,:)=R*cos(thret5).*sin(sret5/R);
C5ret(2,:)=R*sin(thret5).*sin(sret5/R);
C5ret(3,:)=R-R*cos(sret5/R);
C6ret(1,:)=R*cos(thret6).*sin(sret6/R);
C6ret(2,:)=R*sin(thret6).*sin(sret6/R);
C6ret(3,:)=R-R*cos(sret6/R);
C7ret(1,:)=R*cos(thret7).*sin(sret7/R);
C7ret(2,:)=R*sin(thret7).*sin(sret7/R);
C7ret(3,:)=R-R*cos(sret7/R);
C8ret(1,:)=R*cos(thret8).*sin(sret8/R);
C8ret(2,:)=R*sin(thret8).*sin(sret8/R);
C8ret(3,:)=R-R*cos(sret8/R);
C9ret(1,:)=R*cos(thret9).*sin(sret9/R);
C9ret(2,:)=R*sin(thret9).*sin(sret9/R);
C9ret(3,:)=R-R*cos(sret9/R);
C10ret(1,:)=R*cos(thret10).*sin(sret10/R);
C10ret(2,:)=R*sin(thret10).*sin(sret10/R);
C10ret(3,:)=R-R*cos(sret10/R);
C11ret(1,:)=R*cos(thret11).*sin(sret11/R);
C11ret(2,:)=R*sin(thret11).*sin(sret11/R);
C11ret(3,:)=R-R*cos(sret11/R);
C12ret(1,:)=R*cos(thret12).*sin(sret12/R);
C12ret(2,:)=R*sin(thret12).*sin(sret12/R);
C12ret(3,:)=R-R*cos(sret12/R);
C13ret(1,:)=R*cos(thret13).*sin(sret13/R);
C13ret(2,:)=R*sin(thret13).*sin(sret13/R);
C13ret(3,:)=R-R*cos(sret13/R);
C14ret(1,:)=R*cos(thret14).*sin(sret14/R);
C14ret(2,:)=R*sin(thret14).*sin(sret14/R);
C14ret(3,:)=R-R*cos(sret14/R);
C15ret(1,:)=R*cos(thret15).*sin(sret15/R);
C15ret(2,:)=R*sin(thret15).*sin(sret15/R);
C15ret(3,:)=R-R*cos(sret15/R);
C16ret(1,:)=R*cos(thret16).*sin(sret16/R);
C16ret(2,:)=R*sin(thret16).*sin(sret16/R);
C16ret(3,:)=R-R*cos(sret16/R);
C17ret(1,:)=R*cos(thret17).*sin(sret17/R);
C17ret(2,:)=R*sin(thret17).*sin(sret17/R);
C17ret(3,:)=R-R*cos(sret17/R);
C18ret(1,:)=R*cos(thret18).*sin(sret18/R);
C18ret(2,:)=R*sin(thret18).*sin(sret18/R);
C18ret(3,:)=R-R*cos(sret18/R);
C19ret(1,:)=R*cos(thret19).*sin(sret19/R);
C19ret(2,:)=R*sin(thret19).*sin(sret19/R);
C19ret(3,:)=R-R*cos(sret19/R);
C20ret(1,:)=R*cos(thret20).*sin(sret20/R);
C20ret(2,:)=R*sin(thret20).*sin(sret20/R);
C20ret(3,:)=R-R*cos(sret20/R);


%% Plotting curves on retina

hold on;
plot3(C1ret(1,:),C1ret(2,:),C1ret(3,:),lc,'LineWidth',1);
plot3(C2ret(1,:),C2ret(2,:),C2ret(3,:),lc,'LineWidth',1);
plot3(C3ret(1,:),C3ret(2,:),C3ret(3,:),lc,'LineWidth',1);
plot3(C4ret(1,:),C4ret(2,:),C4ret(3,:),lc,'LineWidth',1);
plot3(C5ret(1,:),C5ret(2,:),C5ret(3,:),lc,'LineWidth',1);
plot3(C6ret(1,:),C6ret(2,:),C6ret(3,:),lc,'LineWidth',1);
plot3(C7ret(1,:),C7ret(2,:),C7ret(3,:),lc,'LineWidth',1);
plot3(C8ret(1,:),C8ret(2,:),C8ret(3,:),lc,'LineWidth',1);
plot3(C9ret(1,:),C9ret(2,:),C9ret(3,:),lc,'LineWidth',1);
plot3(C10ret(1,:),C10ret(2,:),C10ret(3,:),lc,'LineWidth',1);
plot3(C11ret(1,:),C11ret(2,:),C11ret(3,:),lc,'LineWidth',1);
plot3(C12ret(1,:),C12ret(2,:),C12ret(3,:),lc,'LineWidth',1);
plot3(C13ret(1,:),C13ret(2,:),C13ret(3,:),lc,'LineWidth',1);
plot3(C14ret(1,:),C14ret(2,:),C14ret(3,:),lc,'LineWidth',1);
plot3(C15ret(1,:),C15ret(2,:),C15ret(3,:),lc,'LineWidth',1);
plot3(C16ret(1,:),C16ret(2,:),C16ret(3,:),lc,'LineWidth',1);
plot3(C17ret(1,:),C17ret(2,:),C17ret(3,:),lc,'LineWidth',1);
plot3(C18ret(1,:),C18ret(2,:),C18ret(3,:),lc,'LineWidth',1);
plot3(C19ret(1,:),C19ret(2,:),C19ret(3,:),lc,'LineWidth',1);
plot3(C20ret(1,:),C20ret(2,:),C20ret(3,:),lc,'LineWidth',1); 

hold on;


end
