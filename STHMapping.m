function []= STHMapping(typeStr)

% typestr can be either 'ONDS', 'ONOFFDS', 'OFFDS', or 'retroONDS'. 

%STHMapping

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This script compares the vector fields obtained from the elastic model
%   with that from Shai's data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Physical Parameters.
% This data is imported from data exported from the mapping script. These
% values correspond to measured values on the flattened retina.

%clear all;
%close all;

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

%% Sector Data
%We import the numerical solutions to the mapping problem from each sector.

%Sector1;
Sec1Data=dlmread('Sector1'); %radial stretching factor.
n=size(Sec1Data,1);
RHO1=Sec1Data(1:n,1:n);
F1=Sec1Data(1:n,n+1:2*n);
U1=RHO1.*cos(F1);
V1=RHO1.*sin(F1);
s1=linspace(0,M,n+1);
s1=linspace(s1(2),M,n);
th1=linspace(a1,a2,n);


ds=s1(2)-s1(1);
DS=DiffX(n)/ds;
RHOS1=(DS*RHO1')';
FS1=(DS*F1')';

dth=th1(2)-th1(1);
DTH=DiffX(n)/dth;
RHOTH1=(DTH*RHO1);
FTH1=(DTH*F1);

[S1,TH1]=meshgrid(s1,th1);

C=abs(R.^2*sin(S1./R).^2-S1.^2);
surf(RHO1.*cos(F1),RHO1.*sin(F1),zeros(n,n),C);
view(0,90);
colormap(jet);
caxis([0,1]);
colorbar;
colormap('bone');
shading interp;

%Sector2;
Sec2Data=dlmread('Sector2'); %radial stretching factor.
n=size(Sec2Data,1);
RHO2=Sec2Data(1:n,1:n);
F2=Sec2Data(1:n,n+1:2*n);
U2=RHO2.*cos(F2);
V2=RHO2.*sin(F2);
s2=linspace(0,M,n+1);
s2=linspace(s2(2),M,n);
th2=linspace(a2,a3,n);

ds=s2(2)-s2(1);
DS=DiffX(n)/ds;
RHOS2=(DS*RHO2')';
FS2=(DS*F2')';

dth=th2(2)-th2(1);
DTH=DiffX(n)/dth;
RHOTH2=(DTH*RHO2);
FTH2=(DTH*F2);

[S2,TH2]=meshgrid(s2,th2);

C=abs(R.^2*sin(S2./R).^2-S2.^2);
surf(RHO2.*cos(F2),RHO2.*sin(F2),zeros(n,n),C);
view(0,90);
colormap(jet);
caxis([0,1]);
colorbar;
colormap('bone');
shading interp;

%Sector3;
Sec3Data=dlmread('Sector3'); %radial stretching factor.
n=size(Sec3Data,1);
RHO3=Sec3Data(1:n,1:n);
F3=Sec3Data(1:n,n+1:2*n);
U3=RHO3.*cos(F3);
V3=RHO3.*sin(F3);
s3=linspace(0,M,n+1);
s3=linspace(s3(2),M,n);
th3=linspace(a3,a4,n);

[S3,TH3]=meshgrid(s3,th3);

ds=s3(2)-s3(1);
DS=DiffX(n)/ds;
RHOS3=(DS*RHO3')';
FS3=(DS*F3')';

dth=th3(2)-th3(1);
DTH=DiffX(n)/dth;
RHOTH3=(DTH*RHO3);
FTH3=(DTH*F3);

C=abs(R.^2*sin(S3./R).^2-S3.^2);
surf(RHO3.*cos(F3),RHO3.*sin(F3),zeros(n,n),C);
view(0,90);
colormap(jet);
caxis([0,1]);
colorbar;
colormap('bone');
shading interp;

%Sector4;
Sec4Data=dlmread('Sector4'); %radial stretching factor.
n=size(Sec4Data,1);
RHO4=Sec4Data(1:n,1:n);
F4=Sec4Data(1:n,n+1:2*n);
U4=RHO4.*cos(F4);
V4=RHO4.*sin(F4);
s4=linspace(0,M,n+1);
s4=linspace(s4(2),M,n);
th4=linspace(a4,a1+2*pi,n);

ds=s4(2)-s4(1);
DS=DiffX(n)/ds;
RHOS4=(DS*RHO4')';
FS4=(DS*F4')';

dth=th4(2)-th4(1);
DTH=DiffX(n)/dth;
RHOTH4=(DTH*RHO4);
FTH4=(DTH*F4);

[S4,TH4]=meshgrid(s4,th4);

C=abs(R.^2*sin(S4./R).^2-S4.^2);
surf(RHO4.*cos(F4),RHO4.*sin(F4),zeros(n,n),C);
view(0,90);
colormap(jet);
caxis([0,1]);
colorbar;
colormap('bone');
shading interp;

%%

[filedir,go]=bbdirselector('select data folder',cd);
if ~go disp('no folder selected'); return; end

namesMap=getFileList(filedir,'map_',0,'anywhere');
for i=1:size(namesMap,2)   
    if strfind(namesMap{i},'.mat')
    load(num2str(namesMap{i}));
    end
end

% load('map_May_27_2014_2.mat');
data=map.SzMapOutputDS;

U=data(:,1);
V=data(:,2);
vecU=data(:,3);
vecV=data(:,4);
Sec=data(:,5);


switch typeStr
    case 'ONDS' 
type=map.ONDS;
    case 'ONOFFDS' 
type=map.ONOFFDS;
    case 'OFFDS' 
type=map.OFFDS;
    case 'retroONDS' 
type=map.retroONDS;
end


n=size(U,1);

hold on;
j=0;
for i=1:n,
    
    if type(i)==1
        j=j+1;
        if Sec(i)==1,
            
            temp=(U(i)-U1).^2+(V(i)-V1).^2;
            [minval,ind]=min(temp(:));
            [a,b]=ind2sub([size(temp,1),size(temp,2)],ind);
            (U(i)-U1(a,b)).^2+(V(i)-V1(a,b)).^2
            
            [g1,g2] = InverseVec(vecU(i),vecV(i),RHO1(a,b),F1(a,b),RHOS1(a,b),FS1(a,b),RHOTH1(a,b),FTH1(a,b));
            vecUN(j)=vecU(i);
            vecVN(j)=vecV(i);
            UN(j)=U(i);
            VN(j)=V(i);
            Sn(j)=S1(a,b);
            THn(j)=TH1(a,b);
            vecSn(j)=g1;
            vecTHn(j)=g2;
            
            
            
        elseif Sec(i)==2,
            
            temp=(U(i)-U2).^2+(V(i)-V2).^2;
            [minval,ind]=min(temp(:));
            [a,b]=ind2sub([size(temp,1),size(temp,2)],ind);
            (U(i)-U2(a,b)).^2+(V(i)-V2(a,b)).^2
            
            [g1,g2] = InverseVec(vecU(i),vecV(i),RHO2(a,b),F2(a,b),RHOS2(a,b),FS2(a,b),RHOTH2(a,b),FTH2(a,b));
            vecUN(j)=vecU(i);
            vecVN(j)=vecV(i);
            UN(j)=U(i);
            VN(j)=V(i);
            Sn(j)=S2(a,b);
            THn(j)=TH2(a,b);
            vecSn(j)=g1;
            vecTHn(j)=g2;
            
            
        elseif Sec(i)==3,
            
            temp=(U(i)-U3).^2+(V(i)-V3).^2;
            [minval,ind]=min(temp(:));
            [a,b]=ind2sub([size(temp,1),size(temp,2)],ind);
            (U(i)-U3(a,b)).^2+(V(i)-V3(a,b)).^2
            
            [g1,g2] = InverseVec(vecU(i),vecV(i),RHO3(a,b),F3(a,b),RHOS3(a,b),FS3(a,b),RHOTH3(a,b),FTH3(a,b));
            vecUN(j)=vecU(i);
            vecVN(j)=vecV(i);
            UN(j)=U(i);
            VN(j)=V(i);
            Sn(j)=S3(a,b);
            THn(j)=TH3(a,b);
            vecSn(j)=g1;
            vecTHn(j)=g2;
        else
            temp=(U(i)-U4).^2+(V(i)-V4).^2;
            [minval,ind]=min(temp(:));
            [a,b]=ind2sub([size(temp,1),size(temp,2)],ind);
            (U(i)-U4(a,b)).^2+(V(i)-V4(a,b)).^2
            
            [g1,g2] = InverseVec(vecU(i),vecV(i),RHO4(a,b),F4(a,b),RHOS4(a,b),FS4(a,b),RHOTH4(a,b),FTH4(a,b));
            vecUN(j)=vecU(i);
            vecVN(j)=vecV(i);
            UN(j)=U(i);
            VN(j)=V(i);
            Sn(j)=S4(a,b);
            THn(j)=TH4(a,b);
            vecSn(j)=g1;
            vecTHn(j)=g2;
            
        end
    end
end

quiver(UN,VN,vecUN,vecVN,.5,'linewidth',2);
axis([-1,1,-1,1]);

hold off;

figure;
vecSntemp=vecSn./(vecSn.^2+vecTHn.^2).^(1/2);
vecTHntemp=vecTHn./(vecSn.^2+vecTHn.^2).^(1/2);

vecSn=vecSntemp;
vecTHn=vecTHntemp;

Sn=Sn*2*pi/M;

quiver(Sn,THn,vecSn,vecTHn,.3,'linewidth',2);
axis([0,2*pi,0,2*pi]);

Sn=Sn';
THn=THn';
THn=mod(THn,2*pi);
vecSn=vecSn';
vecTHn=vecTHn';

ind=strfind(filedir,'\');
filename=filedir;
filename(ind(end):length(filename))=[];
filename(1:ind(end-1))=[];
dlmwrite(['STHData_',typeStr,'_',filename],[Sn,THn,vecSn,vecTHn]);
hold off;

end

