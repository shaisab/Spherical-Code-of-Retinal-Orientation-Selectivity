function []= standardMapping(n)

% n = number of discretization points



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This script computes the optimal mapping for a "standard" retina. The
%   model of the flattened retina is based on minimizing the elastic energy
%   of a flattened sphere with four cuts made along merideans. The retina
%   of course is not a perfect sphere but is close eneough that this should
%   have little impact on the modeling. We assume that the retina is a
%   hyperelastic material with a linear stress strain relationship. 
%
%   Input Parameters:
%   1. R - Radius of the hemisphere.
%   2. M - Length of meridean connecting the south pole of the retina to
%   the boundary of S. We normalize all length scales so that this value is
%   1.
%   3. mi - Length of meridan connecting the south pole to cut i. A
%   critical assumption in the modeling is that mi/R<<1.
%   4. nu - Poisson ratio of the material.
%   5. ai - angular position of the cuts.
%   6. n - number of points used in discretization.
%
%   Independent Variables:
%   1. s - The arclength coordinate along meridians.
%   2. th - The azimuthal angle along the sphere.
%
%   Dependent Variables:
%   We assume that the flattening map is of the form:
%   F(s,th)=(rho(s,th)*cos(f(s,th)),rho(s,th)*sin(f(s,th))),
%   where the dependent variables - i.e. functions - are given by
%   1. rho - radial coordinate in the flattened coordinates.
%   2. f - polar angle in the flattened coordinates.
%
%   Ouput:
%   1. rho - this is the radial coordinate of the mapping. We output this
%   as a data file for each sector of the retina.
%   2. f - the is the angular coordiante of the mapping. We output this as
%   a data file for each sector of the retina as well.
%
%   The code determines the functions rho(s,th) and f(s,th) given the
%   values of the physical parameters by minimizing the elastic energy
%   using Matlab's built in lsqnonlin algorithm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all;

f1=figure;

set(f1,'position',[5 150 950 800]);
set(0, 'CurrentFigure', f1) %make f1 the active  figure

phys=dlmread('phys');

R=phys(1);
M=phys(2); %This should be normalized to 1.

m1=1/(2.5*M);
m2=1/(2.5*M);
m3=1/(2.5*M);
m4=1/(2.5*M);

nu=phys(7);

a1=0;
a2=pi/2;
a3=pi;
a4=3*pi/2;

%We discretize s so that the coordinate singularity at s=0 is removed.
s=linspace(0,M,n+1);
s=linspace(s(2),M,n);

%% Sector 1
%This part of the code minimizes the energy in the first sector. 

%Radial coordinates of cuts. This will be used to define boundary
%conditions in minimization routine.
[~,a]=min(abs(s-m1));
[~,b]=min(abs(s-m2));

th=linspace(a1,a2,n);

%Construction of initial guess. We assume perfect rotational symmetry for
%the initial guess. This matches the lowest order approximation to the
%solution near the center of the retina.
[S,TH]=meshgrid(s,th);



%Construction of derivative operators.
ds=s(2)-s(1);
dth=th(2)-th(1);
DS=DiffX(n)/ds;
DTH=DiffX(n)/dth;


guess(:,1:n)=S;
guess(:,n+1:2*n)=TH;
options=optimset('MaxFunEvals',1e+9,'MaxIter',1000,'Display','iter');
 

output = lsqnonlin(@(guess) SectorEnergy(guess,DS,DTH,nu,S,s,R,n,ds,dth,a,b,a1,a2),guess,[],[],options);

rho=output(1:n,1:n);
f=output(1:n,n+1:2*n);

%Construction of strain energy density for coloring the mapping.
DTrho=DTH*rho;
DSrho=(DS*rho')';
DTf=DTH*f;
DSf=(DS*f')';

gamma11=DSrho.^2+rho.^2.*DSf.^2-1;
gamma12=(DSrho.*DTrho+rho.^2.*DTf.*DSf)./(R*sin(S/R));
gamma22=(DTrho.^2+rho.^2.*DTf.^2-R^2.*sin(S/R).^2)./(R*sin(S/R)).^2;

C=nu*(gamma11+gamma22).^2+(1-nu)*(gamma11.^2+2*gamma12.^2+gamma22.^2);
C=C.*(R*sin(S/R));

%Plot of flatened retina
set(0, 'CurrentFigure', f1) %make f1 the active  figure
surf(rho.*cos(f),rho.*sin(f),zeros(n,n),C);
view(0,90);
colormap(jet);
caxis([0,.05]);
colorbar;
shading interp;
hold on;

%Outputting data to ASCII file. We record both the radial displacement rho
%and the angular displacment f. We also record the values of S and TH so
%that they do not have to be reconstructed.
dlmwrite(['StandardSector1_',num2str(n)],[rho,f,S,TH]);

%% Sector 2
%This part of the code minimizes the energy in the second sector. 

%Radial coordinates of cuts
[~,a]=min(abs(s-m2));
[~,b]=min(abs(s-m3));

th=linspace(a2,a3,n);

%Construction of initial guess
[S,TH]=meshgrid(s,th);

%Construction of derivative operators.
ds=s(2)-s(1);
dth=th(2)-th(1);
DS=DiffX(n)/ds;
DTH=DiffX(n)/dth;


% Energy minimization routine. We take as an initial guess for the
% algorithm that rhoguess=s. The initial guess that we take for a minimizer 
% will simply be the lowest order term in the asympotic expansion assuming 
% radial symmetry and matching the boundary conditions.

guess(:,1:n)=S;
guess(:,n+1:2*n)=TH;
options=optimset('MaxFunEvals',1e+9,'MaxIter',1000,'Display','iter');
 

output = lsqnonlin(@(guess) SectorEnergy(guess,DS,DTH,nu,S,s,R,n,ds,dth,a,b,a2,a3),guess,[],[],options);

rho=output(1:n,1:n);
f=output(1:n,n+1:2*n);

%Construction of strain energy density for coloring the mapping.
DTrho=DTH*rho;
DSrho=(DS*rho')';
DTf=DTH*f;
DSf=(DS*f')';

gamma11=DSrho.^2+rho.^2.*DSf.^2-1;
gamma12=(DSrho.*DTrho+rho.^2.*DTf.*DSf)./(R*sin(S/R));
gamma22=(DTrho.^2+rho.^2.*DTf.^2-R^2.*sin(S/R).^2)./(R*sin(S/R)).^2;

C=nu*(gamma11+gamma22).^2+(1-nu)*(gamma11.^2+2*gamma12.^2+gamma22.^2);
C=C.*(R*sin(S/R));

%Plot of flatened retina
surf(rho.*cos(f),rho.*sin(f),zeros(n,n),C);
view(0,90);
colormap(jet);
caxis([0,.05]);
colorbar;
shading interp;

%Outputting data to ASCII file. We record both the radial displacement rho
%and the angular displacment f. We also record the values of S and TH so
%that they do not have to be reconstructed.
dlmwrite(['StandardSector2_',num2str(n)],[rho,f,S,TH]);

%% Sector 3
%This part of the code minimizes the energy in the third sector. 

%Radial coordinates of cuts
[~,a]=min(abs(s-m3));
[~,b]=min(abs(s-m4));

th=linspace(a3,a4,n);

%Construction of initial guess
[S,TH]=meshgrid(s,th);

%Construction of derivative operators.
ds=s(2)-s(1);
dth=th(2)-th(1);
DS=DiffX(n)/ds;
DTH=DiffX(n)/dth;


% Energy minimization routine. We take as an initial guess for the
% algorithm that rhoguess=s. The initial guess that we take for a minimizer 
% will simply be the lowest order term in the asympotic expansion assuming 
% radial symmetry and matching the boundary conditions.

guess(:,1:n)=S;
guess(:,n+1:2*n)=TH;
options=optimset('MaxFunEvals',1e+9,'MaxIter',1000,'Display','iter');
 

output = lsqnonlin(@(guess) SectorEnergy(guess,DS,DTH,nu,S,s,R,n,ds,dth,a,b,a3,a4),guess,[],[],options);

rho=output(1:n,1:n);
f=output(1:n,n+1:2*n);

%Construction of strain energy density for coloring the mapping.
DTrho=DTH*rho;
DSrho=(DS*rho')';
DTf=DTH*f;
DSf=(DS*f')';

gamma11=DSrho.^2+rho.^2.*DSf.^2-1;
gamma12=(DSrho.*DTrho+rho.^2.*DTf.*DSf)./(R*sin(S/R));
gamma22=(DTrho.^2+rho.^2.*DTf.^2-R^2.*sin(S/R).^2)./(R*sin(S/R)).^2;

C=nu*(gamma11+gamma22).^2+(1-nu)*(gamma11.^2+2*gamma12.^2+gamma22.^2);
C=C.*(R*sin(S/R));

%Plot of flatened retina
surf(rho.*cos(f),rho.*sin(f),zeros(n,n),C);
view(0,90);
colormap(jet);
caxis([0,.05]);
colorbar;
shading interp;

%Outputting data to ASCII file. We record both the radial displacement rho
%and the angular displacment f. We also record the values of S and TH so
%that they do not have to be reconstructed.
dlmwrite(['StandardSector3_',num2str(n)],[rho,f,S,TH]);

%% Sector 4
%This part of the code minimizes the energy in the first sector. 

%Radial coordinates of cuts
[~,a]=min(abs(s-m4));
[~,b]=min(abs(s-m1));

th=linspace(a4,a1+2*pi,n);

%Construction of initial guess
[S,TH]=meshgrid(s,th);

%Construction of derivative operators.
ds=s(2)-s(1);
dth=th(2)-th(1);
DS=DiffX(n)/ds;
DTH=DiffX(n)/dth;


% Energy minimization routine. We take as an initial guess for the
% algorithm that rhoguess=s. The initial guess that we take for a minimizer 
% will simply be the lowest order term in the asympotic expansion assuming 
% radial symmetry and matching the boundary conditions.

guess(:,1:n)=S;
guess(:,n+1:2*n)=TH;
options=optimset('MaxFunEvals',1e+9,'MaxIter',1000,'Display','iter');
 

output = lsqnonlin(@(guess) SectorEnergy(guess,DS,DTH,nu,S,s,R,n,ds,dth,a,b,a4,a1+2*pi),guess,[],[],options);

rho=output(1:n,1:n);
f=output(1:n,n+1:2*n);

%Construction of strain energy density for coloring the mapping.
DTrho=DTH*rho;
DSrho=(DS*rho')';
DTf=DTH*f;
DSf=(DS*f')';

gamma11=DSrho.^2+rho.^2.*DSf.^2-1;
gamma12=(DSrho.*DTrho+rho.^2.*DTf.*DSf)./(R*sin(S/R));
gamma22=(DTrho.^2+rho.^2.*DTf.^2-R^2.*sin(S/R).^2)./(R*sin(S/R)).^2;

C=nu*(gamma11+gamma22).^2+(1-nu)*(gamma11.^2+2*gamma12.^2+gamma22.^2);
C=C.*(R*sin(S/R));

%Plot of flatened retina
surf(rho.*cos(f),rho.*sin(f),zeros(n,n),C);
view(0,90);
colormap(jet);
caxis([0,.05]);
colorbar;
shading interp;

%Outputting data to ASCII file. We record both the radial displacement rho
%and the angular displacment f. We also record the values of S and TH so
%that they do not have to be reconstructed.
dlmwrite(['StandardSector4_',num2str(n)],[rho,f,S,TH]);

saveas(f1,['standardMapping_',num2str(n)],'fig');
