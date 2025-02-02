function [alpha1,beta1,alpha2,beta2] = NormToVecField1(normal)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This function computes the values of alpha and beta given the
%   components of a normal and the distance from the south pole of the
%   retina to the focal point of the lens. Here we are assuming that the 
%   axis of rotation is about a point on the retina that is on a line
%   parallel to the normal of the SCC and passes through the focal point of
%   the lens. 
%
%   Input Parameters:
%   1. normal -- the normal vector corresponding to the particular
%   semicircular canal. The vector should contain three components.
%   
%   2. d -- the distance from the south pole of the retina to the focal
%   point of the lens. 
%
%   3. R -- radius of the retina.
%
%
%   Ouput:
%   1. alpha -- aziumuthal coordinate of the center of the axis of
%   rotation.
%   2. beta -- altitude coordinate of the center of the axis of rotation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

normal

[filedir,go]=bbdirselector('select data folder',cd);
if ~go disp('no folder selected'); return; end

namesMR=getFileList(filedir,'meridian_radius',0,'anywhere');
if ~isempty(namesMR)
    MR=load(num2str(namesMR{1}));
end
M=MR(1);


d=1715;    %According to the 'eye anatomy' paper (Remtulla and Hallett, 1985), this posterior nodal distance is 1715 um.
d=d/M;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opticaxis=[-0.713,0.480,-0.512];       %(lens-optic nerve orgin) 
rot=vrrotvec2mat(vrrotvec(opticaxis',[0,0,1]'));
normal=rot*normal';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=normal./norm(normal);
Nx=n(1);
Ny=n(2);
Nz=n(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   The line formed by the normal is given by s*normal+d. The equation for
%   the values of s where this lines passes through the retina is given by
%   s^2*normal^2+d^2=1. This gives a quadratic equation the can be solved.
%   This will give two roots for s and hence we will have two outputs for
%   the values of alpha and beta. Once s is determined the corresponding
%   values of x,y,z can be calculated which allow alpha and beta to be
%   determined.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% First Solution
s=(-Nz*d+sqrt(Nz^2*d^2+1+d^2))/(1+d)^2;

x=s*Nx;
y=s*Ny;
z=s*Nz+d;


alpha1=real(atan(y/x))
beta1=real(acos(-z));


%% Second Solution

s=(-Nz*d-sqrt(Nz^2*d^2+1+d^2))/(1+d)^2;

x=s*Nx;
y=s*Ny;
z=s*Nz+d;

alpha2=real(atan(y/x))
beta2=real(acos(-z));
end

