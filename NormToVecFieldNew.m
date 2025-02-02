function [alpha,beta] = NormToVecFieldNew(normal)
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
%
%
%   Ouput:
%   1. alpha -- aziumuthal coordinate of the center of the axis of
%   rotation.
%   2. beta -- altitude coordinate of the center of the axis of rotation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%opticaxis=[-0.712,0.480,-0.512];    % Nathan's data
opticaxis=[0.4636,-0.8030,0.3746];      % our made up data

rot=vrrotvec2mat(vrrotvec(opticaxis',[0,0,1]'));
normal=rot*normal';


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

%% Solution
alpha=real(atan2(Ny,Nx));
beta=real(acos(2*Nz^2));



