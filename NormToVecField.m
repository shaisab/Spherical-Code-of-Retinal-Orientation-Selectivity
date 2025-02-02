function [alpha,beta] = NormToVecField(normal)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This function computes the values of alpha and beta given the
%   components of a normal. Here we are assuming that the axis of rotation
%   of the semicircular canels is about the center of the eye.
%
%   Input Parameters:
%   1. normal -- the normal vector corresponding to the particular
%   semicircular canal. The vector should contain three components.
%
%   Ouput:
%   1. alpha -- aziumuthal coordinate of the center of the axis of
%   rotation.
%   2. beta -- altitude coordinate of the center of the axis of rotation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%opticaxis=[-0.712,0.480,-0.512];    % Nathan's data
opticaxis=[0.406,-0.833,0.375];      % optL

rot=vrrotvec2mat(vrrotvec(opticaxis',[0,0,1]'));
normal=rot*normal';

n=normal./norm(normal);
alpha=atan2(n(2),n(1));
beta=-acos(n(3));

