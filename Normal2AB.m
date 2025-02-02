function [alpha,beta] = Normal2AB(Normal)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This codes takes a normal and converts it to an alpha beta in ANY
%   coordinate system.
%
%   Inputs
%   1. Normal - This is a normal measured from the origin outward.
%
%   Outputs
%   1. alpha - The alpha coordinate as measured counterclockwsie from the
%   x-axis.
%   2. beta - the azimuth angle as measured relative to the equator.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha=atan2(Normal(2),Normal(1));
beta=acos(-Normal(3));
end

