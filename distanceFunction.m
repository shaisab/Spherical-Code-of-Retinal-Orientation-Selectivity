function [cutoff]=distanceFunction(alpha,beta,x,y,z,R,loc,width)

% This a function which computes a weight (cutoff). 
% The code takes in alpha beta for a canal and x, y, z coordinates of the point you want to measure the geodesic distance from. 
% The parameter loc is the radial distance from the canal beyond which you want points to receive larger weights. 
% The parameter width controls the width of the transition region. I attached an image the illustrates this function.

%% Inputs
%
%   1. alpha coordinate of canal
%   2. beta coordinate of canal
%   3. x coordiante of point
%   4. y coordinate of point
%   5. z coordinate of point
%   6. R radius of sphere
%   7. loc location of the center of the transition layer
%   8. width of the transition layer

x0=R.*cos(alpha).*sin(beta);
y0=R.*sin(alpha).*sin(beta);
z0=R-R.*cos(beta);

% x0=repmat(x0,1,size(x,2));
% y0=repmat(y0,1,size(y,2));
% z0=repmat(z0,1,size(z,2));

d=sqrt((x-x0).^2+(y-y0).^2+(z-z0).^2);

D=real(2*asin(d./(2*R)));

cutoff=(tanh((D-loc)./width)+1)/2;

end

