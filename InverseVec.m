function [g1,g2] = InverseVec(h1,h2,rho,f,DSrho,DSf,DTHrho,DTHf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   This function computes the components of a vector in the (s,th)
%   coordinates given values in the (u,v) coordinates.
%
%   Input variables
%   1. h1 - component of vector in the u direction.
%   2. h2 - component of vector in the v direction.
%   3. rho - the radial coordinate of the mapping F corresponding to the
%   coordinate (u,v).
%   4. f - the angular coordinate of the mapping F corresponding to the
%   coordinate (u,v).
%   5. DSrho - derivative of rho with respect to s at (u,v).
%   6. DTHrho - derivative of rho with respect to theta at (u,v).
%   7. DSf - derivative of f with respect to s at (u,v).
%   8. DTHf - deriviatve of f with respect to th at (u,v).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Construction of Jacobian matrix
J=zeros(2,2);

det=DTHf*rho*DSrho-DSf*rho*DTHrho;

J(1,1)=1/det*(DTHf*rho*cos(f)+DTHrho*sin(f));
J(1,2)=1/det*(-1*DSf*rho*cos(f)-DSrho*sin(f));
J(2,1)=1/det*(-1*DTHrho*cos(f)+DTHf*rho*sin(f));
J(2,2)=1/det*(DSrho*cos(f)-DSf*rho*sin(f));

g1=J(1,1)*h1+J(2,1)*h2;
g2=J(1,2)*h1+J(2,2)*h2;

temp1=g1/(sqrt(g1^2+g2^2));
temp2=g2/(sqrt(g1^2+g2^2));

g1=temp1;
g2=temp2;

end

