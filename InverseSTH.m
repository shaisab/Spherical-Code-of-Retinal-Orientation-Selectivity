function [t,phi] = InverseSTH(alpha,beta,R,s,theta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This function computes the inverse of the functional relationship
%   between (t,phi) and (s,theta). This calculation is crucial for
%   computing the vector fields in the flattened coordinates.
%
%   Inputs: 
%   1. alpha= azimuthal rotation angle.
%   2. beta= altitude rotation angle.
%   3. R= radius of sphere.
%   4. s= geodesic radius from south pole
%   5. theta= polar angle around south pole
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Rotation to determine t and phi

%x0= center point of axis of rotation
x0=[R*cos(alpha)*sin(beta),R*sin(beta)*sin(alpha),-R*cos(beta)];
%x1= point on sphere we want to convert to t and phi coordinates
x1=[R*sin(s/R)*cos(theta),R*sin(s/R)*sin(theta),-R*cos(s/R)];

%construction of rotation matrices
R1=[cos(alpha), -sin(alpha), 0; sin(alpha), cos(alpha), 0; 0, 0, 1];
R2=[cos(beta), 0, -sin(beta); 0, 1, 0; sin(beta), 0, cos(beta)];

%rotated values of y0 and y1.
y0=x0'\(R1*R2);
y1=x1'\(R1*R2);

%% Calculation of t and phi
if y1(3)<0,
    t=R*asin(sqrt((y1(1)-y0(1))^2+(y1(2)-y0(2))^2));
    phi=atan2(y1(2)-y0(2),y1(1)-y0(1));
else
    t=pi*R-R*asin(sqrt((y1(1)-y0(1))^2+(y1(2)-y0(2))^2));
    phi=atan2(y1(2)-y0(2),y1(1)-y0(1));
end

t=real(t);
