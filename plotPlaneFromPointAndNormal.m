function plotPlaneFromPointAndNormal(point,normal)

% point = [1,2,3];
% normal = [1,1,2];

%# a plane is a*x+b*y+c*z+d=0
%# [a,b,c] is the normal. Thus, we have to calculate
%# d and we're set
d = -point*normal'; %'# dot product for less typing

%# create x,y
[xx,yy]=ndgrid(-1:1,-1:1);

%# calculate corresponding z
z = (-normal(1)*xx - normal(2)*yy - d)/normal(3);
if mean(mean(abs(z)))>1
    z=z./10;
    xx=xx./10;
    yy=yy./10;
end
%# plot the surface
%figure
hold on;
surf(xx,yy,z)
colormap gray
shading interp;
alpha(0.5)
%pause(1)


end