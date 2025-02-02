function [alpha,beta] = GreatCircle(normal,n)

[a,b,~]=ABGlobalVisual2RetinaRight(normal);

a=deg2rad(a);
b=deg2rad(b);

N1=cos(a)*sin(b);
N2=sin(a)*sin(b);
N3=-cos(b);

x=linspace(-1,1,n);
discriment=N2^2+N3^2-x.^2;
k=find(discriment>0);

y1=-(N1*N2*x-N3*sqrt(discriment))/(N2^2+N3^2);
z1=-(N1*N3*x+N2*sqrt(discriment))/(N2^2+N3^2);
y2=-(N1*N2*x+N3*sqrt(discriment))/(N2^2+N3^2);
z2=-(N1*N3*x-N2*sqrt(discriment))/(N2^2+N3^2);

% figure;
% plot3(x(k),y1(k),z1(k),x(k),y2(k),z2(k));

x=[x(k),fliplr(x(k))];
y=[y1(k),fliplr(y2(k))];
z=[z1(k),fliplr(z2(k))];

x=[x,x(1)];
y=[y,y(1)];
z=[z,z(1)];

% figure;
% plot3(x,y,z);

beta=real(acos(-1*z));
alpha=mod(real(atan2(y,x)),2*pi);

% figure;
% plot(alpha,beta,'.');
end

