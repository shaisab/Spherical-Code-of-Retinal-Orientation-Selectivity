
alpha=[10,95,185,275];
beta=[120,105,60,70];

for i=1:4
x(i)=cos(alpha(i)).*sin(beta(i));
y(i)=sin(alpha(i)).*sin(beta(i));
z(i)=-cos(beta(i));
end

figure;
quiver3(0,0,0,0,0,1,'k','LineWidth',1.5);
scatter3(x,y,z);
xlabel('x');ylabel('y');zlabel('z');
xlim([0-1 1]);ylim([0-1 1]);zlim([0-1 1]);