

figure('Color',[1 1 1]);
[x,y,z] = sphere(40);
c=repmat(0.5,size(x,1),size(x,2));
surf(x,y,z,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
axis equal
shading interp;
alpha(0.5)