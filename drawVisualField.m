function []=drawVisualField()

%% draw shaded hemisphere
phys=dlmread('phys');
R=phys(1);
a=100;
b=100;
Alpha=linspace(0,2*pi,b+1);
Alpha=Alpha(2:b+1);
beta=linspace(pi/2,pi,a);
[A,B]=meshgrid(Alpha,beta);close 

%% transform to global space
f=1.03;
% figure;
% hold on
for i=1:100
    for j=1:100
[alphaRG(i,j),betaRG(i,j),Normalr]=rotateNormalToExtraPersonalSpace4('Right',rad2deg(A(i,j)),rad2deg(B(i,j)));
%plot3(Normalr(1),Normalr(2),Normalr(3),'o')
    end
end
%xlabel('X');ylabel('Y');zlabel('Z');


X=f*R*cosd(alphaRG).*sind(betaRG);
Y=f*R*sind(alphaRG).*sind(betaRG);
Z=-f*R*cosd(betaRG);

Xshift=X;
Yshift=Y;
Zshift=Z;

Xshift(:,b+1)=X(:,1);
Yshift(:,b+1)=Y(:,1);
Zshift(:,b+1)=Z(:,1);

hold on;
%figure;
surf(Xshift,Yshift,Zshift,'EdgeColor','none','LineStyle','none','FaceColor','y','FaceAlpha',1);
%axis equal



for i=1:100
    for j=1:100
[alphaRG(i,j),betaRG(i,j),Normalr]=rotateNormalToExtraPersonalSpace4('Left',rad2deg(A(i,j)),rad2deg(B(i,j)));
    end
end

X=f*R*cosd(alphaRG).*sind(betaRG);
Y=f*R*sind(alphaRG).*sind(betaRG);
Z=-f*R*cosd(betaRG);

Xshift=X;
Yshift=Y;
Zshift=Z;

Xshift(:,b+1)=X(:,1);
Yshift(:,b+1)=Y(:,1);
Zshift(:,b+1)=Z(:,1);

surf(Xshift,Yshift,Zshift,'EdgeColor','none','LineStyle','none','FaceColor','b','FaceAlpha',1);
%axis equal

end