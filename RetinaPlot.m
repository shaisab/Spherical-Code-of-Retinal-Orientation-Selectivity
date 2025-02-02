function f4 = RetinaPlot(f4) 

discPoints=40;

%% Physical Parameters.
% This data is imported from data exported from the mapping script. These
% values correspond to measured values on the flattened retina.

phys=dlmread('phys');

% for ON and ONOFF DS cells
R=phys(1)*1.1;
M=phys(2)*1.1; %This should be normalized to 1.


% for Hb9 cells
% R=phys(1);
% M=phys(2); %This should be normalized to 1.

m1=phys(3);
m2=phys(4);
m3=phys(5);
m4=phys(6);

nu=phys(7);

a1=phys(8);
a2=phys(9);
a3=phys(10);
a4=phys(11);

s=linspace(0,M,100);
th=linspace(0,2*pi,100);

[S,TH]=meshgrid(s,th);

X=R*cos(TH).*sin(S/R);
Y=R*sin(TH).*sin(S/R);
Z=R*(1-cos(S/R));
%figure
surf(X,Y,Z,S,'EdgeColor','none','FaceAlpha',0.3);
hold on
surf(X.*1.02,Y.*1.02,Z.*1.02,S,'EdgeColor','none','FaceAlpha',1,'CDataMode','auto');        % for ONOFF and ON cells
%surf(X.*1,Y.*1,Z.*1,S,'EdgeColor','none','FaceAlpha',1,'CDataMode','auto');        % for illustration of two orthogonal translation channels

%surf(X,Y,Z,S, 'FaceColor',[0.8,0.8,0.8],'EdgeColor','none','FaceAlpha',0.7);
axis([-R,R,-R,R,0,2*R]);
colormap(flipud(gray))
axis off;
axis equal;
shading interp;
hold on;


end