
phys=dlmread('phys');

R=phys(1);
M=phys(2);


%load('allPossibleFieldsAlphaCorr_Calcium imaging data_real_translation_40_DS_created_Oct_09_2015_13_13_ONOFF.mat');
load('allPossibleFieldsAlphaCorr_All Hb9 cells_real_translation_40_DS_created_Nov_25_2015_17_28_ephys.mat');

resolution=floor(370/size(allPossibleFields.matchind180,2));
c=10;
cutoff=(['c',num2str(c)]);
k=1;
h=1;
for alpha=0:resolution:360
for beta=0:resolution:180
    data(k,h)=allPossibleFields.matchind180(k,h).(cutoff);
    k=k+1;
end
k=1;
h=h+1;
end

% plot surface
figure;
alpha=0:resolution:360;
beta=0:resolution:180;
contourf(alpha,beta,data);
hold on;

% down sampling by a factor of 2 in order to see vectors bettern
downFactor=2;
alphaDown=downsample(alpha,downFactor);
betaDown=downsample(beta,downFactor);
dataDown=downsample(data,downFactor);
dataDown=downsample(dataDown',downFactor);
dataDown=dataDown';
zeroheight=zeros(size(dataDown,1),size(dataDown,2));

size(alphaDown)
size(betaDown)
size(dataDown)
% calculating normals
[X,Y]=meshgrid(0:resolution*downFactor:360,0:resolution*downFactor:180);
[xn,yn,zn]=surfnorm(X,Y,dataDown);

% plot normals
h=quiver3(alphaDown,betaDown,zeroheight,xn,yn,dataDown,'AutoScale','off','Color','r','LineWidth',1.5);
h.ShowArrowHead='off';

figure;
% plot normals
hold on;
[alphaDown,betaDown]=meshgrid(alphaDown,betaDown);
dataDown=dataDown/20;
Vx=cosd(alphaDown).*sind(betaDown).*dataDown;
Vy=sind(alphaDown).*sind(betaDown).*dataDown;
Vz=-cosd(betaDown).*dataDown;

surf(R*cosd(alphaDown).*sind(betaDown),R*sind(alphaDown).*sind(betaDown),-R*cosd(betaDown));
h=quiver3(R*cosd(alphaDown).*sind(betaDown),R*sind(alphaDown).*sind(betaDown),-R*cosd(betaDown),Vx,Vy,Vz,'AutoScale','off','Color','r','LineWidth',1.5);
h.ShowArrowHead='off';
axis equal
hold off;

