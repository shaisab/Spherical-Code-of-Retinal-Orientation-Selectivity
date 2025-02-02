
mat180=generateMatchindMatrix(allPossibleFields);

resolution=5;
f2=figure;
alpha=0:resolution:360;
beta=0:resolution:180;
contourf(alpha,beta,mat180,10);
hold on;
scatter(allPossibleFields.alphasing,allPossibleFields.betasing,'kx','LineWidth',3);

%% translation
% draw intersection of the saccule plane with the retina
hold on
n=1000;
normal=[0.818464310718389;0.544000000000000;-0.163270242482432];          % normal vector of saccule in global space
[alphaCirc,betaCirc]=GreatCircle(normal,n);
plot(rad2deg(alphaCirc),rad2deg(betaCirc),'k.','LineWidth',2);           % sacula
 
% draw intersection of the utricle plane with the retina
normal=[0.023237162559260;0.229457818946242;0.964288931597075];          % normal vector of utricle in global space
[alphaCirc,betaCirc]=GreatCircle(normal,n);
plot(rad2deg(alphaCirc),rad2deg(betaCirc),'r.','LineWidth',2);            % utricle

%% rotation

%normalASC=[0.596,0.766,0.241];          %ASC
normalASC=[0.646079040832055;0.758000000000000;-0.080292421793915];
%normalASC=[0.604,0.758,0.243]; % ASC symmetric
[alphaASC,betaASC,~]=ABGlobalVisual2RetinaRight(normalASC);
betaiASC=90-abs(90-betaASC);
if alphaASC<=180
    alphaiASC=180+alphaASC;
else alphaASC>180 
    alphaiASC=alphaASC-180;
end

%normalPSC=[0.471,-0.766,0.438];         %PSC
normalPSC=[0.618330720986842;-0.760000000000000;0.193486742398264];
%normalPSC=[0.447,-0.76,0.469];    % PSC symmetric
[alphaPSC,betaPSC,~]=ABGlobalVisual2RetinaRight(normalPSC);
betaiPSC=90-abs(90-betaPSC);
if alphaPSC<=180
    alphaiPSC=180+alphaPSC;
else alphaPSC>180 
    alphaiPSC=alphaPSC-180;
end

%normalLSC=[0.421,-0.060,-0.905];        %LSC
normalLSC=[-0.036837075032666;-0.085000000000000;-0.992592580016110];
%normalLSC=[0.449,-0.085,-0.886];   % LSC symmetric
[alphaLSC,betaLSC,~]=ABGlobalVisual2RetinaRight(normalLSC);
betaiLSC=90-abs(90-betaLSC);
if alphaLSC<=180
    alphaiLSC=180+alphaLSC;
else alphaLSC>180 
    alphaiLSC=alphaLSC-180;
end


hold on;
scatter(alphaASC,betaASC,60,'bx','LineWidth',2);
scatter(alphaiASC,betaiASC,50,'bo','LineWidth',2);
scatter(alphaPSC,betaPSC,60,'gx','LineWidth',2);
scatter(alphaiPSC,betaiPSC,50,'go','LineWidth',2);
scatter(alphaLSC,betaLSC,60,'rx','LineWidth',2);
scatter(alphaiLSC,betaiLSC,50,'ro','LineWidth',2);

