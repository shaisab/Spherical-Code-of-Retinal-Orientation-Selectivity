
function VestibularOcularInGlobalSpace

%% input vectors

OptL=[0.406450649610760,-0.833347328309341,0.374606593415912];    % opt axis
ASC=[0.596,0.766,0.241];          %ASC
[alphaASC,betaASC]=AB(normalASC);
alphaASC=rad2deg(alphaASC);
betaASC=rad2deg(betaASC);
betaiASC=90-abs(90-betaASC);
if alphaASC<=180
    alphaiASC=180+alphaASC;
else alphaASC>180 
    alphaiASC=alphaASC-180;
end

normalPSC=[0.471,-0.766,0.438];         %PSC
[alphaPSC,betaPSC]=AB(normalPSC);
alphaPSC=rad2deg(alphaPSC);
betaPSC=rad2deg(betaPSC);
betaiPSC=90-abs(90-betaPSC);
if alphaPSC<=180
    alphaiPSC=180+alphaPSC;
else alphaPSC>180 
    alphaiPSC=alphaPSC-180;
end

normalLSC=[0.421,-0.060,-0.905];        %LSC
[alphaLSC,betaLSC]=AB(normalLSC);
alphaLSC=rad2deg(alphaLSC);
betaLSC=rad2deg(betaLSC);
betaiLSC=90-abs(90-betaLSC);
if alphaLSC<=180
    alphaiLSC=180+alphaLSC;
else alphaLSC>180 
    alphaiLSC=alphaLSC-180;
end

allAlpha=[alphaOptL,alphaASC,alphaPSC,alphaLSC,alphaiASC,alphaiPSC,alphaiLSC];
allBeta=[betaOptL,betaASC,betaPSC,betaLSC,betaiASC,betaiPSC,betaiLSC];

%% process beta to confrom to global space convetions
for i=1:size(allBeta,2)
    if  allBeta(i)==0
        allBeta(i)=180;
    else
        allBeta(i)=mod(180-allBeta(i),180);
    end
end


%% calculate azimuth and elevation of canals in global space
for i=1:size(allAlpha,2)
[azimuthRight(i),elevationRight(i)]=rotateNormalToExtraPersonalSpace3('Right',allAlpha(i),allBeta(i));
[azimuthLeft(i),elevationLeft(i)]=rotateNormalToExtraPersonalSpace3('Left',allAlpha(i),allBeta(i));
end

%% plot canals in global space
figure;
scatter(azimuthRight(1),elevationRight(1),'k'); % opt axis Right
hold on
scatter(azimuthLeft(1),elevationLeft(1),'k');   % opt axis Left
hold on
scatter(azimuthRight(2:4),elevationRight(2:4),'b'); % 3 canals Right
hold on
scatter(azimuthLeft(2:4),elevationLeft(2:4),'r');   % 3 canals Left
xlim([0 360]); ylim([0 180]);





end