%%% code for generating matchind matrices and plotting contour plots
c=10;
cutoff=(['c',num2str(c)]);
k=1;
h=1;
% for alpha=0:5:355
% for beta=90:5:180
for alpha=0:10:350
for beta=90:10:180
    mat180(k,h)=allPossibleFields.matchind180(k,h).(cutoff);
    mat90(k,h)=allPossibleFields.matchind90(k,h).(cutoff);
    k=k+1;
end
k=1;
h=h+1;
end

normalASC=[0.592,0.764,0.247];     %ASC
[alphaASC,betaASC]=AB(normalASC);
normalPSC=[0.470,-0.764,0.438];       %PSC
[alphaPSC,betaPSC]=AB(normalPSC);
normalLSC=[0.421,-0.056,-0.901];       %LSC
[alphaLSC,betaLSC]=AB(normalLSC);

%%%%%%%%%%%%% standardize ASC %%%%%%%%%%%%%%%%%%%%%%%%%%
if rad2deg(betaASC)>90
    betaASC=rad2deg(betaASC);
    if rad2deg(alphaASC)<=0
    alphaASC=360+rad2deg(alphaASC);
    else
        alphaASC=rad2deg(alphaASC);
    end
else
    betaASC=180-rad2deg(betaASC);
    if rad2deg(alphaASC)<=180
        alphaASC=180+rad2deg(alphaASC);
    else
        alphaASC=rad2deg(alphaASC)-180;
    end
end
% standardize PSC
if rad2deg(betaPSC)>90
    betaPSC=rad2deg(betaPSC);
    if rad2deg(alphaPSC)<=0
    alphaPSC=360+rad2deg(alphaPSC);
    else
        alphaPSC=rad2deg(alphaPSC);
    end
else
    betaPSC=180-rad2deg(betaPSC);
    if rad2deg(alphaPSC)<=180
        alphaPSC=180+rad2deg(alphaPSC);
    else
        alphaPSC=rad2deg(alphaPSC)-180;
    end
end
% standardize LSC
if rad2deg(betaLSC)>90
    betaLSC=rad2deg(betaLSC);
    if rad2deg(alphaLSC)<=0
    alphaLSC=360+rad2deg(alphaLSC);
    else
        alphaLSC=rad2deg(alphaLSC);
    end
else
    betaLSC=180-rad2deg(betaLSC);
    if rad2deg(alphaLSC)<=180
        alphaLSC=180+rad2deg(alphaLSC);
    else
        alphaLSC=rad2deg(alphaLSC)-180;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f2=figure;
alpha=0:10:350;
beta=90:10:180;
subplot(1,2,1);
contourf(alpha,beta,mat180);
hold on;
scatter(alphaASC,betaASC,'bx','LineWidth',3);
scatter(alphaPSC,betaPSC,'gx','LineWidth',3);
scatter(alphaLSC,betaLSC,'rx','LineWidth',3);
hold off;
subplot(1,2,2);
contourf(alpha,beta,mat90);
hold on;
scatter(alphaASC,betaASC,'bx','LineWidth',3);
scatter(alphaPSC,betaPSC,'gx','LineWidth',3);
scatter(alphaLSC,betaLSC,'rx','LineWidth',3);
hold off;


currenttime=datestr(now,'mmm_dd_yyyy_HH_MM');
savefig(f2,['allPossibleFields_','_cutoff_',num2str(c),'_DS_created_',currenttime]);
