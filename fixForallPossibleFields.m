
filename=('allPossibleFieldsAlphaCorr_Calcium imaging data_single channel probed with_rotation_alpha_10_beta_120_40_DS_Nov_06_2015_16_56_global space_Left_sing1');
load(filename);

resolution=1;
elements=floor(size(allPossibleFields.meanOutput,2)/2);
dOutput=circshift(allPossibleFields.meanOutput,[0 elements]);   %shift azimuth 0 to middle of graph
dOutput(:,elements+1)=[];
dOutput=[dOutput(:,end),dOutput];
allPossibleFields.meanOutput=fliplr(dOutput);
f2=figure;
hold on
alpha=-180:resolution:180;
beta=-90:resolution:90;
contourf(alpha,beta,allPossibleFields.meanOutput); 
xlim([-180 180]);ylim([-90 90]);

save(filename,'allPossibleFields');
savefig(f2,filename);
clear all