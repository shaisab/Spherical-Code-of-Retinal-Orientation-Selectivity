
filename=('allPossibleFieldsAlphaCorr_single channel probed with_rotation_40_DS_created_Nov_03_2015_08_43_uniform_Right_sing1.mat');
load(filename,'allPossibleFields');

f1=figure;
alpha=-180:5:180;
beta=-90:5:90;
contourf(alpha,beta,allPossibleFields.meanOutput); 
xlim([-180 180]);ylim([-90 90]);

allPossibleFields.meanOutput=fliplr(allPossibleFields.meanOutput);
allPossibleFields.mat180single=fliplr(allPossibleFields.mat180single);

f2=figure;
alpha=-180:5:180;
beta=-90:5:90;
contourf(alpha,beta,allPossibleFields.meanOutput); 
xlim([-180 180]);ylim([-90 90]);

save(filename,'allPossibleFields');
clear all
