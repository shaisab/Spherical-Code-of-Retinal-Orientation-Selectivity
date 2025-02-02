
Icolmax=[2,20,38,56];
Irowmax=[25,23,15,17];
[Icolmax,J]=sort(Icolmax);      % sort Icolmax to assend along alpha
Irowmax=Irowmax(J);             % sort Irowmax to assend along alpha

resolution=5;
alpha=0:resolution:360;
beta=0:resolution:180;
allPossibleFields.alphasing=alpha(Icolmax);
allPossibleFields.betasing=beta(Irowmax);
allPossibleFields.Icolmax=Icolmax;
allPossibleFields.Irowmax=Irowmax;

save('allPossibleFieldsAlphaCorr_Calcium imaging data_real_rotation_40_DS_created_Nov_12_2015_11_44_ON_smoothed.mat','allPossibleFields');
