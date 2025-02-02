
clear all
[names,PathName] = uigetfile('*.mat','Select allPossibleFields ONOFF, uniform translation, and unifrom rotation probed with translation and rotation','MultiSelect','on');


for i=1:size(names,2)  
    if strfind(names{i},'ONOFF') * strfind(names{i},'real_translation')
        load(num2str(names{i}));
        ONOFF.trans=generateMatchindMatrix(allPossibleFields);  
    elseif strfind(names{i},'ONOFF') * strfind(names{i},'real_rotation')
        load(num2str(names{i}));
        ONOFF.rot=generateMatchindMatrix(allPossibleFields); 
    elseif strfind(names{i},'translation_uniform') * strfind(names{i},'real_translation')
        load(num2str(names{i}));
        uniform.trans.probedTrans=generateMatchindMatrix(allPossibleFields); 
    elseif strfind(names{i},'translation_uniform') * strfind(names{i},'single channel probed with_rotation')
        load(num2str(names{i}));
        uniform.trans.probedRot=generateMatchindMatrix(allPossibleFields); 
    elseif strfind(names{i},'rotation_uniform') * strfind(names{i},'real_rotation')
        load(num2str(names{i}));
        uniform.rot.probedRot=generateMatchindMatrix(allPossibleFields); 
    elseif strfind(names{i},'rotation_uniform') * strfind(names{i},'single channel probed with_translation')
        load(num2str(names{i}));
        uniform.rot.probedTrans=generateMatchindMatrix(allPossibleFields); 
    end
    
end

resolution=370/size(allPossibleFields.matchind180,2);
alpha=0:resolution:360;
beta=0:resolution:180;

figure;
subplot(2,2,1)
ONOFF_transProbedTrans=ONOFF.trans.*uniform.trans.probedTrans;
contourf(alpha,beta,ONOFF_transProbedTrans);
%caxis([-26 20]);
title('ONOFF trans - transProbedTrans');
subplot(2,2,3)
ONOFF_rotProbedTrans=ONOFF.trans.*uniform.rot.probedTrans;
contourf(alpha,beta,ONOFF_rotProbedTrans);
%caxis([-26 20]);
title('ONOFF trans - rotProbedTrans');
subplot(2,2,2)
ONOFF_transProbedRot=ONOFF.rot.*uniform.trans.probedRot;
contourf(alpha,beta,ONOFF_transProbedRot);
%caxis([-26 20]);
title('ONOFF rot - transProbedRot');
subplot(2,2,4)
ONOFF_rotProbedRot=ONOFF.rot.*uniform.rot.probedRot;
contourf(alpha,beta,ONOFF_rotProbedRot);
%caxis([-26 20]);
title('ONOFF rot - rotProbedRot');


    
    
    
    
    
    