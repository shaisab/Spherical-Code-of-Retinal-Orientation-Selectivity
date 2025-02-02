


% find max singularities in global maps using the final azimuth and
% elevation conventions
[allPossibleFieldsName,PathName] = uigetfile('*.mat','Select a allPossibleFields translation or rotation in reitinal space');
load(allPossibleFieldsName);

if strfind(allPossibleFieldsName,'inverse polarity')
% for inverse polarity global maps
if allPossibleFields.maxAziR>180
allPossibleFields.maxAziG=180-allPossibleFields.maxAziR
else
allPossibleFields.maxAziG=180-allPossibleFields.maxAziR    
end
allPossibleFields.maxElevG=90-allPossibleFields.maxElevR    
    
else
% for regular polarity global maps
if allPossibleFields.maxAziR>180
allPossibleFields.maxAziG=360-allPossibleFields.maxAziR
else
allPossibleFields.maxAziG=-allPossibleFields.maxAziR    
end
allPossibleFields.maxElevG=allPossibleFields.maxElevR-90
end

save(allPossibleFieldsName,'allPossibleFields')
