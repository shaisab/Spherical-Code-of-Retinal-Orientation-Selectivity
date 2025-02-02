
s=size(map.isDS,1);
map.isDS=ones(s,1);
map.DSI=repmat(0.5,s,1);
save('mapAlphaCorr_Aug_25_2015.mat','map');


s=size(ONorONOFForOFFtmp,1);
ONorONOFForOFFtmp=repmat(1,s,1);