
nTotal=size(map.ONorONOFForOFF,1);
nONDS=sum(logical((map.isDS==1) .* (map.ONorONOFForOFF==1)));
nONOFFDS=sum(logical((map.isDS==1) .* (map.ONorONOFForOFF==2)));
nOFFDS=sum(logical((map.isDS==1) .* (map.ONorONOFForOFF==3)));
nONtDS=sum(logical((map.isDS==1) .* (map.ONorONOFForOFF==4)));

summary=[nTotal,nONDS,nONOFFDS,nOFFDS,nONtDS];
