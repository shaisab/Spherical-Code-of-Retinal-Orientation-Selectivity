
% need to load map.mat file

total=numel(map.isDS);
ONDS=sum(map.isDS.*(map.ONorONOFForOFF==1));
ONOFFDS=sum(map.isDS.*(map.ONorONOFForOFF==2));
OFFDS=sum(map.isDS.*(map.ONorONOFForOFF==3));
ONtDS=sum(map.isDS.*(map.ONorONOFForOFF==4));
output=[total,ONDS,ONOFFDS];

%%
dates={
'Aug_05_2015';
'Aug_20_2014';
'Aug_23_2015';
'Aug_25_2015';
'Aug_27_2014';
'Dec_03_2014';
'Dec_09_2014';
'Feb_17_2015';
'Feb_24_2015';
'Jan_02_2015';
'Jan_06_2015';
'July_14_20152';
'July_14_2015';
'July_17_2015';
'July_22_2015';
'June_10_2015';
'June_18_2015';
'Mar_10_2015';
'Mar_16_2015';
'Nov_04_2014';
'Nov_10_2014';
'Nov_24_2014';
'Nov_26_2014';
'Oct_09_2014'; 
'Sep_12_2014';
'Sep_23_2014'};

j=26;
   
for i=1:size(allCells.cellID,1)
ind(i)=strncmp(allCells.cellID{i},dates(j),11);
end

ONcount=sum(allCells.featureTypeClust(ind==1)==1);
ONOFFcount=sum(allCells.featureTypeClust(ind==1)==2);
dates(j)
[ONcount,ONOFFcount]


