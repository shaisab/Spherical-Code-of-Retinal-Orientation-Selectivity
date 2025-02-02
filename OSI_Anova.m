clear all
close all

filedir=cd;
%% Load files
[FileName,PathName] = uigetfile('*.mat','Select pooledmap file for all clusters');

load(FileName);

OSI_ON_SUS=pooledmap.OSI(pooledmap.typeClust==3);
OSI_ON_TRANS=pooledmap.OSI(pooledmap.typeClust==8);
OSI_ONOFF_SUS=pooledmap.OSI(pooledmap.typeClust==4);
OSI_ONOFF_TRANS=pooledmap.OSI(pooledmap.typeClust==6);

[p,tbl,stats] = anova1(pooledmap.OSI,pooledmap.typeClust)
[c,m,h,gnames] = multcompare(stats)