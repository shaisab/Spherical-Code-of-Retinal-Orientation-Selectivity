clear all

filedir=cd;
%% Load files
[FileName,PathName] = uigetfile('*.mat','Select pooledmap file for all cluster');

load(FileName);

dr=[0,45,90,135,180,225,270,315];
ar=mean(pooledmap.meanrs);
normalized_ar=ar/max(ar);