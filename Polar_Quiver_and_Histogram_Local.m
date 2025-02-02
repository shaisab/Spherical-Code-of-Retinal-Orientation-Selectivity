function [f_quiver_polar_plot_Local, f_polar_histogram_Local]=Polar_Quiver_and_Histogram_Local(pooledmap)
% clear all
% close all
% 
% load('pooledMapAlphaCorr_OS data_40_OSI0.30_1.00_SI_0.00_created_Dec_04_2022_12_53.mat');

%% plot a quiver polar plot for all cells in pooledmap
f_quiver_polar_plot_Local=figure; 
h=quiver(repmat(0,size(pooledmap.clusterdata{1, 56}.UVStvecComp.VecUSt,1),1),repmat(0,size(pooledmap.clusterdata{1, 56}.UVStvecComp.VecUSt,1),1),pooledmap.clusterdata{1, 56}.UVStvecComp.VecUSt,pooledmap.clusterdata{1, 56}.UVStvecComp.VecVSt,1,'k');
set(findobj(gcf, 'type','axes'), 'Visible','off')
set(gcf, 'color', [1 1 1]);
axis equal;
h.ShowArrowHead='off';
%% plot a polar histogram for all cells in pooledmap
f_polar_histogram_Local=figure;
[angle,~]=cart2pol(pooledmap.clusterdata{1, 56}.UVStvecComp.VecUSt,pooledmap.clusterdata{1, 56}.UVStvecComp.VecVSt);
rose(angle,20);
set(f_polar_histogram_Local, 'color', [1 1 1]);
set(f_polar_histogram_Local,'position',[50 50 500 500]);
polarhistogram(angle,20,'FaceColor',[0.29,0.66,0.03],'FaceAlpha',1,'EdgeColor',[0.30,0.75,0.93],'LineWidth',1);
%polarhistogram(angle,40,'FaceColor',[0.6350 0.0780 0.1840],'FaceAlpha',1,'EdgeColor',[0.9290 0.6940 0.1250],'LineWidth',1);
ax = gca;
ax.RTick = [];
ax.ThetaTick =[];
ax.ThetaColor = 'w';
end