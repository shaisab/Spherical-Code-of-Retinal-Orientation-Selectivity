function [f_quiver_polar_plot_Global, f_polar_histogram_Global]=Polar_Quiver_and_Histogram_Global_local_manual(pooledmap)
filedir=cd;
%% Load files
[FileName,PathName] = uigetfile('*.mat','Select pooledmap file for all clusters');

load(FileName);

tc=pooledmap.clusterdata{1, 1}.typeClust;
pooledmap.clusterdata{1, 1}.UVStvecComp.VecUSt=pooledmap.clusterdata{1, 1}.UVStvecComp.VecUSt(tc~=6);
pooledmap.clusterdata{1, 1}.UVStvecComp.VecVSt=pooledmap.clusterdata{1, 1}.UVStvecComp.VecVSt(tc~=6);
%% plot a quiver polar plot for all cells in pooledmap
f_quiver_polar_plot_Global=figure; 
h=quiver(repmat(0,size(pooledmap.clusterdata{1, 1}.UVStvecComp.VecUSt,1),1),repmat(0,size(pooledmap.clusterdata{1, 1}.UVStvecComp.VecUSt,1),1),pooledmap.clusterdata{1, 1}.UVStvecComp.VecUSt,pooledmap.clusterdata{1, 1}.UVStvecComp.VecVSt,1,'k');
set(findobj(gcf, 'type','axes'), 'Visible','off')
set(gcf, 'color', [1 1 1]);
axis equal;
h.ShowArrowHead='off';

f_polar_histogram_Global=figure;
[angle,~]=cart2pol(pooledmap.clusterdata{1, 1}.UVStvecComp.VecUSt,pooledmap.clusterdata{1, 1}.UVStvecComp.VecVSt);
rose(angle,20);
set(f_polar_histogram_Global, 'color', [1 1 1]);
set(f_polar_histogram_Global,'position',[50 50 500 500]);
polarhistogram(angle,20,'FaceColor',[0.29,0.66,0.03],'FaceAlpha',1,'EdgeColor',[0.30,0.75,0.93],'LineWidth',1);
%polarhistogram(angle,40,'FaceColor',[0.6350 0.0780 0.1840],'FaceAlpha',1,'EdgeColor',[0.9290 0.6940 0.1250],'LineWidth',1);
ax = gca;
ax.RTick = [];
ax.ThetaTick =[];
ax.ThetaColor = 'w';
end