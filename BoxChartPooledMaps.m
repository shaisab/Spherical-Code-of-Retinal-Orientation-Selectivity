clear all
close all
clc

load('pooledMapAlphaCorr_OS data_40_OSI0.20_1.00_SI_0.00_created_Dec_08_2022_13_38_All_Cells_20_bins.mat');
pooledmap_all_cells=pooledmap;
load('pooledMapAlphaCorr_OS data_40_OSI0.20_1.00_SI_0.00_created_Dec_04_2022_09_02_On_Sus.mat');
pooledmap_on_sus=pooledmap;
load('pooledMapAlphaCorr_OS data_40_OSI0.20_1.00_SI_0.00_created_Dec_06_2022_10_26_On_Trans.mat');
pooledmap_on_trans=pooledmap;
load('pooledMapAlphaCorr_OS data_40_OSI0.20_1.00_SI_0.00_created_Dec_06_2022_11_13_On_Off_Sus.mat');
pooledmap_on_off_sus=pooledmap;
load('pooledMapAlphaCorr_OS data_40_OSI0.20_1.00_SI_0.00_created_Dec_06_2022_11_55_On_Off_Trans.mat');
pooledmap_on_off_trans=pooledmap;


OSI_Group=padcat(pooledmap_all_cells.OSI,pooledmap_on_sus.OSI,pooledmap_on_trans.OSI,pooledmap_on_off_sus.OSI,pooledmap_on_off_trans.OSI);
SI_Group=padcat(pooledmap_all_cells.symmetry_ratio,pooledmap_on_sus.symmetry_ratio,pooledmap_on_trans.symmetry_ratio,pooledmap_on_off_sus.symmetry_ratio,pooledmap_on_off_trans.symmetry_ratio);
R2_Group=padcat(pooledmap_all_cells.oriR2,pooledmap_on_sus.oriR2,pooledmap_on_trans.oriR2,pooledmap_on_off_sus.oriR2,pooledmap_on_off_trans.oriR2);
type = ["All types","ON Sustained","ON Transient","ON-OFF Sustained","ON-OFF Transient"];
% data.OSI_Group=OSI_Group;
% data.SI_Group=SI_Group;
% data.type=type;

figure
hold on

boxplot(OSI_Group,type,'symbol', '')
ylabel('OSI')
xlabel('Type')

hold off 

figure
hold on

boxplot(SI_Group,type,'symbol', '')
ylabel('SI')
xlabel('Type')

hold off

% figure
% hold on
% 
% boxplot(R2_Group,type,'symbol', '')
% ylabel('R^2')
% xlabel('Type')
% 
% hold off

for i=1:5
    figure
    histogram(OSI_Group(:,i));
    title(['OSI',type(i)]);
end

for i=1:5
    figure
    histogram(SI_Group(:,i));
    title(['SI',type(i)]);
end

% for i=1:5
%     figure
%     histogram(R2_Group(:,i));
%     title(['R^2',type(i)]);
% end

for i=1:5
    figure
    scatter(OSI_Group(:,i),SI_Group(:,i));
    title(type(i));
    xlabel('OSI')
    ylabel('SI')
end

% x=mean(OSI_Group,'omitnan');
% y=mean(SI_Group,'omitnan');
% z=mean(R2_Group,'omitnan');

% scatter3(x,y,z);
% text(x(1),y(1),z(1),'All types')
% text(x(2),y(2),z(2),'ON Sustained')
% text(x(3),y(3),z(3),'ON Transient')
% text(x(4),y(4),z(4),'ON-OFF Sustained')
% text(x(5),y(5),z(5),'ON-OFF Transient')
% xlabel('OSI')
% ylabel('SI')
% zlabel('R^2')

% x=mean(OSI_Group,'omitnan');
% y=mean(SI_Group,'omitnan');
% 
% scatter3(x,y);
% text(x(1),y(1),'All types')
% text(x(2),y(2),'ON Sustained')
% text(x(3),y(3),'ON Transient')
% text(x(4),y(4),'ON-OFF Sustained')
% text(x(5),y(5),'ON-OFF Transient')
% xlabel('OSI')
% ylabel('SI')