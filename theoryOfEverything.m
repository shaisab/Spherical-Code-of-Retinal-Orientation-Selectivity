

function []= theoryOfEverything()

%% Load files and extract meanOutput matrices
filedir=cd;
names=getFileList(filedir,'.mat',0,'anywhere');     % detect all .mat files
Right=[];
type={'translation','rotation'};
for t=1:2
    typestr=type{t};
    for i=1:size(names,2)
        if strfind(names{i},typestr)
            if strfind(names{i},'Right')
                if strfind(names{i},'sing1')
                    load(num2str(names{i}));
                    Right.(typestr){1,1}=allPossibleFields.meanOutput;
                elseif strfind(names{i},'sing2')
                    load(num2str(names{i}));
                    Right.(typestr){1,2}=allPossibleFields.meanOutput;
                elseif strfind(names{i},'sing3')
                    load(num2str(names{i}));
                    Right.(typestr){1,3}=allPossibleFields.meanOutput;
                elseif strfind(names{i},'sing4')
                    load(num2str(names{i}));
                    Right.(typestr){1,4}=allPossibleFields.meanOutput;
                end
            elseif strfind(names{i},'Left')
                if strfind(names{i},'sing1')
                    load(num2str(names{i}));
                    Left.(typestr){1,1}=allPossibleFields.meanOutput;
                elseif strfind(names{i},'sing2')
                    load(num2str(names{i}));
                    Left.(typestr){1,2}=allPossibleFields.meanOutput;
                elseif strfind(names{i},'sing3')
                    load(num2str(names{i}));
                    Left.(typestr){1,3}=allPossibleFields.meanOutput;
                elseif strfind(names{i},'sing4')
                    load(num2str(names{i}));
                    Left.(typestr){1,4}=allPossibleFields.meanOutput;
                end
            end
        end
    end
end

%% calculate standard deviation for pre-synaptic channels for the right eye (sd for the left eye will be identical) 
for c=1:size(Right.translation,2)
sdTright(c)=std(Right.translation{1,c}(:),'omitnan');
sdRright(c)=std(Right.rotation{1,c}(:),'omitnan');
rangeTright(c)=max(Right.translation{1,c}(:),[],'omitnan')-min(Right.translation{1,c}(:),[],'omitnan');
rangeRright(c)=max(Right.rotation{1,c}(:),[],'omitnan')-min(Right.rotation{1,c}(:),[],'omitnan');
end
range=[rangeTright;rangeRright];
%% plot bar graphs of std for pre-synaptic channels
groupedBars(1:4,[rangeTright(1),rangeTright(3),rangeTright(2),rangeTright(4)],[rangeRright(1),rangeRright(3),rangeRright(2),rangeRright(4)]);
set(gcf, 'color', [1 1 1]);
box off
ax = gca;
ax.FontSize=14;
ax.XTick = 1:4;
ylim([0 1.4]);
ax.XTickLabel={'N','T','D','V'};
xlabel('Pre-synaptic channel','FontSize',16);
ylabel('Tuning index','FontSize',16); 
legend('during translation','during rotation','northeast');
legend('boxoff');


%% generate postsynaptic cells for translation

cellT{1,1}=Right.translation{1,2}-Right.translation{1,4}+Left.translation{1,2}-Left.translation{1,4};   % cellT t1 - translate up - vertical axis
cellT{1,2}=-Right.translation{1,2}+Right.translation{1,4}-Left.translation{1,2}+Left.translation{1,4};  % cellT t2 - translate down - vertical axis
cellT{1,3}=Right.translation{1,1}-Right.translation{1,3}+Left.translation{1,1}-Left.translation{1,3};                                               % cellT t3 - translate forward L - Left 45 axis
cellT{1,4}=-Right.translation{1,1}+Right.translation{1,3}-Left.translation{1,1}+Left.translation{1,3};                                              % cellT t4 - translate backward R - Left 45 axis
cellT{1,5}=0.*(Left.translation{1,1}-Left.translation{1,3});                                                 % cellT t5 - translate forward R - Right 45 axis
cellT{1,6}=0.*(-Left.translation{1,1}+Left.translation{1,3});                                                % cellT t6 - translate backward L - Right 45 axis

cellT{1,7}=Right.translation{1,1}-Right.translation{1,3}-Left.translation{1,1}+Left.translation{1,3};   % cellT r1 - rotate left - vertical axis
cellT{1,8}=-Right.translation{1,1}+Right.translation{1,3}+Left.translation{1,1}-Left.translation{1,3};  % cellT r2 - rotate right - vertical axis
cellT{1,9}=Right.translation{1,2}-Left.translation{1,2};%+cellT{1,4};                         % cellT r3 - pitch backward - Left 45 axis
cellT{1,10}=Right.translation{1,4}-Left.translation{1,4};%+cellT{1,3};                       % cellT r4 - pitch forward - Left 45 axis
cellT{1,11}=-Right.translation{1,2}+Left.translation{1,2};                        % cellT r5 - pitch backward - Right 45 axis
cellT{1,12}=-Right.translation{1,4}+Left.translation{1,4};                        % cellT r6 - pitch forward - Right 45 axis

cellID={'t1 - translate up - vertical axis';'t2 - translate down - vertical axis';'t3 - translate forward L - Left 45 axis';...
    't4 - translate backward R - Left 45 axis';'t5 - translate forward R - Right 45 axis';'t6 - translate backward L - Right 45 axis';...
    'r1 - rotate left - vertical axis';'r2 - rotate right - vertical axis';'r3 - pitch backward - Left 45 axis';...
    'r4 - pitch forward - Left 45 axis';'r5 - pitch backward - Right 45 axis';'r6 - pitch forward - Right 45 axis'};

%% calculate standard deviation and CV for translation cells
for c=1:size(cellT,2)
minT{1,c}=repmat(min(cellT{1,c}(:)),size(cellT{1,c},1),size(cellT{1,c},2));         % find the minimum of each map
cellTa{1,c}=cellT{1,c}-minT{1,c};       %adjust map values to start from 0
cvT(c)=std(cellTa{1,c}(:),'omitnan')/mean(cellTa{1,c}(:),'omitnan');
sdT(c)=std(cellT{1,c}(:),'omitnan');
rangeT(c)=max(cellT{1,c}(:),[],'omitnan')-min(cellT{1,c}(:),[],'omitnan');
end
    
%% plot translation cells
resolution=1;
azi=-180:resolution:180;
elev=-90:resolution:90;
figure('units','normalized','outerposition',[0 0 1 1]);
hold on
for j=1:size(cellT,2)
ax=subplot(4,3,j);
[C,h]=contourf(azi,elev,cellT{1,j}); 
h.LevelList=linspace(-2,2,9)';
hold on
plot([-180 180],[0 0],'m','LineWidth',1);
plot([0 0],[-90 90],'m','LineWidth',1);
xlim([-180 180]);ylim([-90 90]);
colormap jet
caxis([-2 2]);
grid on
ax.XTick = -180:45:180;
ax.YTick = -90:45:90;
ax.XTickLabel = {'behind','-135','left','-45','ahead','45','right','135','behind'};
ax.FontSize=12;
%xlabel('Azimuth (deg.)','FontSize',16);
%ax.YTickLabel = {'','','','',''};
hold on
%ylabh=ylabel('Elevation (deg.)','FontSize',16);
%set(ylabh, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
ax.YTickLabel = {'nadir','-45', 'horizon','45','zenith'};
title(num2str(cellID{j,1}));
hold off
end

%% generate postsynaptic cells
cellR{1,1}=Right.rotation{1,2}-Right.rotation{1,4}+Left.rotation{1,2}-Left.rotation{1,4};       % cellR t1 - translate up - vertical axis
cellR{1,2}=-Right.rotation{1,2}+Right.rotation{1,4}-Left.rotation{1,2}+Left.rotation{1,4};      % cellR t2 - translate down - vertical axis
cellR{1,3}=Right.rotation{1,1}-Right.rotation{1,3}+Left.rotation{1,1}-Left.rotation{1,3};;                                              % cellR t3 - translate forward L - Left 45 axis
cellR{1,4}=-Right.rotation{1,1}+Right.rotation{1,3}-Left.rotation{1,1}+Left.rotation{1,3};;                                             % cellR t4 - translate backward R - Left 45 axis
cellR{1,5}=0.*(Left.rotation{1,1}-Left.rotation{1,3});                                                % cellR t5 - translate forward R - Right 45 axis
cellR{1,6}=0.*(-Left.rotation{1,1}+Left.rotation{1,3});                                               % cellR t6 - translate backward L - Right 45 axis

cellR{1,7}=Right.rotation{1,1}-Right.rotation{1,3}-Left.rotation{1,1}+Left.rotation{1,3};       % cellR r1 - rotate left - vertical axis
cellR{1,8}=-Right.rotation{1,1}+Right.rotation{1,3}+Left.rotation{1,1}-Left.rotation{1,3};      % cellR r2 - rotate right - vertical axis
cellR{1,9}=Right.rotation{1,2}-Right.rotation{1,4}-Left.rotation{1,2};%+cellR{1,4};                          % cellR r3 - pitch backward - Left 45 axis
cellR{1,10}=-Right.rotation{1,2}+Right.rotation{1,4}-Left.rotation{1,4};%+cellR{1,3};                        % cellR r4 - pitch forward - Left 45 axis
cellR{1,11}=-Right.rotation{1,2}+Left.rotation{1,2}-Left.rotation{1,4};                        % cellT r5 - pitch backward - Right 45 axis
cellR{1,12}=-Right.rotation{1,4}-Left.rotation{1,2}+Left.rotation{1,4};                        % cellT r6 - pitch forward - Right 45 axis

%% calculate standard deviation and CV for rotation cells
for c=1:size(cellR,2)
minR{1,c}=repmat(min(cellR{1,c}(:)),size(cellR{1,c},1),size(cellR{1,c},2));         % find the minimum of each map
cellRa{1,c}=cellR{1,c}-minR{1,c};       %adjust map values to start from 0
cvR(c)=std(cellRa{1,c}(:),'omitnan')/mean(cellRa{1,c}(:),'omitnan');
sdR(c)=std(cellR{1,c}(:),'omitnan');
rangeR(c)=max(cellR{1,c}(:),[],'omitnan')-min(cellR{1,c}(:),[],'omitnan');
end

%% plot bar graphs of std
%groupedBars(1:10,sdT([1:4,7:12]),sdR([1:4,7:12]));
groupedBars(1:10,rangeT([1:4,7:12]),rangeR([1:4,7:12]));
set(gcf, 'color', [1 1 1]);
box off
ax = gca;
ax.FontSize=14;
ax.XTick = 1:10;
ax.XTickLabel={'T1','T2','T3','T4','R1','R2','R3','R4','R5','R6'};
xlabel('Post-synaptic cell','FontSize',16);
ylabel('Tuning index','FontSize',16); 
legend('during translation','during rotation','northeast');
legend('boxoff');

%groupedBars(1:10,cvT([1:4,7:12]),cvR([1:4,7:12]))

%% plot rotation cells
resolution=1;
azi=-180:resolution:180;
elev=-90:resolution:90;
figure('units','normalized','outerposition',[0 0 1 1]);
hold on
for j=1:size(cellR,2)
ax=subplot(4,3,j);
[C,h]=contourf(azi,elev,cellR{1,j}); 
h.LevelList=linspace(-2,2,9)';
hold on
plot([-180 180],[0 0],'m','LineWidth',1);
plot([0 0],[-90 90],'m','LineWidth',1);
xlim([-180 180]);ylim([-90 90]);
colormap jet
caxis([-2 2]);
grid on
ax.XTick = -180:45:180;
ax.YTick = -90:45:90;
ax.XTickLabel = {'behind','-135','left','-45','ahead','45','right','135','behind'};
ax.FontSize=12;
%xlabel('Azimuth (deg.)','FontSize',16);
%ax.YTickLabel = {'','','','',''};
hold on
%ylabh=ylabel('Elevation (deg.)','FontSize',16);
%set(ylabh, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
ax.YTickLabel = {'nadir','-45', 'horizon','45','zenith'};
title(num2str(cellID{j,1}));
hold off
end

save('summary_postsynaptic cells.mat','cellR','cellT');

end