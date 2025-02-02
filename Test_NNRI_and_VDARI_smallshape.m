clear all 
close all

% Load files - map overlap file
[fileName, PathName] = uigetfile('*.mat', 'Select map overlap file');
load(fileName);

% Define one polyshapes 
smallShape = map.cluster_polygons{1, 5}; 

% Define the number of points and required valid iterations
numPoints = 1000;                 % Number of points to generate
requiredValidIterations = 10000;  % Number of valid iterations required

% Call the function to get the NNRI distribution
% [nnri_square,vdari_square,nnriDistribution, vdariDistribution] = calculateNNRIAndVDARIForSmallShape(smallShape, numPoints, requiredValidIterations);
[nnri_square,vdari_square,nnriDistribution, vdariDistribution,voronoi_areas, voronoi_vertices,transformedSmallShape] = calculateNNRIAndVDARIForSmallShape(smallShape, numPoints, requiredValidIterations);

% Value to find the percentile of
random_nnri = 1.91;

% Sort the list
sorted_nnriDistribution = sort(nnriDistribution);

% Calculate the percentile rank
percentile_random_nnri = 100 * sum(sorted_nnriDistribution < random_nnri) / numel(nnriDistribution);

% Visualize the NNRI distribution
figure;
edges=[0:0.25:15];
histogram(nnriDistribution, edges);
title('NNRI Distribution');
xlabel('NNRI');
ylabel('Frequency');

valuesToMark = [median(nnriDistribution),random_nnri,nnri_square];  % Values to mark on the histogram
labels = {'Median','Random','Simulated Field'};  % Labels for the values
lineColors = {'red','green','black'};  % Colors for the lines

% Add vertical lines, labels, and customize colors
ylimits = ylim;  % Get the y-axis limits
middleY = (ylimits(2) - ylimits(1)) / 2;  % Middle of the y-axis

for i = 1:length(valuesToMark)
    % Add vertical line
    xline(valuesToMark(i), 'Color', lineColors{i}, 'LineWidth', 2);  
    
     % Create label text with the value included
    labelText = sprintf('%s = %.2f', labels{i}, valuesToMark(i));  % Format the label with the value
    
    % Add label near the middle of the line
    text(valuesToMark(i), middleY, labelText, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'Color', lineColors{i}, 'FontSize', 12, 'FontWeight', 'bold');
end

hold off;

% Value to find the percentile of
random_vdari = 1.91;

% Sort the list
sorted_vdariDistribution = sort(vdariDistribution);

% Calculate the percentile rank
percentile_random_vdari = 100 * sum(sorted_vdariDistribution < random_vdari) / numel(vdariDistribution);

% Visualize the VDARI distribution
figure;
edges=[0:0.25:15];
histogram(vdariDistribution, edges);
title('VDARI Distribution');
xlabel('VDARI');
ylabel('Frequency');

valuesToMark = [median(vdariDistribution),random_vdari,vdari_square];  % Values to mark on the histogram
labels = {'Median','Random','Simulated Field'};  % Labels for the values
lineColors = {'red','green','black'};  % Colors for the lines

% Add vertical lines, labels, and customize colors
ylimits = ylim;  % Get the y-axis limits
middleY = (ylimits(2) - ylimits(1)) / 2;  % Middle of the y-axis

for i = 1:length(valuesToMark)
    % Add vertical line
    xline(valuesToMark(i), 'Color', lineColors{i}, 'LineWidth', 2);  
    
     % Create label text with the value included
    labelText = sprintf('%s = %.2f', labels{i}, valuesToMark(i));  % Format the label with the value
    
    % Add label near the middle of the line
    text(valuesToMark(i), middleY, labelText, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'Color', lineColors{i}, 'FontSize', 12, 'FontWeight', 'bold');
end

hold off;