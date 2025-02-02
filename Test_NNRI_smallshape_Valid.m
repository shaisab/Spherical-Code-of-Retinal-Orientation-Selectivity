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
nnriDistribution = calculateNNRIForSmallShapeValid(smallShape, numPoints, requiredValidIterations);

% Visualize the NNRI distribution
figure;
histogram(nnriDistribution, 20);
title('NNRI Distribution');
xlabel('NNRI');
ylabel('Frequency');

% Statistical analysis (optional)
% For example, comparing the distribution to a target NNRI of 4.5
[p, h] = signrank(nnriDistribution, 4.5);
disp(['P-value from Wilcoxon signed-rank test: ', num2str(p)]);
