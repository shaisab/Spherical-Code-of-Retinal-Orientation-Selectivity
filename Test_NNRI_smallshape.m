clear all 
close all

% Load files - map overlap file
[fileName, PathName] = uigetfile('*.mat', 'Select map overlap file');
load(fileName);

% Define one polyshapes 
smallShape = map.cluster_polygons{1, 5}; 

% Define the number of points and iterations
numPoints = 1000;        % Number of points to generate
numIterations = 10000;   % Number of iterations

% Call the function to get the NNRI distribution
nnriDistribution = calculateNNRIForSmallShape(smallShape, numPoints, numIterations);

% Visualize the NNRI distribution
figure;
histogram(nnriDistribution, 20);
title('NNRI Distribution');
xlabel('NNRI');
ylabel('Frequency');

% Perform the Wilcoxon signed-rank test to check if the NNRI values are significantly different from 4.5
[p, h] = signrank(nnriDistribution, 4.5);

% Display the results
if h == 0
    disp(['The NNRI distribution is not significantly different from 4.5 (p = ', num2str(p), ')']);
else
    disp(['The NNRI distribution is significantly different from 4.5 (p = ', num2str(p), ')']);
end

