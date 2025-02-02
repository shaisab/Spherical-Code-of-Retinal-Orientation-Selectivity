clear all 
close all

% Load files - map overlap file
[fileName, PathName] = uigetfile('*.mat', 'Select map overlap file');
load(fileName);

% Define two polyshapes (these could be any arbitrary polygons)
largeShape = map.combinedShape;  % Larger square shape
smallShape = map.cluster_polygons{1, 5};  % Smaller square shape

% Set parameters
largeArea = 17;   % mm^2 for large shape
smallArea = 0.28; % mm^2 for small shape
numPoints = 700;  % Number of points to generate
nnriTarget = 4.5; % Target NNRI value
numIterations = 10;  % Number of random placements

% Call the function
nnriDistribution = calculateNNRIForShapes(largeShape, smallShape, largeArea, smallArea, numPoints, nnriTarget, numIterations);
