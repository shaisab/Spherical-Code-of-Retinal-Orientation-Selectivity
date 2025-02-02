close all
clear all

numPoints = 1765;
squareSideLength = sqrt(16.5);
alpha = 0.3;  % Weight for NNRI
beta = 0.7;   % Weight for VDARI
points = generatePointsForBalancedNNRIandVDARI(squareSideLength, numPoints, alpha, beta);
nnri = calculateNNRI(points);
vdari = calculateVDARISquare(points, squareSideLength);

% Display NNRI and VDARI
disp(['NNRI: ', num2str(nnri)]);
disp(['VDARI: ', num2str(vdari)]);

% Plot the points
figure;
scatter(points(:, 1), points(:, 2), 'filled');
title('Generated Points with Balanced NNRI and VDARI');
xlabel('X');
ylabel('Y');
axis equal;

function points = generatePointsForBalancedNNRIandVDARI(squareSideLength, numPoints, alpha, beta)
    % Step 1: Create an initial regular grid of points
    gridSize = ceil(sqrt(numPoints));
    [xGrid, yGrid] = meshgrid(linspace(0, squareSideLength, gridSize));
    points = [xGrid(:), yGrid(:)];
    
    % Trim excess points
    points = points(1:numPoints, :);
    
    % Step 2: Iteratively perturb points to balance NNRI and VDARI
    maxIterations = 100;
    for iter = 1:maxIterations
        for i = 1:numPoints
            % Select a random direction and magnitude for perturbation
            perturbation = (rand(1, 2) - 0.5) * squareSideLength / gridSize;
            candidatePoint = points(i, :) + perturbation;
            
            % Keep candidate point within bounds
            candidatePoint = max(0, min(candidatePoint, squareSideLength));
            
            % Calculate scores
            nnriScore = calculateNNRIForPoint(candidatePoint, points, i);
            vdariScore = calculateLocalVDARIScore(candidatePoint, points, squareSideLength, i);
            
            % Combine scores
            totalScore = alpha * nnriScore + beta * vdariScore;
            
            % Accept the candidate point if it improves the score
            currentNNRIScore = calculateNNRIForPoint(points(i, :), points, i);
            currentVDARIScore = calculateLocalVDARIScore(points(i, :), points, squareSideLength, i);
            currentScore = alpha * currentNNRIScore + beta * currentVDARIScore;
            
            if totalScore > currentScore
                points(i, :) = candidatePoint;
            end
        end
    end
end

function nnriScore = calculateNNRIForPoint(candidatePoint, points, excludeIndex)
    % Calculate NNRI-like score for a candidate point
    tempPoints = points;
    tempPoints(excludeIndex, :) = []; % Exclude the current point
    distances = pdist2(candidatePoint, tempPoints);
    nnriScore = min(distances); % Use minimum distance as the NNRI proxy
end

function localVDARIScore = calculateLocalVDARIScore(candidatePoint, points, squareSideLength, excludeIndex)
    % Calculate VDARI-like score for a candidate point
    tempPoints = points;
    tempPoints(excludeIndex, :) = candidatePoint; % Temporarily replace the current point
    
    [~, v] = voronoin(tempPoints);
    validAreas = [];
    for i = 1:length(v)
        if any(v{i} <= 0) || any(v{i} > size(tempPoints, 1))
            continue;
        end
        domainVertices = tempPoints(v{i}, :);
        if all(domainVertices(:, 1) >= 0 & domainVertices(:, 1) <= squareSideLength & ...
               domainVertices(:, 2) >= 0 & domainVertices(:, 2) <= squareSideLength)
            validAreas = [validAreas; polyarea(domainVertices(:, 1), domainVertices(:, 2))];
        end
    end
    
    if isempty(validAreas) || std(validAreas) == 0
        localVDARIScore = 0;
    else
        localVDARIScore = mean(validAreas) / std(validAreas); % VDARI proxy
    end
end
