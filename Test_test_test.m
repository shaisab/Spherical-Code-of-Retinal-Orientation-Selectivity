close all
clear all

numPoints = 1765;
squareSideLength = sqrt(16.5);
alpha = 0.3;  % Weight for NNRI
beta = 0.7;   % Weight for VDARI
points = generatePointsWithBalancedNNRIandVDARI(squareSideLength, numPoints, alpha, beta);
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

% Function to generate points balancing NNRI and VDARI
function points = generatePointsWithBalancedNNRIandVDARI(squareSideLength, numPoints, alpha, beta)
    points = [];
    for i = 1:numPoints
        if i == 1
            points = [rand(1, 1) * squareSideLength, rand(1, 1) * squareSideLength];
        else
            bestPoint = [];
            bestScore = -inf;  % Start with a very low score
            
            for attempt = 1:10  % Max attempts per point
                candidatePoint = [rand * squareSideLength, rand * squareSideLength];
                
                % Calculate the NNRI score
                nnriScore = calculateNNRIForPoint(candidatePoint, points);
                
                % Calculate the VDARI score
                vdariScore = calculateVDARIScoreForPoint(candidatePoint, points, squareSideLength);
                
                % Combine the NNRI and VDARI scores using alpha and beta
                totalScore = alpha * nnriScore + beta * vdariScore;
                
                % If the total score of this candidate is better, select it
                if totalScore > bestScore
                    bestScore = totalScore;
                    bestPoint = candidatePoint;
                end
            end
            points = [points; bestPoint];
        end
    end
end

% Function to calculate NNRI score for a candidate point
function nnriScore = calculateNNRIForPoint(candidatePoint, points)
    distances = pdist2(candidatePoint, points);
    minDist = min(distances);
    nnriScore = minDist;  % NNRI score is based on the minimum distance to nearest neighbor
end

% Function to calculate VDARI score for a candidate point
function vdariScore = calculateVDARIScoreForPoint(candidatePoint, points, squareSideLength)
    tempPoints = [points; candidatePoint];  % Add candidate point to the list
    vdari = calculateVDARISquare(tempPoints, squareSideLength);
    if isnan(vdari)
        vdariScore = 0;  % If VDARI is NaN, assign a score of 0
    else
        vdariScore = vdari;  % Use the VDARI value as the score
    end
end

% Function to calculate NNRI
function nnri = calculateNNRI(points)
    distances = pdist2(points, points);
    distances(1:size(distances, 1)+1:end) = inf;
    nearestDistances = min(distances, [], 2);
    meanDist = mean(nearestDistances);
    stdDist = std(nearestDistances);
    if stdDist == 0
        nnri = NaN;
    else
        nnri = meanDist / stdDist;
    end
end

% Helper function to calculate VDARI
function vdari = calculateVDARISquare(points, squareSideLength)
    if size(points, 1) < 3
        vdari = NaN;
        return;
    end
    [~, v] = voronoin(points);  
    validAreas = [];
    for i = 1:length(v)
        if any(v{i} <= 0) || any(v{i} > size(points, 1))
            continue;
        end
        domainVertices = points(v{i}, :);
        if all(domainVertices(:, 1) >= 0 & domainVertices(:, 1) <= squareSideLength & ...
               domainVertices(:, 2) >= 0 & domainVertices(:, 2) <= squareSideLength)
            validAreas = [validAreas; polyarea(domainVertices(:, 1), domainVertices(:, 2))];
        end
    end
    if isempty(validAreas) || std(validAreas) == 0
        vdari = NaN;
    else
        meanArea = mean(validAreas);
        stdArea = std(validAreas);
        vdari = meanArea / stdArea;  
    end
end
