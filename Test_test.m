close all
clear all

numPoints = 1765;
squareSideLength = sqrt(16.5);
points = generatePointsForMaximizingNNRI(squareSideLength, numPoints);
nnri = calculateNNRI(points);
vdari = calculateVDARISquare(points, squareSideLength);

% Display NNRI and VDARI
disp(['NNRI: ', num2str(nnri)]);
disp(['VDARI: ', num2str(vdari)]);

% Plot the points
figure;
scatter(points(:, 1), points(:, 2), 'filled');
title('Generated Points with Maximized NNRI');
xlabel('X');
ylabel('Y');
axis equal;

% Function to generate points for maximizing NNRI
function points = generatePointsForMaximizingNNRI(squareSideLength, numPoints)
    points = [];
    for i = 1:numPoints
        if i == 1
            points = [rand(1, 1) * squareSideLength, rand(1, 1) * squareSideLength];
        else
            bestPoint = [];
            bestMinDist = -inf;
            for attempt = 1:5  % Max attempts per point
                candidatePoint = [rand * squareSideLength, rand * squareSideLength];
                distances = pdist2(candidatePoint, points);
                minDist = min(distances);
                if minDist > bestMinDist
                    bestMinDist = minDist;
                    bestPoint = candidatePoint;
                end
            end
            points = [points; bestPoint];
        end
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
