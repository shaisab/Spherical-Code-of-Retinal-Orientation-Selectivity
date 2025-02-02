 close all
 clear all

numPoints = 1765; 
squareSideLength = sqrt(16.5);
points = generatePointsForMaximizingNNRI(squareSideLength, numPoints);

% Remove points with Voronoi domain area below threshold
% areaThreshold = 1.4;
% points = removePointsWithSmallVoronoiAreas(points, squareSideLength, areaThreshold);

% Recalculate metrics after removal
nnri = calculateNNRI(points);
[~,vdari] = calculateVDARISquare(points, squareSideLength);

function points = removePointsWithSmallVoronoiAreas(points, squareSideLength, areaThreshold)
    
    % Calculate Voronoi areas
    [validAreas,vdari] = calculateVDARISquare(points, squareSideLength);
    
    % Identify points to keep
    pointsToKeep = validAreas <= areaThreshold;
    
    % Remove points with small Voronoi areas
    points = points(pointsToKeep, :);
end
 
function points = generatePointsForMaximizingNNRI(squareSideLength, numPoints)
    % Generate points inside a square, trying to maximize NNRI (Mean/Std dev)
    
    % Step 1: Start with an empty point array
    points = [];
    
    % Step 2: Iteratively place points in the square
    for i = 1:numPoints
        if i == 1
            % Place the first point randomly in the square
            points = [rand(1, 1) * squareSideLength, rand(1, 1) * squareSideLength];
        else
            % For subsequent points, maximize the minimum distance
            bestPoint = [];
            bestMinDist = -inf;
            
            % Try placing a point at several locations inside the square
            for attempt = 1:5  % Maximum attempts to find a good location
                candidatePoint = [rand * squareSideLength, rand * squareSideLength];
                
                % Calculate the distances to all other points already placed
                distances = pdist2(candidatePoint, points);
                
                % Find the minimum distance from this candidate point to any other point
                minDist = min(distances);
                
                % Only accept the point if the minimum distance is larger than the best found
                if minDist > bestMinDist
                    bestMinDist = minDist;
                    bestPoint = candidatePoint;
                end
            end
            
            % Add the best point found in this iteration
            points = [points; bestPoint];
        end
    end
end

function nnri = calculateNNRI(points)
    % Calculate the Nearest Neighbor Regularity Index (NNRI)
    
    % Calculate pairwise distances between all points
    distances = pdist2(points, points);
    
    % Set diagonal to inf to exclude self-distances
    distances(1:size(distances, 1)+1:end) = inf;
    
    % Find the nearest neighbor distance for each point
    nearestDistances = min(distances, [], 2);
    
    % Calculate NNRI: mean distance / standard deviation of distances
    meanDist = mean(nearestDistances);
    stdDist = std(nearestDistances);
    
    if stdDist == 0
        nnri = NaN;  % If standard deviation is zero, return NaN
    else
        nnri = meanDist / stdDist;
    end
end

function [validAreas, vdari] = calculateVDARISquare(points, squareSideLength)
    % Ensure there are at least three points to compute the Voronoi diagram
    if size(points, 1) < 3
        vdari = NaN;
        validAreas = [];
        return;
    end
    
    % Compute the Voronoi diagram for the points
    [v, voronoiCells] = voronoin(points);
    
    validAreas = [];
    for i = 1:length(voronoiCells)
        cellVertices = voronoiCells{i};
        
        % Skip invalid or infinite Voronoi cells
        if any(cellVertices == 1) || isempty(cellVertices)
            continue;
        end
        
        % Get the coordinates of the cell vertices
        vertexCoords = v(cellVertices, :);
        
        % Clip the vertices to the boundaries of the square
        vertexCoords(:, 1) = max(0, min(squareSideLength, vertexCoords(:, 1)));
        vertexCoords(:, 2) = max(0, min(squareSideLength, vertexCoords(:, 2)));
        
        % Check if all vertices are within bounds after clipping
        if all(vertexCoords(:, 1) >= 0 & vertexCoords(:, 1) <= squareSideLength & ...
               vertexCoords(:, 2) >= 0 & vertexCoords(:, 2) <= squareSideLength)
            % Calculate the area of the Voronoi cell
            validAreas = [validAreas; polyarea(vertexCoords(:, 1), vertexCoords(:, 2))];
        end
    end
    
    % Calculate VDARI based on the areas of valid Voronoi cells
    if isempty(validAreas) || std(validAreas) == 0
        vdari = NaN;
    else
        meanArea = mean(validAreas);
        stdArea = std(validAreas);
        vdari = meanArea / stdArea;
    end
end
