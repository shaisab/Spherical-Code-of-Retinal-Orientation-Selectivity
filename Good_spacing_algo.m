 close all
 clear all

 numPoints=1765;
 squareSideLength=sqrt(16.5);
 points = generatePointsForMaximizingNNRI(squareSideLength, numPoints);
 nnri = calculateNNRI(points);
 vdari = calculateVDARISquare(points, squareSideLength);
 
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

% Helper function to calculate VDARI (Voronoi Domain Area Regularity Index)
function vdari = calculateVDARISquare(points, squareSideLength)
    % Ensure there are at least three points to compute the Voronoi diagram
    if size(points, 1) < 3
        vdari = NaN;
        return;
    end
    
    % Compute the Voronoi diagram for the points
    [~, v] = voronoin(points);  
    
    validAreas = [];
    for i = 1:length(v)
        % Avoid invalid Voronoi cells (those outside the square)
        if any(v{i} <= 0) || any(v{i} > size(points, 1))
            continue;
        end
        
        % Get the vertices of the Voronoi cell
        domainVertices = points(v{i}, :);
        
        % Check if the Voronoi cell lies inside the square boundary
        if all(domainVertices(:, 1) >= 0 & domainVertices(:, 1) <= squareSideLength & ...
               domainVertices(:, 2) >= 0 & domainVertices(:, 2) <= squareSideLength)
            % Calculate the area of the Voronoi cell
            validAreas = [validAreas; polyarea(domainVertices(:, 1), domainVertices(:, 2))];
        end
    end
    
    if isempty(validAreas) || std(validAreas) == 0
        vdari = NaN;
    else
        % VDARI is the ratio of the mean area to the standard deviation of the areas
        meanArea = mean(validAreas);
        stdArea = std(validAreas);
        vdari = meanArea / stdArea;  
    end
end