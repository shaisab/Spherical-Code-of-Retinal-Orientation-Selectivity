function [nnri_square,vdari_square,nnriDistribution, vdariDistribution,voronoi_areas, voronoi_vertices,transformedSmallShape] = calculateNNRIAndVDARIForSmallShape(smallShape, numPoints, requiredValidIterations)
    % This function calculates the NNRI and VDARI distributions for a small shape
    % placed randomly inside a square. The points are generated once and
    % then the shape is moved inside the square in each iteration.
    
    % Parameters
    squareSideLength = sqrt(16.5);  % Side length of the square (e.g., 4 mm)
    targetArea = 0.28;     % Target area of the small shape (in mm^2)

    % Preallocate arrays for results (overestimating required iterations)
    maxIterations = 20 * requiredValidIterations; % Allow for many invalid cases
    nnriResults = NaN(maxIterations, 1);         % NNRI results
    vdariResults = NaN(maxIterations, 1);        % VDARI results

    % Step 1: Generate points with the desired NNRI and VDARI
    points = generatePointsForMaximizingNNRI(squareSideLength, numPoints);
    nnri_square = calculateNNRISquare(points);
    vdari_square = calculateVDARISquare(points, squareSideLength);
    visualizeVoronoiSquare(points, squareSideLength);

    %voronoi based on square approximation
%     shape = [0 0; squareSideLength 0; squareSideLength squareSideLength; 0 squareSideLength]; % Square with side length 4
%     [vdari_square, nearestIndices] = calculateVDARI_new(points, shape, 100);
%     visualizeVoronoiFromVDARI(points, shape, 100, nearestIndices);
%     vdari_square = calculateVDARI(points, shape,100);

    % Step 2: Scale the shape to the target area
    scaledSmallShape = scaleShapeToArea(smallShape, targetArea);

    % Loop to achieve the required number of valid iterations
    validCount = 0;  % Count of valid NNRI and VDARI results
    totalIterations = 0;  % Total iterations attempted
    
    while validCount < requiredValidIterations
        % Step 3: Randomly place the small shape inside the square
        transformedSmallShape = randomPlacementInsideSquare(scaledSmallShape, squareSideLength);

        % Step 4: Find points inside the transformed shape
        pointsInsideShape = pointsInsideTransformedShape(transformedSmallShape, points);

        % Only proceed if there are points inside the shape
        if size(pointsInsideShape, 1) < 2  % Need at least 2 points to calculate NNRI or VDARI
            totalIterations = totalIterations + 1;
            continue;
        end

        % Step 5: Calculate NNRI for the points inside the shape
        nnri = calculateNNRI(pointsInsideShape);

        % Step 6: Calculate VDARI for the points inside the shape
        
        %voronoi based on square approximation
        %[vdari,nearestIndices] = calculateVDARI_new(pointsInsideShape, transformedSmallShape,100); 
        [voronoi_areas, voronoi_vertices, vdari] = calculate_vdari(pointsInsideShape, transformedSmallShape, false);
        
        % Store the results, only if NNRI and VDARI are valid
        if ~isnan(nnri) && ~isnan(vdari)  % Check if NNRI and VDARI are not NaN
            validCount = validCount + 1;
            nnriResults(validCount) = nnri;
            vdariResults(validCount) = vdari;
        end

        totalIterations = totalIterations + 1;
        
        % Safety check to prevent infinite loops
        if totalIterations > maxIterations
            warning('Exceeded max iterations without achieving the required valid iterations.');
            break;
        end
    end

    % Remove NaN values from results (cases where no valid points were inside the shape)
    nnriDistribution = nnriResults(~isnan(nnriResults));
    vdariDistribution = vdariResults(~isnan(vdariResults));
    
    visualizeVoronoiSquareSmall(points, squareSideLength, transformedSmallShape);
end

% Function to generate points with large NNRI
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
            for attempt = 1:2  % Maximum attempts to find a good location
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

% Function to scale the shape to the target area
function scaledShape = scaleShapeToArea(shape, targetArea)
    % Calculate the current area of the shape
    currentArea = polyarea(shape.Vertices(:, 1), shape.Vertices(:, 2));
    
    % Calculate the scaling factor to match the target area
    scaleFactor = sqrt(targetArea / currentArea);
    
    % Scale the shape by the scale factor
    scaledShape = scaleShape(shape, scaleFactor);
end

% Function to scale the shape to a specific factor
function scaledShape = scaleShape(shape, scaleFactor)
    scaledShape = shape;
    scaledShape.Vertices = scaledShape.Vertices * scaleFactor;
end

% Function to randomly place the shape inside the square
function transformedShape = randomPlacementInsideSquare(shape, squareSideLength)
    % Get the bounding box of the shape
    boundingBox = [min(shape.Vertices(:, 1)), min(shape.Vertices(:, 2)); 
                   max(shape.Vertices(:, 1)), max(shape.Vertices(:, 2))];
    
    % Calculate the width and height of the bounding box
    width = boundingBox(2, 1) - boundingBox(1, 1);
    height = boundingBox(2, 2) - boundingBox(1, 2);
    
    % Check if the shape fits in the square
    if width > squareSideLength || height > squareSideLength
        error('The shape does not fit inside the square!');
    end
    
    % Randomly offset the shape within the square
    offsetX = rand * (squareSideLength - width);
    offsetY = rand * (squareSideLength - height);
    
    % Apply the offset to the shape
    transformedShape = shape;
    transformedShape.Vertices = transformedShape.Vertices + [offsetX, offsetY];
end

% Function to check which points fall inside the transformed shape
function pointsInside = pointsInsideTransformedShape(shape, points)
    % Check if each point is inside the shape
    insideMask = inpolygon(points(:, 1), points(:, 2), shape.Vertices(:, 1), shape.Vertices(:, 2));
    
    % Return the points that are inside the shape
    pointsInside = points(insideMask, :);
end

% Function to calculate the NNRI (Mean nearest neighbor distance / Standard deviation)
function nnri = calculateNNRI(points)
    % Calculate pairwise distances between all points
    distances = pdist2(points, points);

    % Set the diagonal to infinity to exclude self-distances
    distances(1:size(distances, 1)+1:end) = inf;

    % Find the nearest neighbor distance for each point
    nearestDistances = min(distances, [], 2);

    % Check if there are valid distances
    if isempty(nearestDistances) || std(nearestDistances) == 0
        nnri = NaN;  % Return NaN if the NNRI cannot be computed
    else
        % Calculate the mean and standard deviation of the nearest neighbor distances
        meanDistance = mean(nearestDistances);
        stdDistance = std(nearestDistances);

        % Calculate the NNRI
        nnri = meanDistance / stdDistance;
    end
end

% Helper function to calculate NNRI (Nearest Neighbor Regularity Index)
% for intial square
function nnri = calculateNNRISquare(points)
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
% for intial square
function [vdari] = calculateVDARISquare(points, squareSideLength)
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

function [vdari] = calculateVDARI(points, shape,gridResolution)
    % Ensure there are at least three points to compute the Voronoi diagram
    if size(points, 1) < 3
        vdari = NaN;
        return;
    end

    % Validate and convert shape to numeric array if it is a polyshape
    if isa(shape, 'polyshape')
        shape = shape.Vertices; % Extract vertices from polyshape
    elseif ~isnumeric(shape) || size(shape, 2) ~= 2
        error('Shape must be an Mx2 numeric matrix or a polyshape object.');
    end

    % Precompute shape bounding box
    shapeBounds = [min(shape(:, 1)), max(shape(:, 1)), ...
                   min(shape(:, 2)), max(shape(:, 2))];
    xMin = shapeBounds(1);
    xMax = shapeBounds(2);
    yMin = shapeBounds(3);
    yMax = shapeBounds(4);

    % Generate grid
    if nargin < 3 || gridResolution <= 0
        gridResolution = 100; % Default: 100x100 grid
    end
    xGrid = linspace(xMin, xMax, gridResolution);
    yGrid = linspace(yMin, yMax, gridResolution);
    [X, Y] = meshgrid(xGrid, yGrid);

    % Precompute distances from grid points to input points
    gridPoints = [X(:), Y(:)];
    distances = pdist2(gridPoints, points);

    % Assign each grid point to its nearest point
    [~, nearestIndices] = min(distances, [], 2);

    % Approximate Voronoi cell areas by counting grid points per region
    regionCounts = accumarray(nearestIndices, 1, [size(points, 1), 1]);
    gridArea = (xMax - xMin) * (yMax - yMin) / (gridResolution^2);
    validAreas = regionCounts * gridArea;

    % Filter out zero areas (nonexistent Voronoi cells)
    validAreas = validAreas(validAreas > 0);

    % Calculate VDARI based on the areas of valid Voronoi cells
    if isempty(validAreas) || std(validAreas) == 0
        vdari = NaN;
    else
        meanArea = mean(validAreas);
        stdArea = std(validAreas);
        vdari = meanArea / stdArea;
    end
end