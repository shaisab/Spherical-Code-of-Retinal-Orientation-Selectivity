function nnriDistribution = calculateNNRIForSmallShapeValid(smallShape, numPoints, requiredValidIterations)
    % This function calculates the NNRI distribution for a small shape
    % placed randomly inside a square. The points are generated once and
    % then the shape is moved inside the square in each iteration.
    
    % Parameters
    squareSideLength = 4;  % Side length of the square (e.g., 4 mm)
    targetNNRI = 4.5;      % Target NNRI for the points (mean distance / std deviation)
    targetArea = 0.28;     % Target area of the small shape (in mm^2)

    % Preallocate an array for NNRI results (overestimating required iterations)
    maxIterations = 10 * requiredValidIterations; % Allow for many invalid cases
    nnriResults = NaN(maxIterations, 1);  % Initialize with NaN to track invalid cases

    % Step 1: Generate points with the desired NNRI of 4.5
    points = generatePointsWithTargetNNRI(squareSideLength, numPoints, targetNNRI);

    % Step 2: Scale the shape to the target area
    scaledSmallShape = scaleShapeToArea(smallShape, targetArea);

    % Loop to achieve the required number of valid iterations
    validCount = 0;  % Count of valid NNRI results
    totalIterations = 0;  % Total iterations attempted
    
    while validCount < requiredValidIterations
        % Step 3: Randomly place the small shape inside the square
        transformedSmallShape = randomPlacementInsideSquare(scaledSmallShape, squareSideLength);

        % Step 4: Find points inside the transformed shape
        pointsInsideShape = pointsInsideTransformedShape(transformedSmallShape, points);

        % Only proceed if there are points inside the shape
        if size(pointsInsideShape, 1) < 2  % Need at least 2 points to calculate NNRI
            totalIterations = totalIterations + 1;
            continue;
        end

        % Step 5: Calculate NNRI for the points inside the shape
        nnri = calculateNNRI(pointsInsideShape);

        % Store the result, only if NNRI is valid
        if ~isnan(nnri)  % Check if NNRI is not NaN
            validCount = validCount + 1;
            nnriResults(validCount) = nnri;
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
end

% Function to generate points with a specific NNRI
function points = generatePointsWithTargetNNRI(squareSideLength, numPoints, targetNNRI)
    % Estimate the distance between points based on the target NNRI
    approxMeanNNDist = targetNNRI;  % Target NNRI is the ratio of mean distance to std deviation
    
    % Calculate the average distance between points required to achieve the target NNRI
    targetPointSpacing = approxMeanNNDist * sqrt(2);  % Approximate distance for uniform spacing
    
    % Calculate the total area available in the square
    totalArea = squareSideLength^2;
    
    % Estimate the number of points required based on the desired spacing
    estimatedPointCount = round(totalArea / (targetPointSpacing^2));
    
    % If the required points are greater than the desired number, limit the number of points
    if estimatedPointCount > numPoints
        numPoints = estimatedPointCount;
    end

    % Generate points in a regular grid that approximately achieves the desired NNRI
    points = generatePointsInGrid(squareSideLength, numPoints, targetPointSpacing);
end

% Function to generate points in a regular grid inside the square
function points = generatePointsInGrid(squareSideLength, numPoints, targetSpacing)
    % Calculate the grid size based on the target spacing
    gridSize = ceil(sqrt(numPoints));  % Create a square grid
    
    % Create points in a regular grid pattern
    [xGrid, yGrid] = meshgrid(linspace(0, squareSideLength, gridSize), linspace(0, squareSideLength, gridSize));
    
    % Flatten the grid and select the first 'numPoints' points
    points = [xGrid(:), yGrid(:)];
    points = points(1:numPoints, :);
    
    % Adjust spacing if necessary
    distances = pdist2(points, points);
    minDistance = min(distances(distances > 0));
    if minDistance < targetSpacing
        % Spread points apart if necessary (this is a simple example, more complex logic may be needed)
        points = points + rand(size(points)) * (targetSpacing - minDistance);
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
