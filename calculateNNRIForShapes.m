function nnriDistribution = calculateNNRIForShapes(largeShape, smallShape, largeArea, smallArea, numPoints, nnriTarget, numIterations)
    % Scaling factor to adjust small shape area to match target size
    scalingFactor = sqrt(largeArea / smallArea);
    
    % Scale the small shape to the desired area
    smallShapeScaled = scale(smallShape, scalingFactor);
    
    % Precompute the vertices for fast processing in 'inpolygon'
    largeShapeVertices = largeShape.Vertices;
    smallShapeScaledVertices = smallShapeScaled.Vertices;

    % Initialize the array for NNRI values
    nnriValues = NaN(numIterations, 1);

    % Step 1: Generate random points inside the large shape with the target NNRI
    points = generatePointsWithNNRI(largeShape, numPoints, nnriTarget);

    % Step 2 - 4: Run the iterations
    for i = 1:numIterations
        % Step 2: Randomly place the smaller shape inside the large shape
        transformedSmallShape = randomPlacement(largeShape, smallShapeScaled);

        % Step 3: Extract points inside the transformed small shape
        pointsInsideSmallShape = extractPointsInsideShape(points, transformedSmallShape);

        % Step 4: Calculate the NNRI for points inside the small shape
        nnriValues(i) = calculateNNRI(pointsInsideSmallShape);
    end

    % Return the distribution of NNRI values
    nnriDistribution = nnriValues;
end

function transformedSmallShape = randomPlacement(largeShape, smallShapeScaled)
    % Get the bounding box of both shapes
    largeBoundingBox = [min(largeShape.Vertices, [], 1); max(largeShape.Vertices, [], 1)];
    smallBoundingBox = [min(smallShapeScaled.Vertices, [], 1); max(smallShapeScaled.Vertices, [], 1)];

    % Randomly place the small shape inside the large shape
    while true
        offsetX = randi([0, round(largeBoundingBox(2, 1) - smallBoundingBox(2, 1))]);
        offsetY = randi([0, round(largeBoundingBox(2, 2) - smallBoundingBox(2, 2))]);

        % Apply the offset to the small shape
        transformedSmallShape = smallShapeScaled;
        transformedSmallShape.Vertices = transformedSmallShape.Vertices + [offsetX, offsetY];

        % Check if the small shape is completely inside the large shape
        if all(inpolygon(transformedSmallShape.Vertices(:, 1), transformedSmallShape.Vertices(:, 2), largeShape.Vertices(:, 1), largeShape.Vertices(:, 2)))
            break;
        end
    end
end

function pointsInsideSmallShape = extractPointsInsideShape(points, transformedSmallShape)
    % Extract the X and Y coordinates of the points
    x = points(:, 1);
    y = points(:, 2);

    % Check which points are inside the transformed small shape
    [in, ~] = inpolygon(x, y, transformedSmallShape.Vertices(:, 1), transformedSmallShape.Vertices(:, 2));

    % Return only the points inside the small shape
    pointsInsideSmallShape = points(in, :);
end

function nnri = calculateNNRI(points)
    % If there are fewer than 2 points, return NaN
    if size(points, 1) < 2
        nnri = NaN;
        return;
    end

    % Calculate Euclidean distances between points
    distMatrix = pdist2(points, points);
    distMatrix(distMatrix == 0) = Inf;  % Ignore self-distances
    nearestDistances = min(distMatrix, [], 2);  % Nearest neighbor distances

    % Calculate mean and standard deviation of nearest neighbor distances
    meanDistance = mean(nearestDistances);
    stdDistance = std(nearestDistances);

    % Return NNRI (mean / std)
    if stdDistance == 0
        nnri = NaN;
    else
        nnri = meanDistance / stdDistance;
    end
end

function points = generatePointsWithNNRI(shape, numPoints, nnriTarget)
    % Generate points inside the shape with the desired NNRI
    points = [];
    while true
        % Generate random points inside the shape
        points = generateRandomPointsInside(shape, numPoints);
        
        % Calculate the NNRI of the generated points
        nnri = calculateNNRI(points);
        
        % Check if the NNRI is close to the target
        if abs(nnri - nnriTarget) < 0.1
            break;
        end
    end
end

function points = generateRandomPointsInside(shape, numPoints)
    % Get the bounding box of the shape
    minX = min(shape.Vertices(:, 1));
    maxX = max(shape.Vertices(:, 1));
    minY = min(shape.Vertices(:, 2));
    maxY = max(shape.Vertices(:, 2));

    % Generate random points inside the bounding box
    while true
        x = (maxX - minX) * rand(numPoints, 1) + minX;
        y = (maxY - minY) * rand(numPoints, 1) + minY;

        % Check if points are inside the shape using inpolygon
        [in, ~] = inpolygon(x, y, shape.Vertices(:, 1), shape.Vertices(:, 2));

        if sum(in) == numPoints
            break;
        end
    end

    % Return the points that are inside the shape
    points = [x(in), y(in)];
end
