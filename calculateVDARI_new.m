function [vdari, nearestIndices] = calculateVDARI_new(points, shape, gridResolution)
    % Ensure there are at least three points to compute the Voronoi diagram
    if size(points, 1) < 3
        vdari = NaN;
        nearestIndices = [];
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
