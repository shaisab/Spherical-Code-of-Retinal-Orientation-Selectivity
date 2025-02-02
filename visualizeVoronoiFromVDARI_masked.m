function visualizeVoronoiFromVDARI_masked(points, shape, gridResolution, nearestIndices)
    % Validate input
    if nargin < 4 || isempty(nearestIndices)
        error('nearestIndices must be provided to visualize the Voronoi domains.');
    end
    if size(points, 1) < 3
        error('At least three points are required to compute Voronoi domains.');
    end
    if isa(shape, 'polyshape')
        shape = shape.Vertices; % Extract vertices from polyshape
    elseif ~isnumeric(shape) || size(shape, 2) ~= 2
        error('Shape must be an Mx2 numeric matrix or a polyshape object.');
    end

    % Shape bounds
    shapeBounds = [min(shape(:, 1)), max(shape(:, 1)), ...
                   min(shape(:, 2)), max(shape(:, 2))];
    xMin = shapeBounds(1);
    xMax = shapeBounds(2);
    yMin = shapeBounds(3);
    yMax = shapeBounds(4);

    % Generate grid
    xGrid = linspace(xMin, xMax, gridResolution);
    yGrid = linspace(yMin, yMax, gridResolution);
    [X, Y] = meshgrid(xGrid, yGrid);

    % Reshape nearestIndices to grid shape
    voronoiMap = reshape(nearestIndices, gridResolution, gridResolution);

    % Create a mask for points outside the polygon
    inPolygon = inpolygon(X, Y, shape(:, 1), shape(:, 2));
    voronoiMap(~inPolygon) = NaN; % Mask out values outside the polygon

    % Plot the Voronoi domains
    figure;
    h = imagesc(xGrid, yGrid, voronoiMap); % Display the regions
    axis xy; % Correct axis orientation
    colormap(jet(max(nearestIndices))); % Use a colormap for distinct colors
    colorbar; % Add a color bar
    
    % Set NaN values to white color (background area)
    set(h, 'AlphaData', ~isnan(voronoiMap)); % Only show non-NaN regions
    colormap(gca, [1 1 1; jet(max(nearestIndices))]); % Set first color as white

    % Overlay input points
    hold on;
    scatter(points(:, 1), points(:, 2), 50, 'k', 'filled', 'MarkerEdgeColor', 'w');

    % Overlay the shape boundary
    plot([shape(:, 1); shape(1, 1)], [shape(:, 2); shape(1, 2)], 'k-', 'LineWidth', 2);

    % Limit the axes to match the bounding box of the shape (polygon or grid)
    axis([xMin, xMax, yMin, yMax]);

    % Customize plot
    title('Voronoi Domains Visualization (Masked by Polygon)');
    xlabel('X');
    ylabel('Y');
    axis equal;
    hold off;
end
