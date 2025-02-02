function visualizeVoronoiFromVDARI(points, shape, gridResolution, nearestIndices)
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

    % Plot the Voronoi domains
    figure;
    imagesc(xGrid, yGrid, voronoiMap); % Display the regions
    axis xy; % Correct axis orientation
    colormap(jet(max(nearestIndices))); % Unique color for each region
    colorbar; % Add a color bar

    % Overlay input points
    hold on;
    scatter(points(:, 1), points(:, 2), 50, 'k', 'filled', 'MarkerEdgeColor', 'w');

    % Overlay the shape boundary
    plot([shape(:, 1); shape(1, 1)], [shape(:, 2); shape(1, 2)], 'k-', 'LineWidth', 2);

    xlim([0 sqrt(16.5)])
    ylim([0 sqrt(16.5)])
    
    title('Voronoi Domains Visualization');
    xlabel('X');
    ylabel('Y');
    hold off;
end
