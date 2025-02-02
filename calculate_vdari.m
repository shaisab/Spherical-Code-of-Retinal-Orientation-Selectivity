function [voronoi_areas, voronoi_vertices, vdari] = calculate_vdari(points, bounding_polygon, visualize)
    % This function calculates the Voronoi diagram and the areas of the Voronoi cells 
    % for a randomly selected subset of points from the input coordinates within a 
    % specified bounding polygon, and calculates the Voronoi Domain Area Regularity Index (VDARI).
    %
    % Inputs:
    %   points - A matrix of size (n_points x 2) containing the coordinates of points.
    %   bounding_polygon - A polyshape object or a set of coordinates defining the bounding area.
    %   visualize - A boolean flag to determine whether to visualize the Voronoi diagram (default: false).
    %
    % Outputs:
    %   voronoi_areas - A vector containing the areas of the Voronoi cells for each selected point.
    %   voronoi_vertices - A cell array containing the vertices of the Voronoi cells for each selected point.
    %   vdari - The Voronoi Domain Area Regularity Index, calculated as the coefficient of variation of the Voronoi areas.

    if nargin < 3
        visualize = false;  % Set default value for visualization flag
    end
    
    n_points = size(points, 1);
    
    % Try to compute the Voronoi diagram for the points
    try
        [v, c] = voronoin(points);
    catch
        % If there's an error (e.g., not enough unique points), return NaN for vdari and exit
        voronoi_areas = NaN(n_points, 1);
        voronoi_vertices = cell(n_points, 1);
        vdari = NaN;
        return;
    end

    % Check if bounding_polygon is a polyshape object or a set of coordinates
    if isa(bounding_polygon, 'polyshape')
        bounding_vertices = bounding_polygon.Vertices;
    elseif ismatrix(bounding_polygon) && size(bounding_polygon, 2) == 2
        bounding_vertices = bounding_polygon;
    else
        error('Invalid bounding_polygon format. It should be a polyshape object or a set of coordinates.');
    end

    % Convert bounding_vertices to a polyshape
    bounding_polygon = polyshape(bounding_vertices);

    % Initialize outputs
    voronoi_areas = zeros(n_points, 1);
    voronoi_vertices = cell(n_points, 1);

    % Calculate the areas of the Voronoi cells and store the vertices
    for i = 1:n_points
        % Extract the vertices of the Voronoi cell for point i
        cell_indices = c{i};
        if isempty(cell_indices)
            voronoi_vertices{i} = [];
            voronoi_areas(i) = 0;
            continue;
        end

        cell_vertices = v(cell_indices, :);

        % Remove infinite vertices
        finite_vertices = cell_vertices(~any(isinf(cell_vertices), 2), :);
        if size(finite_vertices, 1) < 3
            voronoi_vertices{i} = [];
            voronoi_areas(i) = 0;
            continue;
        end

        % Create a polyshape for the Voronoi cell
        cell_polygon = polyshape(finite_vertices);

        % Clip the Voronoi cell to the bounding polygon
        try
            clipped_polygon = intersect(cell_polygon, bounding_polygon);
            if ~isempty(clipped_polygon)
                voronoi_vertices{i} = clipped_polygon.Vertices;
                voronoi_areas(i) = area(clipped_polygon);
            else
                voronoi_vertices{i} = [];
                voronoi_areas(i) = 0;
            end
        catch
            % Handle any errors that occur during clipping
            voronoi_vertices{i} = [];
            voronoi_areas(i) = 0;
        end
    end

    % Check if the number of Voronoi cells (areas) is equal to the number of points
    if sum(voronoi_areas > 0) ~= n_points
        vdari = NaN;  % Return NaN if the number of valid areas does not match the number of points
        return;
    end

    % Calculate VDARI
    non_zero_areas = voronoi_areas(voronoi_areas > 0); % Exclude zero areas
    if numel(non_zero_areas) > 1
        mean_area = mean(non_zero_areas);
        std_area = std(non_zero_areas);
        vdari = mean_area / std_area; 
    else
        vdari = NaN; % Assign NaN if there are not enough non-zero areas
    end

    % Optional: Visualize the Voronoi diagram
    if visualize
        figure;
        hold on;

        % Plot the bounding polygon using the polyshape plot method
        plot(bounding_polygon, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 2);

        % Plot the Voronoi cells
        for i = 1:n_points
            if ~isempty(voronoi_vertices{i})
                plot(polyshape(voronoi_vertices{i}), 'FaceColor', 'none', 'EdgeColor', 'r', 'LineWidth', 1);
            end
        end

        % Plot the original points
        scatter(points(:, 1), points(:, 2), 50, 'b', 'filled');

        % Title and axis settings
        xlabel('X');
        ylabel('Y');
        title('Voronoi Diagram with Clipped Cells');
        axis equal;
        hold off;
    end
end
