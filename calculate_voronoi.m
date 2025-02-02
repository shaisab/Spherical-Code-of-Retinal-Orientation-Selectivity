function [voronoi_areas, voronoi_vertices] = calculate_voronoi(x_coords, y_coords, n_random_points, bounding_polygon, visualize)
    % This function calculates the Voronoi diagram and the areas of the Voronoi cells 
    % for a randomly selected subset of points from the input coordinates within a 
    % specified bounding polygon.
    %
    % Inputs:
    %   x_coords - A vector of x coordinates for the points.
    %   y_coords - A vector of y coordinates for the points.
    %   n_random_points - The number of random points to select for the calculation.
    %   bounding_polygon - A polyshape object or a set of coordinates defining the bounding area.
    %   visualize - A boolean flag to determine whether to visualize the Voronoi diagram (default: false).
    %
    % Outputs:
    %   voronoi_areas - A vector containing the areas of the Voronoi cells for each selected point.
    %   voronoi_vertices - A cell array containing the vertices of the Voronoi cells for each selected point.

    if nargin < 5
        visualize = false;  % Set default value for visualization flag
    end

    % Ensure that the number of random points doesn't exceed the available points
    n_total_points = length(x_coords);
    if n_random_points > n_total_points
        error('Number of random points exceeds the total number of available points.');
    end

    % Randomly select n_random_points indices from the available points
    selected_indices = randperm(n_total_points, n_random_points);

    % Extract the randomly selected coordinates
    x_coords_selected = x_coords(selected_indices);
    y_coords_selected = y_coords(selected_indices);

    % Compute the Voronoi diagram for the selected points
    [v, c] = voronoin([x_coords_selected, y_coords_selected]);

    % Initialize outputs
    voronoi_areas = zeros(n_random_points, 1);
    voronoi_vertices = cell(n_random_points, 1);

    % Check if bounding_polygon is a polyshape object or a set of coordinates
    if isa(bounding_polygon, 'polyshape')
        bounding_vertices = bounding_polygon.Vertices;
    elseif ismatrix(bounding_polygon) && size(bounding_polygon, 2) == 2
        bounding_vertices = bounding_polygon;
    else
        error('Invalid bounding_polygon format. It should be a polyshape object or a set of coordinates.');
    end

    % Calculate the areas of the Voronoi cells and store the vertices
    for i = 1:n_random_points
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

        % Clip the Voronoi cell to the bounding polygon using inpolygon
        if ~isempty(finite_vertices)
            % Check which vertices are inside the bounding polygon
            [in, ~] = inpolygon(finite_vertices(:, 1), finite_vertices(:, 2), bounding_vertices(:, 1), bounding_vertices(:, 2));
            clipped_vertices = finite_vertices(in, :);

            % Create a polyshape for the clipped Voronoi cell
            if size(clipped_vertices, 1) >= 3
                try
                    clipped_polygon = polyshape(clipped_vertices);
                    voronoi_vertices{i} = clipped_polygon.Vertices;
                    voronoi_areas(i) = area(clipped_polygon);
                catch
                    % If creation of polyshape fails, set empty results
                    voronoi_vertices{i} = [];
                    voronoi_areas(i) = 0;
                end
            else
                voronoi_vertices{i} = [];
                voronoi_areas(i) = 0;
            end
        else
            voronoi_vertices{i} = [];
            voronoi_areas(i) = 0;
        end
    end

    % Optional: Visualize the Voronoi diagram
    if visualize
        figure;
        hold on;

        % Plot the bounding polygon using the polyshape plot method
        if isa(bounding_polygon, 'polyshape')
            plot(bounding_polygon, 'FaceColor', 'none', 'EdgeColor', 'k');
        else
            plot(polyshape(bounding_polygon), 'FaceColor', 'none', 'EdgeColor', 'k');
        end

        % Plot the Voronoi cells
        for i = 1:n_random_points
            if ~isempty(voronoi_vertices{i})
                plot(polyshape(voronoi_vertices{i}), 'FaceColor', 'none', 'EdgeColor', 'r');
            end
        end

        % Plot the original points
        scatter(x_coords_selected, y_coords_selected, 'filled');

        xlabel('X');
        ylabel('Y');
        title('Voronoi Diagram');
        axis equal;
        hold off;
    end
end

