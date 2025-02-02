function [all_voronoi_areas, all_voronoi_vertices, vdari_values] = calculate_voronoi_distribution(x_coords, y_coords, n_random_points, bounding_polygon, num_iterations, visualize)
    % This function calculates the Voronoi diagram and the areas of the Voronoi cells 
    % for a randomly selected subset of points from the input coordinates within a 
    % specified bounding polygon, repeating the process num_iterations times.
    % It also calculates the Voronoi Domain Area Regularity Index (VDARI) for each iteration.
    %
    % Inputs:
    %   x_coords - A vector of x coordinates for the points.
    %   y_coords - A vector of y coordinates for the points.
    %   n_random_points - The number of random points to select for each calculation.
    %   bounding_polygon - A polyshape object or a set of coordinates defining the bounding area.
    %   num_iterations - The number of iterations to perform (e.g., 1000).
    %   visualize - A boolean flag to determine whether to visualize the Voronoi diagram (default: false).
    %
    % Outputs:
    %   all_voronoi_areas - A matrix containing the areas of the Voronoi cells for each iteration and each selected point.
    %   all_voronoi_vertices - A cell array containing the vertices of the Voronoi cells for each iteration and each selected point.
    %   vdari_values - A vector containing the Voronoi Domain Area Regularity Index (VDARI) for each iteration.

    if nargin < 6
        visualize = false;  % Set default value for visualization flag
    end

    % Ensure that the number of random points doesn't exceed the available points
    n_total_points = length(x_coords);
    if n_random_points > n_total_points
        error('Number of random points exceeds the total number of available points.');
    end

    % Initialize outputs
    all_voronoi_areas = zeros(num_iterations, n_random_points);
    all_voronoi_vertices = cell(num_iterations, n_random_points);
    vdari_values = zeros(num_iterations, 1);

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

    for iter = 1:num_iterations
        % Randomly select n_random_points indices from the available points
        selected_indices = randperm(n_total_points, n_random_points);

        % Extract the randomly selected coordinates
        x_coords_selected = x_coords(selected_indices);
        y_coords_selected = y_coords(selected_indices);

        % Compute the Voronoi diagram for the selected points
        [v, c] = voronoin([x_coords_selected, y_coords_selected]);

        % Calculate the areas of the Voronoi cells and store the vertices
        cell_areas = zeros(n_random_points, 1);
        for i = 1:n_random_points
            % Extract the vertices of the Voronoi cell for point i
            cell_indices = c{i};
            if isempty(cell_indices)
                all_voronoi_vertices{iter, i} = [];
                all_voronoi_areas(iter, i) = 0;
                cell_areas(i) = 0;
                continue;
            end

            cell_vertices = v(cell_indices, :);

            % Remove infinite vertices
            finite_vertices = cell_vertices(~any(isinf(cell_vertices), 2), :);
            if size(finite_vertices, 1) < 3
                all_voronoi_vertices{iter, i} = [];
                all_voronoi_areas(iter, i) = 0;
                cell_areas(i) = 0;
                continue;
            end

            % Create a polyshape for the Voronoi cell
            cell_polygon = polyshape(finite_vertices);

            % Clip the Voronoi cell to the bounding polygon
            try
                clipped_polygon = intersect(cell_polygon, bounding_polygon);
                if ~isempty(clipped_polygon)
                    all_voronoi_vertices{iter, i} = clipped_polygon.Vertices;
                    area_value = area(clipped_polygon);
                    all_voronoi_areas(iter, i) = area_value;
                    cell_areas(i) = area_value;
                else
                    all_voronoi_vertices{iter, i} = [];
                    all_voronoi_areas(iter, i) = 0;
                    cell_areas(i) = 0;
                end
            catch
                % Handle any errors that occur during clipping
                all_voronoi_vertices{iter, i} = [];
                all_voronoi_areas(iter, i) = 0;
                cell_areas(i) = 0;
            end
        end

        % Calculate VDARI for this iteration
        non_zero_areas = cell_areas(cell_areas > 0); % Exclude zero areas
        if numel(non_zero_areas) > 1
            mean_area = mean(non_zero_areas);
            std_area = std(non_zero_areas);
            vdari_values(iter) = std_area / mean_area; % Coefficient of Variation
        else
            vdari_values(iter) = NaN; % Assign NaN if there are not enough non-zero areas
        end
    end

    % Optional: Visualize the Voronoi diagram for the last iteration
    if visualize
        figure;
        hold on;

        % Plot the bounding polygon using the polyshape plot method
        plot(bounding_polygon, 'FaceColor', 'none', 'EdgeColor', 'k');

        % Plot the Voronoi cells for the last iteration
        for i = 1:n_random_points
            if ~isempty(all_voronoi_vertices{num_iterations, i})
                plot(polyshape(all_voronoi_vertices{num_iterations, i}), 'FaceColor', 'none', 'EdgeColor', 'r');
            end
        end

        % Plot the original points
        scatter(x_coords, y_coords, 'filled');

        xlabel('X');
        ylabel('Y');
        title('Voronoi Diagram (Last Iteration)');
        axis equal;
        hold off;
    end
end
