function [nearest_neighbor_distances, NNRI] = calculate_nearest_neighbors(x_coords, y_coords, n_random_points, visualize)
    % This function calculates the nearest neighbor distances for a randomly selected 
    % subset of points from the input coordinates. Optionally, it can visualize the points 
    % and their nearest neighbors.
    %
    % Inputs:
    %   x_coords - A vector of x coordinates for the points.
    %   y_coords - A vector of y coordinates for the points.
    %   n_random_points - The number of random points to select for the calculation.
    %   visualize - A boolean flag to determine whether to visualize the points and
    %               their nearest neighbors (default: false).
    %
    % Outputs:
    %   nearest_neighbor_distances - A vector containing the nearest neighbor distances
    %                                for each selected point.
    %   NNRI - The Nearest Neighbor Ratio Index, calculated as the mean nearest neighbor
    %          distance divided by the standard deviation of the nearest neighbor distances.

    if nargin < 4
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

    % Initialize an array to store the nearest neighbor distance for each point
    nearest_neighbor_distances = zeros(n_random_points, 1);

    % Calculate the nearest neighbor distance for each selected point
    for i = 1:n_random_points
        % Get the current point
        x_i = x_coords_selected(i);
        y_i = y_coords_selected(i);

        % Calculate distances from the current point to all other selected points
        distances = sqrt((x_coords_selected - x_i).^2 + (y_coords_selected - y_i).^2);

        % Ignore the distance from the point to itself (set it to Inf)
        distances(i) = Inf;

        % Find the minimum distance, which is the nearest neighbor distance
        nearest_neighbor_distances(i) = min(distances);
    end

    % Calculate the mean and standard deviation of the nearest neighbor distances
    mean_distance = mean(nearest_neighbor_distances);
    std_distance = std(nearest_neighbor_distances);

    % Calculate the Nearest Neighbor Ratio Index (NNRI)
    NNRI = mean_distance / std_distance;

    % Display the nearest neighbor distances
    disp('Nearest Neighbor Distances:');
    disp(nearest_neighbor_distances);

    % Display the NNRI
    disp('Nearest Neighbor Ratio Index (NNRI):');
    disp(NNRI);

    % Optional: Visualize the points and their nearest neighbors
    if visualize
        figure;
        scatter(x_coords_selected, y_coords_selected, 'filled');
        hold on;

        % Plot lines to nearest neighbors
        for i = 1:n_random_points
            % Find the index of the nearest neighbor
            [~, nearest_idx] = min(sqrt((x_coords_selected - x_coords_selected(i)).^2 + (y_coords_selected - y_coords_selected(i)).^2));

            % Draw a line between the current point and its nearest neighbor
            plot([x_coords_selected(i), x_coords_selected(nearest_idx)], [y_coords_selected(i), y_coords_selected(nearest_idx)], 'r-');
        end

        xlabel('X');
        ylabel('Y');
        title('Nearest Neighbor Visualization');
        axis equal;
        hold off;
    end
end