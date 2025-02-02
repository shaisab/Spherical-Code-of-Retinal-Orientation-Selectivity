function [NNRI_distribution, nearest_neighbor_distances_distribution] = calculate_nnri_distribution(x_coords, y_coords, num_random_points, num_iterations, visualize)
    % This function calculates the Nearest Neighbor Ratio Index (NNRI) 
    % for a specified number of random points out of the total number of points.
    % The process is repeated num_iterations times to produce a distribution of NNRI values.
    %
    % Inputs:
    %   x_coords - A vector of x coordinates for the points.
    %   y_coords - A vector of y coordinates for the points.
    %   num_random_points - The number of random points to select in each iteration.
    %   num_iterations - The number of iterations to perform.
    %   visualize - A boolean flag to determine whether to visualize the points and
    %               their nearest neighbors in the final iteration (default: false).
    %
    % Outputs:
    %   NNRI_distribution - A vector containing the NNRI values for each iteration.
    %   nearest_neighbor_distances_distribution - A matrix containing the nearest 
    %                                             neighbor distances for each iteration.
    
    if nargin < 5
        visualize = false;  % Set default value for visualization flag
    end

    % Initialize an array to store NNRI for each iteration
    NNRI_distribution = zeros(num_iterations, 1);

    % Initialize a matrix to store nearest neighbor distances for each iteration
    nearest_neighbor_distances_distribution = zeros(num_random_points, num_iterations);

    % Perform the process num_iterations times
    for iter = 1:num_iterations
        % Randomly select num_random_points indices from the total number of points
        random_indices = randperm(length(x_coords), num_random_points);
        
        % Get the random subset of points
        x_subset = x_coords(random_indices);
        y_subset = y_coords(random_indices);
        
        % Calculate the nearest neighbor distances for the random subset
        nearest_neighbor_distances = zeros(num_random_points, 1);
        
        for i = 1:num_random_points
            % Get the current point
            x_i = x_subset(i);
            y_i = y_subset(i);

            % Calculate distances from the current point to all other points
            distances = sqrt((x_subset - x_i).^2 + (y_subset - y_i).^2);

            % Ignore the distance from the point to itself (set it to Inf)
            distances(i) = Inf;

            % Find the minimum distance, which is the nearest neighbor distance
            nearest_neighbor_distances(i) = min(distances);
        end
        
        % Store the nearest neighbor distances for this iteration
        nearest_neighbor_distances_distribution(:, iter) = nearest_neighbor_distances;
        
        % Calculate the NNRI for this iteration
        NNRI_distribution(iter) = mean(nearest_neighbor_distances) / std(nearest_neighbor_distances);
    end
    
    % Optional: Visualize the points and their nearest neighbors in the final iteration
    if visualize
        figure;
        scatter(x_subset, y_subset, 'filled');
        hold on;

        % Plot lines to nearest neighbors
        for i = 1:num_random_points
            % Find the index of the nearest neighbor
            [~, nearest_idx] = min(sqrt((x_subset - x_subset(i)).^2 + (y_subset - y_subset(i)).^2));

            % Draw a line between the current point and its nearest neighbor
            plot([x_subset(i), x_subset(nearest_idx)], [y_subset(i), y_subset(nearest_idx)], 'r-');
        end

        xlabel('X');
        ylabel('Y');
        title('Nearest Neighbor Visualization (Final Iteration)');
        axis equal;
        hold off;
    end
end
