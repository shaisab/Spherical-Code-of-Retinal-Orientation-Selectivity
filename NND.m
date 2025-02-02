% Define the number of points
n_points = length(map.umroiX);  % Example number of points

% Generate random x and y coordinates
x_coords = map.umroiX;  % Example x coordinates
y_coords = map.umroiY;    % Example y coordinates

% Initialize an array to store the nearest neighbor distance for each point
nearest_neighbor_distances = zeros(n_points, 1);

% Calculate the nearest neighbor distance for each point
for i = 1:n_points
    % Get the current point
    x_i = x_coords(i);
    y_i = y_coords(i);
    
    % Calculate distances from the current point to all other points
    distances = sqrt((x_coords - x_i).^2 + (y_coords - y_i).^2);
    
    % Ignore the distance from the point to itself (set it to Inf)
    distances(i) = Inf;
    
    % Find the minimum distance, which is the nearest neighbor distance
    nearest_neighbor_distances(i) = min(distances);
end

% Display the nearest neighbor distances
disp('Nearest Neighbor Distances:');
disp(nearest_neighbor_distances);

% Optional: Visualize the points and their nearest neighbors
figure;
scatter(x_coords, y_coords, 'filled');
hold on;

% Plot lines to nearest neighbors
for i = 1:n_points
    % Find the index of the nearest neighbor
    [~, nearest_idx] = min(sqrt((x_coords - x_coords(i)).^2 + (y_coords - y_coords(i)).^2));
    
    % Draw a line between the current point and its nearest neighbor
    plot([x_coords(i), x_coords(nearest_idx)], [y_coords(i), y_coords(nearest_idx)], 'r-');
end

xlabel('X');
ylabel('Y');
title('Nearest Neighbor Visualization');
axis equal;

plot(map.cluster_polygons{1, 7}  )
hold off;