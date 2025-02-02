function [cluster_polygons,cluster_indices,combined_area,combined_perimeter] = process_and_analyze_squares(position_table, map)
    % Rename the variables in the position table
    position_table = renamevars(position_table, ["Var1","Var2","Var3","Var4"], ...
                     ["fov","relXPosition","relYPosition","relZPosition"]);

    % Number of squares
    n_squares = length(position_table.relXPosition);

    % Extract center coordinates from the data
    x_centers = position_table.relXPosition;
    y_centers = position_table.relYPosition;
    
    % Define the widths and heights of the squares
    widths = repmat(256, n_squares, 1);    % Example width for all squares
    heights = repmat(128, n_squares, 1);   % Example height for all squares

    % Combine into a single matrix
    squares = [x_centers, y_centers, widths, heights];

    % Initialize overlap matrix
    overlap_matrix = zeros(n_squares, n_squares);

    % Calculate overlap matrix
    for i = 1:n_squares
        for j = i+1:n_squares
            % Get the details of each square
            x1 = squares(i, 1); y1 = squares(i, 2); w1 = squares(i, 3); h1 = squares(i, 4);
            x2 = squares(j, 1); y2 = squares(j, 2); w2 = squares(j, 3); h2 = squares(j, 4);

            % Calculate overlap area
            overlap_area_ij = overlap_area(x1, y1, w1, h1, x2, y2, w2, h2);

            % Store in matrix
            overlap_matrix(i, j) = overlap_area_ij;
            overlap_matrix(j, i) = overlap_area_ij; % Symmetric matrix
        end
    end

    % Calculate clusters of overlapping squares
    adj_matrix = overlap_matrix > 0;
    G = graph(adj_matrix);
    bins = conncomp(G);

    % Calculate combined area and perimeter for each cluster
    n_clusters = max(bins);
    combined_area = zeros(1, n_clusters);
    combined_perimeter = zeros(1, n_clusters);
    cluster_polygons = cell(1, n_clusters);

    for cluster_idx = 1:n_clusters
        cluster_squares_idx = find(bins == cluster_idx);
        combined_polygon = polyshape();

        for k = 1:length(cluster_squares_idx)
            idx = cluster_squares_idx(k);
            x_center = squares(idx, 1);
            y_center = squares(idx, 2);
            width = squares(idx, 3);
            height = squares(idx, 4);

            half_width = width / 2;
            half_height = height / 2;

            x_coords = [x_center - half_width, x_center + half_width, ...
                        x_center + half_width, x_center - half_width];
            y_coords = [y_center - half_height, y_center - half_height, ...
                        y_center + half_height, y_center + half_height];

            poly_k = polyshape(x_coords, y_coords);
            combined_polygon = union(combined_polygon, poly_k);
        end

        combined_area(cluster_idx) = area(combined_polygon);
        combined_perimeter(cluster_idx) = perimeter(combined_polygon);
        cluster_polygons{cluster_idx} = combined_polygon;
    end

    % Visualize the square clusters with area and perimeter
    figure;
    hold on;

    for cluster_idx = 1:n_clusters
        plot(cluster_polygons{cluster_idx}, 'FaceColor', 'none', 'EdgeColor', 'k');

        [centroid_x, centroid_y] = centroid(cluster_polygons{cluster_idx});
        label_text = sprintf('Cluster %d\nArea = %.2f\nPerimeter = %.2f', ...
             cluster_idx, combined_area(cluster_idx), combined_perimeter(cluster_idx));
        text(centroid_x, centroid_y, label_text, ...
             'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'BackgroundColor', 'white');
    end

    axis equal;
    xlabel('X');
    ylabel('Y');
    title('Square Clusters with Area and Perimeter');
    hold off;

    % Get the coordinates to check from the map data
    x_coords_to_check = map.umroiX;
    y_coords_to_check = map.umroiY;

    % Check which cluster each coordinate belongs to
    cluster_indices = zeros(length(x_coords_to_check), 1); % Initialize output list

    for point_idx = 1:length(x_coords_to_check)
        test_x = x_coords_to_check(point_idx);
        test_y = y_coords_to_check(point_idx);

        found_cluster = false;
        for cluster_idx = 1:n_clusters
            if isinterior(cluster_polygons{cluster_idx}, test_x, test_y)
                cluster_indices(point_idx) = cluster_idx;
                found_cluster = true;
                break; % Stop checking once we find the cluster
            end
        end

        if ~found_cluster
            cluster_indices(point_idx) = 0; % If not in any cluster
        end
    end
end