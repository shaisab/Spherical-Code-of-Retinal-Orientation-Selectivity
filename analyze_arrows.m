function analyze_arrows(center_polygon, arrow_centers, arrow_orientations)
    % Inputs:
    % center_polygon: [x_center, y_center] of the polygon (marked by a point)
    % arrow_centers: N-by-2 matrix where each row is [x_center, y_center] of an arrow
    % arrow_orientations: N-by-1 vector containing the orientation angle (in degrees or radians) of each arrow

    % Extract the center coordinates of the polygon
    x_polygon = center_polygon(1);
    y_polygon = center_polygon(2);

    % Number of arrows
    n_arrows = size(arrow_centers, 1);

    % Preallocate arrays for results
    distances = zeros(n_arrows, 1); % To store distances from polygon center
    theta_angles = zeros(n_arrows, 1); % To store theta angles
    arrow_orientation_angles = zeros(n_arrows, 1); % To store arrow orientation angles

    % Loop over each arrow
    for i = 1:n_arrows
        % Extract arrow center coordinates
        x_arrow = arrow_centers(i, 1);
        y_arrow = arrow_centers(i, 2);

        % a) Calculate the distance between polygon center and arrow center
        distances(i) = sqrt((x_arrow - x_polygon)^2 + (y_arrow - y_polygon)^2);

        % b) Calculate the angle theta between the horizontal line and the line connecting the center of the polygon and the center of the arrow
        delta_x = x_arrow - x_polygon;
        delta_y = y_arrow - y_polygon;
        theta_angles(i) = atan2(delta_y, delta_x); % atan2 returns the angle in radians
        
        % Convert to degrees if needed
        % theta_angles(i) = rad2deg(theta_angles(i));

        % c) Arrow orientation angle (this is provided in arrow_orientations)
        arrow_orientation_angles(i) = arrow_orientations(i); % Assuming input is in degrees or radians
    end

    % Display results
    fprintf('Arrow Analysis Results:\n');
    for i = 1:n_arrows
        fprintf('Arrow %d:\n', i);
        fprintf('  Distance to Polygon Center: %.2f\n', distances(i));
        fprintf('  Theta Angle (relative to horizontal): %.2f radians (%.2f degrees)\n', theta_angles(i), rad2deg(theta_angles(i)));
        fprintf('  Arrow Orientation Angle: %.2f\n', arrow_orientation_angles(i));
        fprintf('\n');
    end

    % Plot for visualization
    figure;
    hold on;

    % Plot the center of the polygon
    plot(x_polygon, y_polygon, 'ro', 'MarkerSize', 10, 'DisplayName', 'Polygon Center');

    % Plot arrows' centers and orientations
    for i = 1:n_arrows
        % Plot arrow center
        plot(arrow_centers(i, 1), arrow_centers(i, 2), 'bx', 'MarkerSize', 8, 'DisplayName', sprintf('Arrow %d Center', i));
        
        % Plot line from polygon center to arrow center
        plot([x_polygon, arrow_centers(i, 1)], [y_polygon, arrow_centers(i, 2)], 'k--', 'DisplayName', sprintf('Arrow %d Line', i));

        % Draw the arrow (as a line indicating its orientation)
        quiver(arrow_centers(i, 1), arrow_centers(i, 2), cos(arrow_orientations(i)), sin(arrow_orientations(i)), 0.1, 'r', 'LineWidth', 2, 'MaxHeadSize', 2);
    end

    % Formatting the plot
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Arrows and Polygon Center Analysis');
    legend('show');
    axis equal;
    grid on;
    hold off;
end
