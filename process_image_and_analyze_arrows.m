function process_image_and_analyze_arrows(image_path)
    % Load and display the image
    img = imread(image_path);
    figure;
    imshow(img);
    title('Click on the center of the polygon');

    % Manually select the center of the polygon
    [x_center, y_center] = ginput(1); % Click on the center of the polygon
    center_polygon = [x_center, y_center];
    hold on;
    plot(center_polygon(1), center_polygon(2), 'bo', 'MarkerSize', 10, 'LineWidth', 2);

    % Initialize empty arrays for manually marked arrow centers
    arrow_centers = [];
    orientation_angles = [];
    
    % Define a threshold distance to avoid duplicate points
    threshold_distance = 15; % Adjust this value based on image scale

    % Start a loop for manually marking arrow centers with feedback
    num_arrows = input('How many arrows do you want to mark? ');
    for i = 1:num_arrows
        title(sprintf('Click on the center of arrow %d', i));
        [x_arrow, y_arrow] = ginput(1); % Click on the center of an arrow

        % Check if the marked point is close to any previously marked point
        if ~isempty(arrow_centers)
            distances = sqrt((arrow_centers(:,1) - x_arrow).^2 + (arrow_centers(:,2) - y_arrow).^2);
            if any(distances < threshold_distance)
                fprintf('Point too close to a previously marked arrow. Please select a different point.\n');
                i = i - 1; % Decrease counter to allow another attempt
                continue;
            end
        end
        
        % If not too close, store the center point
        arrow_centers = [arrow_centers; x_arrow, y_arrow];
        plot(x_arrow, y_arrow, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
        
        % Get arrow orientation manually from user input (optional)
        orientation_angle = input(sprintf('Enter orientation angle for arrow %d (in degrees): ', i));
        orientation_angles = [orientation_angles; orientation_angle];

        % Provide feedback on the marked point
        fprintf('Arrow %d marked at (%.2f, %.2f) with orientation %.2f degrees.\n', ...
                i, x_arrow, y_arrow, orientation_angle);
    end

    % After marking, perform analysis
    analyze_arrows(center_polygon, arrow_centers, orientation_angles);
end

