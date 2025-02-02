function [deviation_angle, distance] = find_tangents_on_concentric_ellipses_x(x_c, y_c, a_base, b_base, spacing, num_ellipses)
    % Input:
    % x_c, y_c: Center of ellipses (Cartesian coordinates)
    % a_base: Semi-major axis of the smallest ellipse
    % b_base: Semi-minor axis of the smallest ellipse
    % spacing: Distance between successive ellipses
    % num_ellipses: Number of concentric ellipses
    %
    % Output:
    % deviation_angle: List of deviation angles from the y-axis for each ellipse
    % distance: List of distances from the center to the intersection points

    % Generate semi-major and semi-minor axes for concentric ellipses
    a_values = a_base + (0:(num_ellipses-1)) * spacing; % Semi-major axes
    b_values = b_base + (0:(num_ellipses-1)) * spacing; % Semi-minor axes
    
    % Define the slopes of the two orthogonal lines
    slopes = [1, -1]; % Lines y = x and y = -x
    
    % Initialize arrays to store deviation angles and distances
    deviation_angle = NaN(1, num_ellipses); 
    distance = NaN(1, num_ellipses);
    
    % Iterate through each ellipse
    for k = 1:num_ellipses
        a = a_values(k); % Semi-major axis
        b = b_values(k); % Semi-minor axis
        
        % Iterate through each line
        for m_line = slopes
            % Define the function for the ellipse equation minus the line equation
            % (x - x_c)^2 / a^2 + (m_line * x - y_c)^2 / b^2 - 1 = 0
            ellipse_eq = @(x) ((x - x_c)^2 / a^2) + ((m_line * x - y_c)^2 / b^2) - 1;
            
            % Use fsolve to find the root numerically (real intersections)
            options = optimset('Display', 'off'); % Suppress output
            x_intersections = fsolve(ellipse_eq, 0, options);  % Initial guess is 0 (can be adjusted)

            % Loop through the intersections (although there should typically be 1 or 2)
            for i = 1:length(x_intersections)
                x_int = x_intersections(i); % Intersection x-coordinate
                y_int = m_line * x_int;     % Intersection y-coordinate
                
                % Tangent slope calculation
                if abs(y_int - y_c) > 1e-6 % Avoid division by zero
                    tangent_slope = -((b^2 / a^2) * (x_int - x_c)) / (y_int - y_c);
                else
                    tangent_slope = Inf; % Vertical tangent
                end
                
                % Deviation from y-axis (angle in degrees)
                if isfinite(tangent_slope) && isreal(tangent_slope)
                    deviation_angle(k) = atand(1 / tangent_slope);
                else
                    deviation_angle(k) = NaN; % Invalid angle
                end
                
                % Distance from the center of the ellipses
                if isreal(x_int) && isreal(y_int)
                    distance(k) = sqrt((x_int - x_c)^2 + (y_int - y_c)^2);
                else
                    distance(k) = NaN; % Invalid distance if complex
                end
            end
        end
    end
end
