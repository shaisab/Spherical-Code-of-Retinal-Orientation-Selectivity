% Ellipse parameters
x_c = -0.1; y_c = -0.1; % Center of the ellipses
a_base = 0.1; % Semi-major axis of the smallest ellipse
b_base = 0.05; % Semi-minor axis of the smallest ellipse
spacing = 0.1; % Distance between successive ellipses
num_ellipses = 4; % Number of concentric ellipses

% Find tangents
[deviation_angle, distance] = find_tangents_on_concentric_ellipses_x(x_c, y_c, a_base, b_base, spacing, num_ellipses);
