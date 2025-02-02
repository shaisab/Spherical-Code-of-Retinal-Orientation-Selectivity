% Number of points
n_points = 20;

% Generate random x coordinates between -128 and 128
x_coords = -128 + (128 - (-128)) * rand(1, n_points);

% Generate random y coordinates between -64 and 64
y_coords = -64 + (64 - (-64)) * rand(1, n_points);

% Combine x and y coordinates into an array of points
points = [x_coords; y_coords];

% Display the points
disp('Generated points (x, y):');
disp(points');