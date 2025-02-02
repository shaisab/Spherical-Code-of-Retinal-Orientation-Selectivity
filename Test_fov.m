clear all 
close all

% Load files - map overlap file
[fileName, PathName] = uigetfile('*.mat', 'Select map overlap file');
load(fileName);

map.updated_indices = [map.updated_indices;map.updated_indices];
type_OS_cluster_5 = (map.OSI>=0.2).* (map.DSI<=0.17).* (map.updated_indices==1).*(map.typeclust==8);
%type_OS_cluster_5 = (map.OSI>=0.2).* (map.DSI<=0.17);
% type_OS_cluster_5 = ones (6500,1);
% type_OS_cluster_5=type_OS_cluster_5;
list_length = length(type_OS_cluster_5);
type_non_OS_cluster_5 = ones (list_length,1);
type_non_OS_cluster_5 = (type_non_OS_cluster_5 - type_OS_cluster_5);

type_OS_cluster_5 = type_OS_cluster_5(1:(list_length)/2);
type_non_OS_cluster_5 = type_non_OS_cluster_5(1:(list_length)/2);

x_coords_OS_cluster_5 = map.umroiX .* type_OS_cluster_5;
x_coords_OS_cluster_5 = x_coords_OS_cluster_5(x_coords_OS_cluster_5 ~= 0);
y_coords_OS_cluster_5 = map.umroiY .* type_OS_cluster_5;
y_coords_OS_cluster_5 = y_coords_OS_cluster_5(y_coords_OS_cluster_5 ~= 0);

n_random_points = length(x_coords_OS_cluster_5);
x_coords = x_coords_OS_cluster_5;
y_coords = y_coords_OS_cluster_5;
[nearest_neighbor_distances, NNRI] = calculate_nearest_neighbors(x_coords, y_coords, n_random_points, true);

bounding_polygon = map.cluster_polygons (1,5);
bounding_polygon = bounding_polygon{1, 1};  

points = [x_coords, y_coords];
[vdari] = calculateVDARI(points, bounding_polygon);
[voronoi_areas, voronoi_vertices, vdari] = calculate_voronoi_vdari(x_coords, y_coords, n_random_points, bounding_polygon,true);

% x_coords_non_OS_cluster_5 = map.umroiX .* type_non_OS_cluster_5;
% x_coords_non_OS_cluster_5 = x_coords_non_OS_cluster_5(x_coords_non_OS_cluster_5 ~= 0);
% y_coords_non_OS_cluster_5 = map.umroiY .* type_non_OS_cluster_5;
% y_coords_non_OS_cluster_5 = y_coords_non_OS_cluster_5(y_coords_non_OS_cluster_5 ~= 0);
% 
% num_random_points = n_random_points;
% num_iterations = 1000;
% x_coords = x_coords_non_OS_cluster_5;
% y_coords = y_coords_non_OS_cluster_5;
% [NNRI_distribution, nearest_neighbor_distances_distribution] = calculate_nnri_distribution(x_coords, y_coords, num_random_points, num_iterations);
% [all_voronoi_areas, all_voronoi_vertices, vdari_values] = calculate_voronoi_distribution(x_coords, y_coords, n_random_points, bounding_polygon, num_iterations);
