clear all
close all

% Load files - all cells file
[fileName, PathName] = uigetfile('*.mat', 'Select all cells file');
load(fileName);

% Load files 2 - map file
[fileName_2, PathName_2] = uigetfile('*.mat', 'Select map file');
load(fileName_2);

% Load files 3 - 
[fileName_3, PathName_3] = uigetfile('*.txt', 'Select xy position file');
position_table = readtable(fileName_3);

% Specify the substring you want to filter by
filterSubstring = 'Aug_23_2015';

% Specify the indices you want to filter by
filterIndices = [3, 4, 6, 8];

% Specify the new list (e.g., some descriptive names or labels)
newList = {'ON Sus', 'ON OFF Sus', 'ON OFF Trans', 'ON Trans'};

filtered_data = filterDataBySubstringAndIndex(fileName, filterSubstring, filterIndices, newList);

long_list_size = length(map.cellID)/2;
long_list = map.cellID (1:long_list_size);
short_list_size = length(filtered_data.name)/2;
short_list = filtered_data.name(1:short_list_size);
value_list = filtered_data.index(1:short_list_size);

% Define the pattern to add
pattern = 'Aug_23_2015_';  % This adds the date portion

[new_value_list, numerical_value_list] = match_strings_and_values_preprocessed_4(long_list, short_list, value_list, pattern);

map.typeclust = [new_value_list;new_value_list];
map.fovtype = [numerical_value_list;numerical_value_list];

[cluster_polygons,cluster_indices,combined_area,combined_perimeter] = process_and_analyze_squares(position_table, map);
map.square_cluster_indices = [cluster_indices;cluster_indices];
map.combined_fov_area = combined_area;
map.combined_perimeter = combined_perimeter;
map.cluster_polygons = cluster_polygons;
