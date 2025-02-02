% Choose the folder where the CSV files are located
folder_path = uigetdir('Select Folder Containing CSV Files');

% Get a list of all CSV files in the folder
csv_files = dir(fullfile(folder_path, '*.csv'));

% Preallocate a cell array to store data from each CSV file
data_tables = cell(length(csv_files), 1);

% Loop through each CSV file and extract the data
for i = 1:length(csv_files)
    % Get the file name
    file_name = fullfile(folder_path, csv_files(i).name);
    
    % Read the CSV file and store it in a table
    data_tables{i} = readtable(file_name);
    
    % Display a message confirming the file has been read
    fprintf('Data from %s has been read and stored.\n', csv_files(i).name);
end

% Optional: Combine all tables into one big table if they have the same structure
% combined_table = vertcat(data_tables{:});

% Now you can access individual tables via data_tables{i}, where i is the index of the file
