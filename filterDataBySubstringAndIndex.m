function filtered_data = filterDataBySubstringAndIndex(fileName, filterSubstring, filterIndices, newList)
    % Function to filter data based on a substring and indices
    % Inputs:
    %   - fileName: Name of the .mat file to load
    %   - filterSubstring: Substring to filter the 'name' field
    %   - filterIndices: Indices to filter the 'index' field
    %   - newList: List of labels corresponding to filterIndices
    % Output:
    %   - filtered_data: Structure with filtered 'name', 'index', and 'label'
    
    % Load the .mat file
    if nargin < 1 || isempty(fileName)
        [FileName, PathName] = uigetfile('*.mat', 'Select map position file');
        if isequal(FileName,0)
            error('No file selected.');
        end
        fileName = fullfile(PathName, FileName);
    end
    load(fileName);
    
    % Define data structure
    data.name = allCells.cellID;
    data.index = allCells.typeClust;
    
    % Initialize new structure to hold the filtered data
    filtered_data.name = {};
    filtered_data.index = [];
    filtered_data.label = {};  % New field for the additional list
    
    % Debugging: Display initial data
    disp('Initial data:');
    disp(table(data.name', data.index', 'VariableNames', {'Name', 'Index'}));
    
    % Filter based on substring and index
    for i = 1:length(data.name)
        % Debugging: Check current entry
        fprintf('Checking entry %d: %s\n', i, data.name{i});
        
        % Check if the name contains the filter substring and index is within the filter
        if contains(data.name{i}, filterSubstring, 'IgnoreCase', true) && ismember(data.index(i), filterIndices)
            % Append to filtered_data if conditions are met
            filtered_data.name{end+1} = data.name{i};
            filtered_data.index(end+1) = data.index(i);
            
            % Find corresponding label from newList
            labelIdx = find(filterIndices == data.index(i), 1);
            if ~isempty(labelIdx) && labelIdx <= length(newList)
                filtered_data.label{end+1} = newList{labelIdx};
            else
                filtered_data.label{end+1} = 'Unknown';  % Default label if index is not found
            end
        end
    end
    
    filtered_data.name = filtered_data.name'; 
    filtered_data.index = filtered_data.index';  
    filtered_data.label = filtered_data.label'; 
    
    % Display the filtered structure
    disp('Filtered data:');
    disp(table(filtered_data.name, filtered_data.index, filtered_data.label, 'VariableNames', {'Name', 'Index', 'Label'}));
end
