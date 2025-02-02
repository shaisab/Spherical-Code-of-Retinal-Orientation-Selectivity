function truncate_workspace_variables_recursive(n)
    % Get all variables in the workspace
    vars = evalin('base', 'whos');
    
    % Iterate through each variable
    for i = 1:length(vars)
        var_name = vars(i).name;  % Variable name
        var_data = evalin('base', var_name);  % Get the variable
        
        % Process the variable with truncation
        truncated_var = truncate_variable(var_data, n);
        
        % Assign the truncated variable back to the base workspace
        assignin('base', var_name, truncated_var);
    end
    
    disp('Workspace variables truncated to the first n values (including nested structures).');
end

function truncated_var = truncate_variable(var_data, n)
    % Base case: Handle numeric, logical, or cell arrays
    if isnumeric(var_data) || islogical(var_data) || iscell(var_data)
        % Truncate the array to the first n rows if it's 2D
        truncated_var = var_data(1:min(n, size(var_data, 1)), :);
    
    % Recursive case: Handle structures
    elseif isstruct(var_data)
        % Iterate over each field in the structure
        field_names = fieldnames(var_data);
        for j = 1:length(field_names)
            field_name = field_names{j};
            field_data = var_data.(field_name);
            
            % Recursively truncate the field
            var_data.(field_name) = truncate_variable(field_data, n);
        end
        truncated_var = var_data;
    
    % Unsupported data types are left unchanged
    else
        truncated_var = var_data;
    end
end
