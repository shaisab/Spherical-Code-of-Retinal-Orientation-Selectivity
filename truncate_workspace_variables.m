function truncate_workspace_variables(n)
    % Get all variables in the workspace
    vars = evalin('base', 'whos');
    
    % Iterate through each variable
    for i = 1:length(vars)
        var_name = vars(i).name;  % Variable name
        var_size = vars(i).size; % Variable size
        var_class = vars(i).class; % Variable type
        
        % Get the variable from the base workspace
        var_data = evalin('base', var_name);
        
        % Process based on type
        if isnumeric(var_data) || islogical(var_data) || iscell(var_data)
            % For arrays or cell arrays: truncate to the first n values
            if numel(var_data) > n
                var_data = var_data(1:min(n, var_size(1)), :); % Take first n rows
            end
            
        elseif isstruct(var_data)
            % For structures: truncate each field
            field_names = fieldnames(var_data);
            for j = 1:length(field_names)
                field_data = var_data.(field_names{j});
                if isnumeric(field_data) || islogical(field_data) || iscell(field_data)
                    if numel(field_data) > n
                        var_data.(field_names{j}) = field_data(1:min(n, size(field_data, 1)), :);
                    end
                end
            end
        end
        
        % Save the updated variable back to the base workspace
        assignin('base', var_name, var_data);
    end
    
    disp('Workspace variables truncated to the first n values.');
end
