function duplicateAllVars()
    % Duplicates all variables in the base workspace by appending their contents
    % along the same row or column, depending on their original orientation.
    %
    % Returns:
    %     None. All variables in the base workspace are updated in place.

    % Get the list of all variables in the base workspace
    vars = evalin('base', 'whos');

    for i = 1:numel(vars)
        % Get the variable name and its value
        varName = vars(i).name;
        varValue = evalin('base', varName);

        % Duplicate the variable contents
        if isstruct(varValue)
            % Recursively duplicate structures
            duplicatedValue = duplicateStruct(varValue);
        else
            % Duplicate non-structure variables
            duplicatedValue = duplicateValue(varValue);
        end

        % Update the variable in the base workspace
        assignin('base', varName, duplicatedValue);
    end
end

function duplicatedValue = duplicateValue(value)
    % Duplicates the contents of a variable by appending it along the same row or column.
    %
    % Args:
    %     value: The variable to duplicate (numeric, cell, char, etc.).
    %
    % Returns:
    %     duplicatedValue: The variable with its contents duplicated in the same row/column.

    if isnumeric(value) || islogical(value)
        if isrow(value)
            % Append along the same row
            duplicatedValue = [value, value];
        elseif iscolumn(value)
            % Append along the same column
            duplicatedValue = [value; value];
        else
            % For matrices, duplicate along rows (vertical concatenation)
            duplicatedValue = [value; value];
        end
    elseif ischar(value)
        % For strings, treat them as rows
        duplicatedValue = [value, value];
    elseif iscell(value)
        if isrow(value)
            duplicatedValue = [value, value];
        elseif iscolumn(value)
            duplicatedValue = [value; value];
        else
            duplicatedValue = [value; value];
        end
    else
        error('Unsupported variable type: %s', class(value));
    end
end

function duplicatedStruct = duplicateStruct(inputStruct)
    % Recursively duplicates the contents of a structure by appending fields.
    %
    % Args:
    %     inputStruct: The structure to duplicate.
    %
    % Returns:
    %     duplicatedStruct: The structure with fields duplicated recursively.

    duplicatedStruct = struct();

    % Iterate over all fields in the structure
    fields = fieldnames(inputStruct);
    for i = 1:numel(fields)
        fieldName = fields{i};
        fieldValue = inputStruct.(fieldName);

        % Check if the field is itself a structure
        if isstruct(fieldValue)
            % Recursively duplicate the nested structure
            duplicatedStruct.(fieldName) = duplicateStruct(fieldValue);
        else
            % Duplicate non-structure fields
            duplicatedStruct.(fieldName) = duplicateValue(fieldValue);
        end
    end
end
