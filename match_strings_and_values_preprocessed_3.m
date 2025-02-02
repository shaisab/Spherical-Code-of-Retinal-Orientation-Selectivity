function new_value_list = match_strings_and_values_preprocessed_3(long_list, short_list, value_list, pattern)
    % Initialize the new value list with zeros
    new_value_list = zeros(size(long_list));
    
    % Preprocess the short list to remove or replace the specified pattern
    for i = 1:length(short_list)
        short_list{i} = regexprep(short_list{i}, pattern, '');  % Remove pattern
    end
    
    % Create a logical array to track which strings in short_list have been matched
    matched = false(length(short_list), 1);
    
    % Loop through the long list and match with preprocessed short list
    for i = 1:length(long_list)
        long_str = long_list{i};  % Get the current string from the long list
        
        % Loop through the short list to find a matching substring
        for j = 1:length(short_list)
            short_str = short_list{j};  % Get the preprocessed string from the short list
            
            % Check if the long string contains the short string as a substring
            % and if the short string has not been matched yet
            if contains(long_str, short_str) && ~matched(j)
                % If a match is found, update the value in the new value list
                new_value_list(i) = value_list(j);
                matched(j) = true;  % Mark this string as matched
                break;  % Exit the loop once a match is found
            end
        end
    end
end
