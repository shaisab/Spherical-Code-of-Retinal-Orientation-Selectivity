% Function to calculate overlap area between two rectangles
function area = overlap_area(x1, y1, w1, h1, x2, y2, w2, h2)
    % Calculate edges of the rectangles
    x1_min = x1 - w1 / 2;
    x1_max = x1 + w1 / 2;
    y1_min = y1 - h1 / 2;
    y1_max = y1 + h1 / 2;
    
    x2_min = x2 - w2 / 2;
    x2_max = x2 + w2 / 2;
    y2_min = y2 - h2 / 2;
    y2_max = y2 + h2 / 2;
    
    % Calculate intersection boundaries
    x_int_min = max(x1_min, x2_min);
    x_int_max = min(x1_max, x2_max);
    y_int_min = max(y1_min, y2_min);
    y_int_max = min(y1_max, y2_max);
    
    % Calculate intersection area
    if x_int_min < x_int_max && y_int_min < y_int_max
        int_width = x_int_max - x_int_min;
        int_height = y_int_max - y_int_min;
        area = int_width * int_height;
    else
        area = 0;
    end
end