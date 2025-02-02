function [new_polygon_vertices, updated_indices] = redrawPolygonAndUpdateIndices(map)
    % Example data (you already have points and initial indices)
    points = horzcat(map.umroiX,map.umroiY);  % Random points in a 10x10 space
    polygon_vertices = map.cluster_polygons{1, 3}.Vertices ;  % Initial polygon vertices

    % Plot the points and initial polygon
    fig = figure;
    scatter(points(:,1), points(:,2), 'k', 'filled'); % Plot points
    hold on;
    plot([polygon_vertices(:,1); polygon_vertices(1,1)], ...
         [polygon_vertices(:,2); polygon_vertices(1,2)], 'r-', 'LineWidth', 2); % Initial polygon
    title('Redraw Polygon to Encompass More Points');
    xlabel('X');
    ylabel('Y');
    grid on;
    hold off;

    % Store points and initial polygon in the figure's UserData
    data = struct('points', points, 'polygon_vertices', polygon_vertices, 'new_polygon_vertices', [], 'updated_indices', []);
    guidata(fig, data);

    % Allow user to draw a new polygon
    disp('Redraw the polygon around the points you want to label.');
    h_polygon = drawpolygon('Position', polygon_vertices, 'Color', 'r', 'LineWidth', 2);

    % Create a custom UI button to finalize polygon drawing
    uicontrol('Style', 'pushbutton', 'String', 'Finish Polygon', 'Position', [20, 20, 100, 30], ...
              'Callback', @(src, event) onFinishButtonPressed(fig, h_polygon));

    % Block until the user finishes drawing the polygon and presses the button
    uiwait(fig);  % Wait for the button press to close the figure

    % After the user finishes, the callback will return the outputs
    data = guidata(fig);
    new_polygon_vertices = data.new_polygon_vertices;  % Updated polygon vertices
    updated_indices = data.updated_indices;  % Updated indices of points inside the polygon
end

function onFinishButtonPressed(fig, h_polygon)
    % Retrieve the points and initial polygon from the figure's UserData
    data = guidata(fig);
    points = data.points;  % Retrieve the points
    polygon_vertices = data.polygon_vertices;  % Retrieve the initial polygon

    % Get the new position of the polygon after redrawing
    new_polygon_vertices = h_polygon.Position;

    % Recompute which points are inside the new polygon
    updated_indices = inpolygon(points(:,1), points(:,2), new_polygon_vertices(:,1), new_polygon_vertices(:,2));

    % Store the updated results in the figure's UserData for output
    data.new_polygon_vertices = new_polygon_vertices;
    data.updated_indices = updated_indices;
    guidata(fig, data);  % Update the figure's UserData

    % Close the figure once finished
    uiresume(fig);  % Allow the main function to continue execution
end

