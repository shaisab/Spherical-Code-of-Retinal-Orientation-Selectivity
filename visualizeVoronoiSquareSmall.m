function visualizeVoronoiSquareSmall(points, squareSideLength, smallShape)
    % Ensure there are at least three points to compute the Voronoi diagram
    if size(points, 1) < 3
        error('At least three points are required to compute Voronoi domains.');
    end

    % Compute the Voronoi diagram for the points
    [v, voronoiCells] = voronoin(points);

    % Create a figure to visualize the Voronoi diagram
    figure;
    hold on;
    
    % Plot the square boundary (large square)
    rectangle('Position', [0, 0, squareSideLength, squareSideLength], ...
              'EdgeColor', 'k', 'LineWidth', 2);

    % Plot the Voronoi cells
    for i = 1:length(voronoiCells)
        cellVertices = voronoiCells{i};
        
        % Skip invalid or infinite Voronoi cells
        if any(cellVertices == 1) || isempty(cellVertices)
            continue;
        end
        
        % Get the coordinates of the cell vertices
        vertexCoords = v(cellVertices, :);
        
        % Clip the vertices to the boundaries of the square
        vertexCoords(:, 1) = max(0, min(squareSideLength, vertexCoords(:, 1)));
        vertexCoords(:, 2) = max(0, min(squareSideLength, vertexCoords(:, 2)));
        
        % Check if all vertices are within bounds after clipping
        if all(vertexCoords(:, 1) >= 0 & vertexCoords(:, 1) <= squareSideLength & ...
               vertexCoords(:, 2) >= 0 & vertexCoords(:, 2) <= squareSideLength)
            % Plot the Voronoi cell boundary
            plot([vertexCoords(:, 1); vertexCoords(1, 1)], ...
                 [vertexCoords(:, 2); vertexCoords(1, 2)], 'b-', 'LineWidth', 1);
        end
    end
    
    % Plot the points
    scatter(points(:, 1), points(:, 2), 50, 'r', 'filled');
    
    smallShape=smallShape.Vertices;

    % Overlay the shape boundary
    plot([smallShape(:, 1); smallShape(1, 1)], [smallShape(:, 2); smallShape(1, 2)], 'k-', 'LineWidth', 2);
    
    % Label the plot
    title('Voronoi Diagram for Points within a Square');
    xlabel('X');
    ylabel('Y');
    axis equal;
    xlim([0, squareSideLength]);
    ylim([0, squareSideLength]);
    hold off;
end
