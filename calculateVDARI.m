% Function to calculate the Voronoi Domain Area Regularity Index (VDARI)
function vdari = calculateVDARI(points, shape)
    % Ensure there are at least three points to compute the Voronoi diagram
    if size(points, 1) < 3
        vdari = NaN; % Not enough points to compute Voronoi domains
        return;
    end

    % Compute the Voronoi diagram for the points
    [~, v] = voronoin(points);

    % Initialize an array to store valid areas
    validAreas = [];

    % Iterate over each Voronoi cell
    for i = 1:length(v)
        % Skip if Voronoi cell is invalid (contains infinite points or indices out of bounds)
        if any(v{i} <= 0) || any(v{i} > size(points, 1))
            continue;
        end

        % Extract domain vertices
        domainVertices = points(v{i}, :);

        % Check if all vertices are within the shape
        if all(inpolygon(domainVertices(:, 1), domainVertices(:, 2), shape.Vertices(:, 1), shape.Vertices(:, 2)))
            % Calculate and store the area of the Voronoi domain
            validAreas = [validAreas; polyarea(domainVertices(:, 1), domainVertices(:, 2))];
        end
    end

    % Compute the VDARI if valid areas are found
    if isempty(validAreas) || std(validAreas) == 0
        vdari = NaN;  % Return NaN if VDARI cannot be computed
    else
        meanArea = mean(validAreas);
        stdArea = std(validAreas);
        vdari = meanArea / stdArea;
    end
end
