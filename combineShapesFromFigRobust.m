function combinedShape = combineShapesFromFigRobust(figFile)
% COMBINESHAPESFROMFIGROBUST Extracts and combines shapes from a .fig file into a single boundary.
%
% Usage:
%   combinedShape = combineShapesFromFigRobust('your_figure_file.fig');
%
% Input:
%   figFile - Name of the .fig file (including extension).
%
% Output:
%   combinedShape - A polyshape representing the outer perimeter of the combined shapes.
%
% This function recursively searches the figure for graphical objects, extracts
% their coordinates, and computes the combined perimeter.

    % Load the .fig file (invisible mode to prevent displaying it)
    fig = openfig(figFile, 'invisible');
    
    % Initialize arrays to store all X and Y coordinates
    allX = [];
    allY = [];
    
    % Recursive function to extract data from all children
    function extractFromObjects(objects)
        for obj = objects'
            if isprop(obj, 'XData') && isprop(obj, 'YData')
                % Line or similar objects with X and Y data
                xData = get(obj, 'XData');
                yData = get(obj, 'YData');
                if ~isempty(xData) && ~isempty(yData)
                    allX = [allX, xData(:)']; %#ok<*AGROW>
                    allY = [allY, yData(:)'];
                end
            elseif isprop(obj, 'Vertices') && ~isempty(get(obj, 'Vertices'))
                % Patch objects with vertices
                vertices = get(obj, 'Vertices');
                allX = [allX, vertices(:, 1)'];
                allY = [allY, vertices(:, 2)'];
            elseif isprop(obj, 'Children') && ~isempty(get(obj, 'Children'))
                % Recursively process child objects
                extractFromObjects(get(obj, 'Children'));
            end
        end
    end

    % Start extracting from the root figure
    extractFromObjects(get(fig, 'Children'));

    % Close the figure file to free resources
    close(fig);

    % Compute the boundary of the combined shape
    if ~isempty(allX) && ~isempty(allY)
        boundaryIndices = boundary(allX(:), allY(:)); % Boundary indices
        boundaryX = allX(boundaryIndices);
        boundaryY = allY(boundaryIndices);
        combinedShape = polyshape(boundaryX, boundaryY); % Create polyshape
    else
        combinedShape = polyshape(); % Return empty polyshape if no data
    end

    % Visualize the combined shape (optional)
    figure;
    plot(combinedShape);
    title('Combined Shape Perimeter');

    % Optionally, save the combined shape to a .mat file
    save('combined_shape.mat', 'combinedShape');
end
