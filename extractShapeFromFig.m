function shapeData = extractShapeFromFig(figFile)
% EXTRACTSHAPEFROMFIG Extracts shapes from a .fig file and returns their data.
%
% Usage:
%   shapeData = extractShapeFromFig('your_figure_file.fig');
%
% Input:
%   figFile - Name of the .fig file (including extension)
%
% Output:
%   shapeData - A structure array containing extracted shape data
%
% The function identifies line, surface, and patch objects within the figure
% and extracts their respective data.

    % Load the .fig file (invisible mode to prevent displaying it)
    fig = openfig(figFile, 'invisible');
    
    % Find all axes in the figure
    axesHandles = findall(fig, 'type', 'axes');

    % Initialize a structure array to store shape data
    shapeData = struct('Type', {}, 'XData', {}, 'YData', {}, 'ZData', {}, ...
                       'Vertices', {}, 'Faces', {});

    % Loop through all axes to extract graphical objects
    for i = 1:length(axesHandles)
        % Find line objects (2D shapes)
        lineObjects = findall(axesHandles(i), 'type', 'line');
        for j = 1:length(lineObjects)
            shapeData(end+1).Type = 'line'; %#ok<*AGROW>
            shapeData(end).XData = get(lineObjects(j), 'XData');
            shapeData(end).YData = get(lineObjects(j), 'YData');
            shapeData(end).ZData = []; % Lines are typically 2D
        end

        % Find surface objects (3D shapes)
        surfaceObjects = findall(axesHandles(i), 'type', 'surface');
        for j = 1:length(surfaceObjects)
            shapeData(end+1).Type = 'surface';
            shapeData(end).XData = get(surfaceObjects(j), 'XData');
            shapeData(end).YData = get(surfaceObjects(j), 'YData');
            shapeData(end).ZData = get(surfaceObjects(j), 'ZData');
        end

        % Find patch objects (polygons)
        patchObjects = findall(axesHandles(i), 'type', 'patch');
        for j = 1:length(patchObjects)
            shapeData(end+1).Type = 'patch';
            shapeData(end).Vertices = get(patchObjects(j), 'Vertices');
            shapeData(end).Faces = get(patchObjects(j), 'Faces');
            shapeData(end).XData = []; % Not applicable for patches
            shapeData(end).YData = [];
            shapeData(end).ZData = [];
        end
    end

    % Close the figure file to free resources
    close(fig);

    % Example visualization (optional, can be commented out if not needed)
    for i = 1:length(shapeData)
        if strcmp(shapeData(i).Type, 'line')
            % Plot 2D lines
            figure; 
            plot(shapeData(i).XData, shapeData(i).YData);
            title(['2D Line Object ', num2str(i)]);
        elseif strcmp(shapeData(i).Type, 'surface')
            % Plot 3D surfaces
            figure; 
            mesh(shapeData(i).XData, shapeData(i).YData, shapeData(i).ZData);
            title(['3D Surface Object ', num2str(i)]);
        elseif strcmp(shapeData(i).Type, 'patch')
            % Visualize patch objects
            figure; 
            patch('Vertices', shapeData(i).Vertices, 'Faces', shapeData(i).Faces, ...
                  'FaceColor', 'cyan', 'EdgeColor', 'black');
            title(['Patch Object ', num2str(i)]);
        end
    end

    % Optionally, save the shape data to a .mat file
    save('extracted_shape_data.mat', 'shapeData');
end
