% Visualize the VDARI distribution
figure;
edges=[0:5000:40000];
histogram(voronoi_areas, edges);
title('Voronoi Domain Area size Distribution');
xlabel('μm^2');
ylabel('Frequency');

valuesToMark = [vdari];  % Values to mark on the histogram
labels = {'vdari'};  % Labels for the values
lineColors = {'black'};  % Colors for the lines

% Add vertical lines, labels, and customize colors
ylimits = ylim;  % Get the y-axis limits
middleY = (ylimits(2) - ylimits(1)) / 2;  % Middle of the y-axis

for i = 1:length(valuesToMark)
    % Add vertical line
    xline(valuesToMark(i), 'Color', lineColors{i}, 'LineWidth', 2);  
    
     % Create label text with the value included
    labelText = sprintf('%s = %.2f', labels{i}, valuesToMark(i));  % Format the label with the value
    
    % Add label near the middle of the line
    text(valuesToMark(i), middleY, labelText, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'Color', lineColors{i}, 'FontSize', 12, 'FontWeight', 'bold');
end

hold off;