% Visualize the NNRI distribution
figure;
edges=[0:0.25:15];
histogram(nnriDistribution, edges);
title('NNRI Distribution');
xlabel('NNRI');
ylabel('Frequency');

valuesToMark = [median(nnriDistribution),random_nnri,nnri_square];  % Values to mark on the histogram
labels = {'Median','Random','Simulated Field'};  % Labels for the values
lineColors = {'red','green','black'};  % Colors for the lines

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

% Value to find the percentile of
random_vdari = 1.91;

% Sort the list
sorted_vdariDistribution = sort(vdariDistribution);

% Calculate the percentile rank
percentile_random_vdari = 100 * sum(sorted_vdariDistribution < random_vdari) / numel(vdariDistribution);

% Visualize the VDARI distribution
figure;
edges=[0:0.1:5];
histogram(vdariDistribution, edges);
title('VDARI Distribution');
xlabel('VDARI');
ylabel('Frequency');

valuesToMark = [median(vdariDistribution),random_vdari,vdari_square];  % Values to mark on the histogram
labels = {'Median','Random','Simulated Field'};  % Labels for the values
lineColors = {'red','green','black'};  % Colors for the lines

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