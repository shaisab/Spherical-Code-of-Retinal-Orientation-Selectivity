function [map,poriSz_bars_ON] = processMapOrientationBars_ON(map)
    % Function to process the map structure and update poriSz and related values

    % Initialize poriSz
    poriSz_bars_ON = zeros(size(map.poriSz_bars_ON, 1), 1);

    % Process poriSz based on conditions
    for i = 1:size(map.poriSz_bars_ON, 1)
        if map.poriSz_bars_ON(i, 1) >= 90
            poriSz_bars_ON(i, 1) = map.poriSz_bars_ON(i, 1) - 90;
        elseif map.poriSz_bars_ON(i, 1) < 90
            poriSz_bars_ON(i, 1) = map.poriSz_bars_ON(i, 1) + 90;
        elseif isnan(map.poriSz_bars_ON(i, 1))
            poriSz_bars_ON(i, 1) = map.poriSz_bars_ON(i, 1);
        end
    end

    % Update map.poriSz
    map.poriSz_bars_ON = poriSz_bars_ON;

    % Correct pori for quiver direction
    imporiSz_bars_ON = 360 - poriSz_bars_ON;

    % Extend poriSz by adding 180 degrees
    poriSzt_bars_ON = poriSz_bars_ON + 180;
    poriSz_bars_ON = [poriSz_bars_ON; poriSzt_bars_ON];

    % Compute aori and bori
    aori_bars_ON = cosd(poriSz_bars_ON(:, 1));
    bori_bars_ON = sind(poriSz_bars_ON(:, 1));

    % Update map with new values
    map.xporiSz_bars_ON = aori_bars_ON;
    map.yporiSz_bars_ON = bori_bars_ON;
    map.SzMapOutputOS_bars_ON = [[map.umroiXSznor; map.umroiXSznor], [map.umroiYSznor; map.umroiYSznor], map.xporiSz_bars_ON, map.yporiSz_bars_ON];
end
