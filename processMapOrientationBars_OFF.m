function [map,poriSz_bars_OFF] = processMapOrientationBars_OFF(map)
    % Function to process the map structure and update poriSz and related values

    % Initialize poriSz
    poriSz_bars_OFF = zeros(size(map.poriSz_bars_OFF, 1), 1);

    % Process poriSz based on conditions
    for i = 1:size(map.poriSz_bars_OFF, 1)
        if map.poriSz_bars_OFF(i, 1) >= 90
            poriSz_bars_OFF(i, 1) = map.poriSz_bars_OFF(i, 1) - 90;
        elseif map.poriSz_bars_OFF(i, 1) < 90
            poriSz_bars_OFF(i, 1) = map.poriSz_bars_OFF(i, 1) + 90;
        elseif isnan(map.poriSz_bars_OFF(i, 1))
            poriSz_bars_OFF(i, 1) = map.poriSz_bars_OFF(i, 1);
        end
    end

    % Update map.poriSz
    map.poriSz_bars_OFF = poriSz_bars_OFF;

    % Correct pori for quiver direction
    imporiSz_bars_OFF = 360 - poriSz_bars_OFF;

    % Extend poriSz by adding 180 degrees
    poriSzt_bars_OFF = poriSz_bars_OFF + 180;
    poriSz_bars_OFF = [poriSz_bars_OFF; poriSzt_bars_OFF];

    % Compute aori and bori
    aori_bars_OFF = cosd(poriSz_bars_OFF(:, 1));
    bori_bars_OFF = sind(poriSz_bars_OFF(:, 1));

    % Update map with new values
    map.xporiSz_bars_OFF = aori_bars_OFF;
    map.yporiSz_bars_OFF = bori_bars_OFF;
    map.SzMapOutputOS_bars_OFF = [[map.umroiXSznor; map.umroiXSznor], [map.umroiYSznor; map.umroiYSznor], map.xporiSz_bars_OFF, map.yporiSz_bars_OFF];
end
