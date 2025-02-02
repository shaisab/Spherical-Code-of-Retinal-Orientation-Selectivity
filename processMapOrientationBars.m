function [map,poriSz_bars] = processMapOrientationBars(map)
    % Function to process the map structure and update poriSz and related values

    % Initialize poriSz
    poriSz_bars = zeros(size(map.poriSz_bars, 1), 1);

    % Process poriSz based on conditions
    for i = 1:size(map.poriSz_bars, 1)
        if map.poriSz_bars(i, 1) >= 90
            poriSz_bars(i, 1) = map.poriSz_bars(i, 1) - 90;
        elseif map.poriSz_bars(i, 1) < 90
            poriSz_bars(i, 1) = map.poriSz_bars(i, 1) + 90;
        elseif isnan(map.poriSz_bars(i, 1))
            poriSz_bars(i, 1) = map.poriSz_bars(i, 1);
        end
    end

    % Update map.poriSz
    map.poriSz_bars = poriSz_bars;

    % Correct pori for quiver direction
    imporiSz_bars = 360 - poriSz_bars;

    % Extend poriSz by adding 180 degrees
    poriSzt_bars = poriSz_bars + 180;
    poriSz_bars = [poriSz_bars; poriSzt_bars];

    % Compute aori and bori
    aori_bars = cosd(poriSz_bars(:, 1));
    bori_bars = sind(poriSz_bars(:, 1));

    % Update map with new values
    map.xporiSz_bars = aori_bars;
    map.yporiSz_bars = bori_bars;
    map.SzMapOutputOS_bars = [[map.umroiXSznor; map.umroiXSznor], [map.umroiYSznor; map.umroiYSznor], map.xporiSz_bars, map.yporiSz_bars];
end
