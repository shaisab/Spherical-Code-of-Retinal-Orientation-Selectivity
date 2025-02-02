function map = processMapOrientationForFit(map, j, ROIn, meanr_bars, A)
    % Check if the OSI and DSI value meets the threshold
        if map.OSI(j, 1) >= 0.2 && map.DSI(j, 1) < 0.17
            % Define orientation angles
            theta_8 = [0, 45, 90, 135, 180, 225, 270, 315];
            
            % Perform double Gaussian fit
            [oriFWHM, u, gof2] = doubleGaussianFit(ROIn, theta_8, meanr_bars);
            
            % Calculate preferred orientation size
            pori_bars = u;
            poriSz_bars = pori_bars + A;
            
            % Adjust for circular range constraints
            if poriSz_bars < 0
                poriSz_bars = 360 + poriSz_bars;
            end
            if poriSz_bars > 180
                poriSz_bars = poriSz_bars - 180;
            end
            
            % Update map with calculated value
            map.poriSz_bars(j, 1) = poriSz_bars;
        else
            % Assign NaN if conditions are not met
            map.poriSz_bars(j, 1) = NaN;
        end
                
        % Calculate degree difference if poriSz_bars is valid
        if abs(map.poriSz_bars(j, 1) - map.poriSz(j, 1)) > 90
            map.degdiff(j, 1) = 180 - abs(map.poriSz_bars(j, 1) - map.poriSz(j, 1));
        else
            map.degdiff(j, 1) = abs(map.poriSz_bars(j, 1) - map.poriSz(j, 1));
        end
end

