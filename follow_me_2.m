

function follow_me_2
    h=gcf;              
    c = camlight('right');    % Create light
    set(c,'style','infinite');    % Set style
    h = rotate3d;                 % Create rotate3d-handle
    h.ActionPostCallback = @RotationCallback; % assign callback-function
    h.Enable = 'on';              % no need to click the UI-button

    % Sub function for callback
    function RotationCallback(~,~)
        c = camlight(c,'right');
    end
end