function coord = getPositionCubic(t, coords, coordsTime, dd)
    % Find the index of the first time greater than t
    xrp = find(coordsTime >= t, 1);
    
    % Check if t is beyond the left endpoint
    if xrp == 1
        coord = coords(:, 1); % Set coord to the left endpoint
    % Check if t is beyond the right endpoint
    elseif xrp == numel(coordsTime) + 1
        coord = coords(:, end-1); % Set coord to the right endpoint
    else
        % Perform cubic spline interpolation
        xri = xrp; % Index of right coordinate
        xli = xri - 1; % Index of left coordinate
        
        % Compute distances from t to left and right coordinates
        xrd = (coordsTime(xri) - t);
        xld = (t - coordsTime(xli));
        
        % Compute weight between left and right coordinates
        iw = (coordsTime(xri) - coordsTime(xli));
        iws = iw ^ 2 / 6.0;
        
        % Compute coefficients for cubic spline interpolation
        A = xrd / iw;
        B = xld / iw;
        C = (A ^ 3 - A) * iws;
        D = (B ^ 3 - B) * iws;
        
        % Perform cubic spline interpolation to find coord
        coord = coords(:, xli) * A + coords(:, xri) * B + dd(xli) * C + dd(xri) * D;
    end
    % Set the time part of the coordinate
    % coord(3) = t;
end
