
% Define the cubic interpolation function
function coord = cubic_interpolation(t, waypoints, coords, dd)
    % Find the appropriate segment
    xrp = find(waypoints > t, 1); % Upper bound
    if isempty(xrp)
        xrp = length(waypoints); % Right endpoint
    end
    
    if xrp == 1
        coord = coords(1); % Left endpoint
    elseif xrp == length(waypoints)
        coord = coords(end); % Right endpoint
    else
        xli = xrp - 1;
        xrd = (waypoints(xrp) - t); 
        xld = (t - waypoints(xli));
        iw = (waypoints(xrp) - waypoints(xli));
        iws = iw * iw / 6;
        
        A = xrd / iw;
        B = xld / iw;
        C = (A * A * A - A) * iws;
        D = (B * B * B - B) * iws;
        
        coord = coords(xli) * A + coords(xrp) * B + dd(xli) * C + dd(xrp) * D;
    end
end