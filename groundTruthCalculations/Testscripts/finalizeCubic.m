function dd = finalizeCubic(coords, coordsTime)
    size = length(coords);
    tmp = zeros(size, 1);
    dd = zeros(size, 1);

    % Set the second derivative at the end points to zero
    dd(1) = 0;
    dd(end) = 0;
    
    % Forward pass of calculating the second derivatives at each point
    for i = 2:size-1
        yrd = coords(i+1) - coords(i);
        yld = coords(i) - coords(i-1);
        xrd = coordsTime(i+1) - coordsTime(i);
        xld = coordsTime(i) - coordsTime(i-1);
        
        dr = yrd / xrd;
        dl = yld / xld;
        iw = coordsTime(i+1) - coordsTime(i-1);
        si = xld / iw;
        p = dd(i-1) * si + 2.0;
        
        dd(i) = (si - 1.0) / p;
        tmp(i) = ((dr - dl) * 6.0 / iw - tmp(i-1) * si) / p;
    end
    
    % Second (backward) pass of calculation
    for i = size-1:-1:1
        dd(i) = dd(i) * dd(i+1) + tmp(i);
    end
end
