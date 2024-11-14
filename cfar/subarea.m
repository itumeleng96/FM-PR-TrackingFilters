function [value, elem_c] = subarea(data, pos_x, pos_y, y_distance)    
    % Multi-Deimonsional CFAR 
    % function [value, elem_c] = subarea(data, pos_x, pos_y, y_distance,x_distance)
    %
    % The function estimates and avg value for the specified submatrix
    % in 1-D for Doppler and Range Map (Range-resolution low for FM)
    
    % get min/max indices %commented to sum only in one dimension
    min_y = max([pos_y - y_distance, 1]); 
    max_y = min([size(data,1), pos_y + y_distance]);
   
    %min_x = max([pos_x - x_distance, 1]);
    %max_x = min([size(data,2), pos_x + x_distance]);

    %submatrix = data(min_y:max_y, min_x:max_x);
    submatrix = data(min_y:max_y,pos_x);
    [r , c] = size(submatrix);
    
    newsubmatrix = submatrix; % convert to linear space
    %submatrix(pos_y, pos_x) = 0.0; % remove the CUT value
    
    value = sum(newsubmatrix, 'all'); % get total sum
    elem_c = r * c;
    
end