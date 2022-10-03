function [value, elem_c] = subarea(data, pos_x, pos_y, distance)
    % the function estimates and avg value for the specified submatrix
    % using one dimension
    
    % get min/max indices %commented to sum only in one dimension
    min_y = max([pos_y - distance, 1]); 
    max_y = min([size(data,1), pos_y + distance]);
   
    %min_x = max([pos_x - distance, 1]);
    %max_x = min([size(data,2), pos_x + distance]);

    %submatrix = data(min_y:max_y, min_x:max_x);
    submatrix = data(min_y:max_y,pos_x);
    [r , c] = size(submatrix);
    
    newsubmatrix = db2pow(submatrix); % convert to linear space
    %submatrix(pos_y, pos_x) = 0.0; % remove the CUT value
    
    value = sum(newsubmatrix, 'all'); % get total sum
    elem_c = r * c;
    
end