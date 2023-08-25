% Read bistatic ranges and Doppler shifts from the HDF5 file
bistatic_ranges = h5read('output_data.h5', '/bistatic_ranges');
doppler_shifts = h5read('output_data.h5', '/doppler_shifts');

% Print the bistatic ranges
fprintf('Bistatic Ranges:\n');
disp(bistatic_ranges);

% Print the Doppler shifts
fprintf('Doppler Shifts:\n');
disp(doppler_shifts);