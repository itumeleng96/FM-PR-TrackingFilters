% This script calculates the Bistatic range(m)  and Doppler(Hz) From FERS
% Simulation
% Author: Itumeleng Malemela

%The Positions are all in the format Below
%[x,y,z] where x,y,z are the coordinates in meters in the Cartesian Coordinate System

wavelength = 299792458/94e6;

%Transmitter Position
Tx_Pos=[6440;10760;1056];
%Reference Reciever Position 
RefRx_Pos=[0;0;1000];

%Surveillance Reciever Position 
SurvRx_Pos=[0;1;1000];

%Target Positions
Target_Pos_1 = [4000;18000;3600];
Target_Pos_2 = [-2000;3000;1600];

Baseline = norm(RefRx_Pos - Tx_Pos);

%Duration of The Simulation and the update intterval in seconds
Simulation_time = 60;
UpdateInterval = 1;

% Calculate the change in position
delta_position = Target_Pos_2 - Target_Pos_1;

% Calculate the change in time
delta_time = 60;

% Calculate the speed in each dimension
speed_x = delta_position(1) / delta_time;
speed_y = delta_position(2) / delta_time;
speed_z = delta_position(3) / delta_time;

speed_xyz = [speed_x;speed_y;speed_z];
% Calculate the overall speed (magnitude)
speed_magnitude = norm(delta_position) / delta_time;

% Display the calculated speed
fprintf('Overall speed: %.2f m/s\n', speed_magnitude);

% Calculate the number of time steps
num_steps = Simulation_time / UpdateInterval;

% Initialize a matrix to store the estimated positions
estimated_positions = zeros(3, num_steps);

% Estimate the positions at each time step
% Initialize a matrix to store the estimated bistatic ranges
bistatic_ranges = zeros(1, num_steps);
bistatic_doppler_shifts = zeros(1, num_steps);

prev_total_range =0;

% Estimate the bistatic ranges at each time step
for step = 1:num_steps+1
    estimated_position = Target_Pos_1 +   speed_xyz* UpdateInterval * (step-1);
    
    % Calculate the distance between target and transmitter
    range_tx_target = norm(estimated_position - Tx_Pos);
    
    % Calculate the distance between target and reference receiver
    range_ref_target = norm(estimated_position - RefRx_Pos);
    
    % Calculate the distance between target and surveillance receiver
    range_surv_target = norm(estimated_position - SurvRx_Pos);
    
    total_range = range_tx_target + range_ref_target;

    % Calculate the bistatic ranges
    bistatic_range_ref = range_tx_target + range_ref_target - Baseline;
    bistatic_range_surv = range_tx_target + range_surv_target - Baseline;
    
    % Store the bistatic ranges in the matrix
    bistatic_ranges(step) = bistatic_range_ref;

    %Calculate the rate of change of the total range with respect to time
    d_total_range_dt = (total_range - prev_total_range) / UpdateInterval;
    prev_total_range = total_range;
    
    % Calculate the Doppler shift
    doppler_shift = - (1 / wavelength) * (d_total_range_dt);
    
    bistatic_doppler_shifts(step)=doppler_shift;

end


% Save bistatic ranges and Doppler shifts to an HDF5 file
h5create('output_data.h5', '/bistatic_ranges', size(bistatic_ranges));
h5write('output_data.h5', '/bistatic_ranges', bistatic_ranges);

h5create('output_data.h5', '/doppler_shifts', size(bistatic_doppler_shifts));
h5write('output_data.h5', '/doppler_shifts', bistatic_doppler_shifts);
