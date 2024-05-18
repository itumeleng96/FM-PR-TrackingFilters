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


Time_for_waypoints =[0,60];

Target_waypoints = [4000,18000,3600;
                    -2000,3000,1600];



Baseline = norm(RefRx_Pos - Tx_Pos);

%Duration of The Simulation and the update intterval in seconds
Simulation_time = Time_for_waypoints(end)-Time_for_waypoints(1);
UpdateInterval = 1;

bistatic_ranges = zeros(1, Simulation_time);
bistatic_doppler_shifts = zeros(1, Simulation_time);

prev_total_range =0;
num_waypoints = size(Target_waypoints, 1);


for waypoint_idx = 2:num_waypoints 
    prev_waypoint = Target_waypoints(waypoint_idx - 1, :);
    current_waypoint = Target_waypoints(waypoint_idx, :);
    
    delta_position = current_waypoint - prev_waypoint;
    delta_time = Time_for_waypoints(waypoint_idx)-Time_for_waypoints(waypoint_idx-1);
    
    speed_xyz = delta_position / delta_time;
    speed_magnitude = norm(delta_position) / delta_time;

    Target_Pos_1 = prev_waypoint;
    num_steps = Time_for_waypoints(waypoint_idx);
    start_step = Time_for_waypoints(waypoint_idx-1);
    
    % Estimate the bistatic ranges at each time step
    for step = start_step+1:num_steps+1
        estimated_position = Target_Pos_1 +   speed_xyz* UpdateInterval * (step-1);
        
        % Calculate the distance between target and transmitter
        range_tx_target = norm(estimated_position - Tx_Pos);
        
        % Calculate the distance between target and reference receiver
        range_ref_target = norm(estimated_position - RefRx_Pos);
        
        % Calculate the distance between target and surveillance receiver
        
        total_range = range_tx_target + range_ref_target;
    
        % Calculate the bistatic ranges
        bistatic_range_ref = range_tx_target + range_ref_target - Baseline;
    
        % Store the bistatic ranges in the matrix
        bistatic_ranges(step) = bistatic_range_ref;
    
        %Calculate the rate of change of the total range with respect to time
        
        d_total_range_dt = (total_range - prev_total_range) / UpdateInterval;
        
        prev_total_range = total_range;
        
        % Calculate the Doppler shift
        doppler_shift = - (1 / wavelength) * (d_total_range_dt);
        
        bistatic_doppler_shifts(step)=doppler_shift;
    
         % Calculate the Doppler shift for the previous position
        if step == 1
            prev_position = Target_Pos_1 + speed_xyz * UpdateInterval * (step - 2);
            prev_range_tx_target = norm(prev_position - Tx_Pos);
            prev_range_ref_target = norm(prev_position - RefRx_Pos);
            prev_total_range = prev_range_tx_target + prev_range_ref_target;
            
            % Calculate the Doppler shift for the previous position
            prev_doppler_shift = - (1 / wavelength) * ((total_range - prev_total_range) / UpdateInterval);
            
            % Store the previous Doppler shift in the matrix
            bistatic_doppler_shifts(step) = prev_doppler_shift;
            prev_total_range = total_range;
        end
    
    end
end


if exist('true_data.h5', 'file')
    delete('true_data.h5');
end

disp(bistatic_ranges);
disp(bistatic_doppler_shifts);
% Save bistatic ranges and Doppler shifts to an HDF5 file
h5create('true_data.h5', '/bistatic_ranges', size(bistatic_ranges));
h5write('true_data.h5', '/bistatic_ranges', bistatic_ranges);

h5create('true_data.h5', '/doppler_shifts', size(bistatic_doppler_shifts));
h5write('true_data.h5', '/doppler_shifts', bistatic_doppler_shifts);
