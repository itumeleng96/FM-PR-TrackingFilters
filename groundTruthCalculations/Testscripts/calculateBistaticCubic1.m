% This script calculates the Bistatic range(m)  and Doppler(Hz) From
% Bistatic setup
% Simulation
% Author: Itumeleng Malemela

%The Positions are all in the format Below
%[x,y,z] where x,y,z are the coordinates in meters in the Cartesian Coordinate System

wavelength = 299792458/94e6;
c=299792458;

%Transmitter Position
Tx_Pos=[6440;10760;1056];
%Reference Reciever Position 
RefRx_Pos=[0;0;1000];

%Surveillance Reciever Position 
SurvRx_Pos=[0;1;1000];

%Target Positions

%Single Target Scenario 
%TargetPos = [[4000;18000;3600;],[-2000;3000;1600]];
%TargetWayPoints =[0,60];

%Lane Change Maneuver Scenario 
%TargetPos = [[3806;20680;3200;],[-3116;14010;3200],[-6884;10971;4000],[-14979;3765;4000]];
%TargetWayPoints =[0,20,30,60];

%Landing Maneuver
%TargetPos = [[4000;18000;2000;],[-1000;8000;1700],[-1500;6000;1500],[-1000;4000;1300],[0;4000;1100],[500;6000;1000]];
%TargetWayPoints =[0,30,38,45,50,60];

%360 Maneuver
TargetPos = [[4000;18000;2000;],[-8528;13851;2000],[-12354;15339;2000],[-7997;17783;2000],[-13311;11938;2000],[-20431;9706;2000]];
TargetWayPoints =[0,35,50,75,95,120];

%Multit-Target Target 1 
%TargetPos = [[4000;18000;3600;],[-2000;3000;1600]];
%TargetWayPoints =[0,60];

%Multit-Target Target 2
%TargetPos = [[4000;4000;1600;],[2000;20000;3600]];
%TargetWayPoints =[0,60];

SimulationTime = TargetWayPoints(end)-TargetWayPoints(1)+1;

% Define the update interval
UpdateInterval = 1;

% Calculate the baseline
Baseline = norm(RefRx_Pos - Tx_Pos);

% Initialize arrays to store bistatic ranges and Doppler shifts
bistatic_ranges = [];
bistatic_doppler_shifts = [];
t_interp = 0:1:SimulationTime;

x_interp = [];
y_interp = [];
z_interp = [];

% Loop through each waypoint section except the last one
for i = 1:length(TargetWayPoints) - 1
    % Get the waypoints for the current segment
    waypoints = TargetWayPoints(i:i+1);
    
    % Get the position data for the current segment
    x_data = TargetPos(1, waypoints == TargetWayPoints(i:i+1));
    y_data = TargetPos(2, waypoints == TargetWayPoints(i:i+1));
    z_data = TargetPos(3, waypoints == TargetWayPoints(i:i+1));
    
    % Define the unique segment time range for interpolation
    if i < length(TargetWayPoints) - 1
        t_segment = t_interp(t_interp >= TargetWayPoints(i) & t_interp < TargetWayPoints(i+1));
    else
        t_segment = t_interp(t_interp >= TargetWayPoints(i) & t_interp <= TargetWayPoints(i+1));
    end
    
    % Perform cubic interpolation for each position
    x_interp_segment = interp1(waypoints, x_data, t_segment, 'spline');
    y_interp_segment = interp1(waypoints, y_data, t_segment, 'spline');
    z_interp_segment = interp1(waypoints, z_data, t_segment, 'spline');
    
    
    % Concatenate with the complete interpolated position
    x_interp = [x_interp, x_interp_segment];
    y_interp = [y_interp, y_interp_segment];
    z_interp = [z_interp, z_interp_segment];
end

% Final interpolated positions
interpolated_posx = [x_interp; y_interp; z_interp];


%code for interpolated_posx

figure(1)
plot3(interpolated_posx(1, :), interpolated_posx(2, :), interpolated_posx(3, :), 'b.-');
grid on;

for position_index = 1:(size(interpolated_posx, 2) - 1)
    % Delta-Time: time taken in position section
    
    % Current position
    Target_Pos_1 = interpolated_posx(:, position_index);
    % Next position
    Target_Pos_end = interpolated_posx(:, position_index + 1);
    
    % Ranges from Tx and RefRx to the target
    range_tx_target = norm(Target_Pos_1 - Tx_Pos);
    range_ref_target = norm(Target_Pos_1 - RefRx_Pos);
    
    % Ranges for the end position
    range_tx_target_end = norm(Target_Pos_end - Tx_Pos);
    range_ref_target_end = norm(Target_Pos_end - RefRx_Pos);

    % Calculate velocities (rate of change of range)
    V_t = (range_tx_target-range_tx_target_end) / UpdateInterval;
    V_r = (range_ref_target-range_ref_target_end) / UpdateInterval;

    % Calculate bistatic Doppler shift using the correct formula
    doppler_shift = (94e6 / c) * (V_t + V_r);  % Assuming similar velocities as the baseline formula
    bistatic_doppler_shifts(end + 1) = doppler_shift;

    % Calculate the bistatic range
    bistatic_range = range_tx_target + range_ref_target - Baseline;

    % Store the bistatic range
    bistatic_ranges(end + 1) = bistatic_range;
end

figure(3);
plot(bistatic_ranges,bistatic_doppler_shifts);
title('Bistatic range vs Doppler shift (Ground Truth)');
xlabel('Bistatic range(m)');
ylabel('Bistatic Doppler shift(Hz)');
xlim([0 7e4]);
ylim([-200 200]);


if exist('./true_data.h5', 'file')
    delete('./true_data.h5');
end

% Save bistatic ranges and Doppler shifts to an HDF5 file
h5create('./true_data.h5', '/bistatic_ranges', size(bistatic_ranges));
h5write('./true_data.h5', '/bistatic_ranges', bistatic_ranges);

h5create('./true_data.h5', '/doppler_shifts', size(bistatic_doppler_shifts));
h5write('./true_data.h5', '/doppler_shifts', bistatic_doppler_shifts);