% This script calculates the Bistatic range(m)  and Doppler(Hz) From
% Bistatic setup
% Simulation
% Author: Itumeleng Malemela


%The Positions are all in the format Below
%[x,y,z] where x,y,z are the coordinates in meters in the Cartesian Coordinate System

% Single Target Scenario: Two points with a 60-unit time interval
%TargetPos = [
%    4000, 18000, 3600;  % Start position (x, y, z)
%    -2000, 3000, 1600;  % End position
%];
%TargetWayPoints = [0, 60];  % Time intervals for waypoints (from 0 to 60 units)

% Lane Change Maneuver Scenario: Four points with varying time intervals
%{
TargetPos = [
    3806, 20680, 3200;  % First point (x, y, z)
    -3116, 14010, 3200; % Second point (lane change)
    -6884, 10971, 4000; % Third point
    -14979, 3765, 4000; % Fourth point (end of the lane change)
];
TargetWayPoints = [0, 20, 30, 60];  % Waypoints indicating time intervals
%}
% Landing Maneuver: Six points simulating a landing trajectory
%{
TargetPos = [
    4000, 18000, 2000;  % Initial point (x, y, z)
    -1000, 8000, 1700;  % Descent begins
    -1500, 6000, 1500;  % Further descent
    -1000, 4000, 1300;  % Approaching the landing area
    0, 4000, 1100;      % Close to landing
    500, 6000, 1000;    % Landing completed
];
TargetWayPoints = [0, 30, 38, 45, 50, 60]; 
%}
% 360 Maneuver: Six points representing a full circle maneuver
%{
TargetPos = [
    4000, 18000, 2000;  % Start point (x, y, z)
    -8528, 13851, 2000; % First quarter
     -12354, 15339, 2000; % Halfway through
     -7997, 17783, 2000; % Three quarters
     -13311, 11938, 2000; % Near completion
     -20431, 9706, 2000;  % Full circle
];
TargetWayPoints = [0, 35, 50, 75, 95, 120];  % Time intervals for each segment of the 360 maneuver
%}
% Multi-Target Scenario - Target 1: Two points representing a different trajectory
%TargetPos = [
%    4000, 18000, 3600;  % Start point (x, y, z)
%    -2000, 3000, 1600;  % End point
%];
%TargetWayPoints = [0, 60];  % Time interval from start to end
% Multi-Target Scenario - Target 2: Another two-point trajectory
% This represents an additional target with a distinct movement pattern

TargetPos = [
    4000, 4000, 1600;  % Start point (x, y, z)
    2000, 20000, 3600; % End point
];
TargetWayPoints = [0, 60];  % Time interval from start to end

wavelength = 299792458/94e6;
c=299792458;

%Transmitter Position
Tx_Pos=[6440;10760;1056];
%Reference Reciever Position 
RefRx_Pos=[0;0;1000];

%Surveillance Reciever Position 
SurvRx_Pos=[0;1;1000];

%Target Positions

SimulationTime = TargetWayPoints(end)-TargetWayPoints(1)+1;

% Define the update interval
UpdateInterval = 1;

% Calculate the baseline
Baseline = norm(RefRx_Pos - Tx_Pos);

% Initialize arrays to store bistatic ranges and Doppler shifts
bistatic_ranges = [];
bistatic_doppler_shifts = [];


% Create interpolation time
t_interp = linspace(TargetWayPoints(1), TargetWayPoints(end), SimulationTime);


% Create cubic spline interpolations with natural boundary conditions
spline_x = csape(TargetWayPoints, TargetPos(:, 1), 'variational');
spline_y = csape(TargetWayPoints, TargetPos(:, 2), 'variational');
spline_z = csape(TargetWayPoints, TargetPos(:, 3), 'variational');

% Evaluate the splines at the interpolation points
x_interp = fnval(spline_x, t_interp);
y_interp = fnval(spline_y, t_interp);
z_interp = fnval(spline_z, t_interp);

% Interpolated positions (as a 3xN matrix)
interpolated_posx = [x_interp; y_interp; z_interp];


%code for interpolated_posx

figure(1)
plot3(interpolated_posx(1, :), interpolated_posx(2, :), interpolated_posx(3, :), 'b.-');
title('Target Path');
xlabel('x');
ylabel('y');
zlabel('z');

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