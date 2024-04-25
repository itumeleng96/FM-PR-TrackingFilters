%Lane Change Maneuver Scenario 
TargetPos = [[3806;20680;3200;],[-3116;14010;3200],[-6884;10971;4000],[-14979;3765;4000]];
TargetWayPoints =[0,20,30,60];

%Landing Maneuver
%TargetPos = [[4000;18000;2000;],[-1000;8000;1700],[-1500;6000;1500],[-1000;4000;1300],[0;4000;1100],[500;6000;1000]];
%TargetWayPoints =[0,30,38,45,50,60];

%360 Maneuver
%TargetPos = [[4000;18000;2000;],[-8528;13851;2000],[-12354;15339;2000],[-7997;17783;2000],[-13311;11938;2000],[-20431;9706;2000]];
%TargetWayPoints =[0,35,50,75,95,120];


SimulationTime = TargetWayPoints(end)-TargetWayPoints(1)+1;

% Perform cubic interpolation to estimate target positions between waypoints

%coeffs = spline(TargetWayPoints, TargetPos);
%interpolated_positions = ppval(coeffs, time_intervals);

interp_positions =[];
interpolated_positions = [];
LastPosition  =[];
LastPositionIndex =[];
% Loop through each waypoint section


% Define the update interval
UpdateInterval = 1;

% Calculate the baseline

% Initialize arrays to store bistatic ranges and Doppler shifts
bistatic_ranges = [];
bistatic_doppler_shifts = [];

% Loop through each time step
prev_total_range =0;

d_start = [500; 100; 100]; % Example slope at the start (tweak as needed)
d_end = [-500; -500; -1000]; % Example slope at the end (tweak as needed)

% Create cubic spline with clamped boundary conditions
% The derivative values are combined with the original data
pp_x = spline(TargetWayPoints, [d_start(1), TargetPos(1, :), d_end(1)]);
pp_y = spline(TargetWayPoints, [d_start(2), TargetPos(2, :), d_end(2)]);
pp_z = spline(TargetWayPoints, [d_start(3), TargetPos(3, :), d_end(3)]);

% Define the new time-based parameterization for interpolation (1-second intervals)
t_interp = 0:1:60;

% Interpolate using cubic splines with clamped boundary conditions
x_interp = ppval(pp_x, t_interp);
y_interp = ppval(pp_y, t_interp);
z_interp = ppval(pp_z, t_interp);

interpolated_posx =[x_interp;y_interp;z_interp];

%code for interpolated_posx

figure(1)
plot3(interpolated_posx(1, :), interpolated_posx(2, :), interpolated_posx(3, :), 'b.-');
grid on;