% Example waypoints
%{ 
Lane Change
x = [0, 6440, 4000, 2000, 1000, -2000];
y = [0, 10760, 18000, 14000, 10000, 2000];
z = [1000, 1056, 3600, 3600, 3300, 3300];
%}

%Landing Maneuver
x = [0, 6440, 4000,-1000,-1500,-1000,0,500];
y = [0, 10760, 18000, 8000,6000,4000,4000,6000];
z = [1000, 1056, 2000, 1700,1500,1300,1100,1000];

%{

x = [0, 6440, 500,0,-1000,-1500,-1000,4000];
y = [0, 10760, 6000, 4000,4000,6000,8000,18000];
z = [1000, 1056, 1000, 1100,1300,1500,1700,2000];
%}
% Plot antennas as discrete points
figure;
plot3(x(1), y(1), z(1), 'o', 'LineWidth', 2, 'MarkerSize', 10, 'Color', 'blue', 'DisplayName', 'Ref and Surv Rx');
hold on;
plot3(x(2), y(2), z(2), 'o', 'LineWidth', 2, 'MarkerSize', 10, 'Color', 'red', 'DisplayName', 'Transmitter');

% Plot the remaining waypoints connected by a line
plot3(x(3:end), y(3:end), z(3:end), '-o', 'LineWidth', 2, 'Color', 'green', 'DisplayName', 'Target Motion');

grid on;
xlabel('X Position (m)');
ylabel('Y Position (m)');
zlabel('Altitude (m)');
title('Target Motion Based on Waypoints');
legend('Location', 'Best');

% Adjust viewpoint for a better 3D perspective
view(3);
