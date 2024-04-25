TargetPos = [[3806;20680;3200;],[-3116;14010;3200],[-6884;10971;4000],[-14979;3765;4000]];
TargetWayPoints =[0,20,30,60];


% Perform cubic spline interpolation on the entire set of coordinates
interp_positions = interp1([0, TargetWayPoints], [TargetPos(:, 1), TargetPos], time_intervals, 'spline')';

% Plot interpolated positions in 3D
figure;
plot3(interp_positions(:, 1), interp_positions(:, 2), interp_positions(:, 3), 'b.-');
hold on;
plot3(TargetPos(1, :), TargetPos(2, :), TargetPos(3, :), 'ro', 'MarkerSize', 10);  % Plot endpoints
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Interpolated Positions');
grid on;