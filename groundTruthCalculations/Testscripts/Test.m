% Sample code for plotting XYZ points in 3D
Coords = [
    3806 -3116 -6884 -14979;
    20680 14010 10971 3765;
    3200 3200 4000 4000;
];

CoordsTime = [0, 20, 30, 60];

% Calculate second derivatives
dd = finalizeCubic(Coords, CoordsTime);

% Plot XYZ points in 3D for each time point
figure;
test =[];
for t = 0:60
    % Get interpolated coordinate at time t
    coord = getPositionCubic(t, Coords, CoordsTime, dd);
    test = [test ,coord];
    % Plot XYZ points
end

figure(1)
plot3(test(1, :), test(2, :), test(3, :), 'b.-');
xlabel('X');
ylabel('Y');
zlabel('Z');
title(['3D Plot of XYZ Points at time t = ' num2str(t)]);
grid on;