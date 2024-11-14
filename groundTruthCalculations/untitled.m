% Define the 3D line
x_line = linspace(-10, 10, 100);
y_line = sin(x_line);
z_line = cos(x_line);

% Plot the 3D line
figure;
plot3(x_line, y_line, z_line, 'r-', 'LineWidth', 2);
hold on;

% Define a plane
[x_plane, y_plane] = meshgrid(-10:1:10, -10:1:10);
z_plane = 3 * ones(size(x_plane)); % A horizontal plane at z = 3

% Plot the plane
mesh(x_plane, y_plane, z_plane, 'FaceAlpha', 0.5); % Mesh plot with transparency

% Labels and view adjustments
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('3D Line with a Plane');
grid on;
view(3); % Set to 3D view
hold off;