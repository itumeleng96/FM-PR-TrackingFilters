clear
% Define the points and their names
x = [0, 0, 6440, 10000, -2000];
y = [1, 0, 10760, 12000,-2000];
z = [1000,1000, 1000, 3600, 1600];

% Plot the points with their names
plot3(x, y, z, 'o')

% Add labels to the first two points

% Connect the last two points with a line
hold on
line([x(4), x(5)], [y(4), y(5)], [z(4), z(5)])


text(x(1), y(1), z(1), 'Ref and Surv Rx') 
text(x(3), y(3), z(3), 'Tx') 
text(x(4), y(4), z(4), 'Target1 Position 1') 
text(x(5), y(5), z(5), 'Target1 Position 2') 


% Add labels and grid
grid on
xlabel('x position(m)')
ylabel('y position(m)')
zlabel('altitude (m)')

