clear
% Define the points and their names
x = [0, 6440, 4000, -4000,50000];
y = [0, 10760, 18000,2000,50000];
z = [1000,1056, 3600, 1600,1600];

% Plot the points with their names
plot3(x, y, z, 'o');

% Add labels to the first two points

% Connect the last two points with a line
hold on
%line([x(3), x(4)], [y(3), y(4)], [z(3), z(4)]);
%line([x(5), x(6)], [y(5), y(6)], [z(5), z(6)]);
line([x(3), x(4)], [y(3), y(4)], [z(3), z(4)]);


text(x(1), y(1), z(1), 'Ref and Surv Rx') ;
text(x(2), y(2), z(2), 'Tx') ;
text(x(3), y(3), z(3), 'Target1 Position 1') ;
text(x(4), y(4), z(4), 'Target1 Position 2') ;
text(x(5), y(5), z(5), 'Noise transmitter') ;




% Add labels and grid
grid on;
xlabel('x position(m)');
ylabel('y position(m)');
zlabel('altitude (m)');
