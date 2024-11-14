% Load and visualize an STL file

% Specify the STL file path
stlFilePath = './plane/plane3d.stl';

% Read the STL file
airplane_model = stlread(stlFilePath);

% Extract vertices and faces
vertices = airplane_model.Points;
faces = airplane_model.ConnectivityList;

% Plot the 3D model
figure;
patch('Vertices', vertices, 'Faces', faces, ...
      'FaceColor', [0.8 0.8 1.0], 'EdgeColor', 'none', 'FaceLighting', 'gouraud');

% Setup lighting and view
camlight('headlight');
lighting gouraud;
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
title('3D Model Visualization');
