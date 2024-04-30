% Load the XML file
xmlFile = 'FERS/flightScenarios/scenario_1_laneChange.fersxml'; % replace with your file path
xmlData = fileread(xmlFile);

% Find position waypoints using regular expressions
pattern = '<positionwaypoint>\s*<x>(.*?)</x>\s*<y>(.*?)</y>\s*<altitude>(.*?)</altitude>\s*<time>(.*?)</time>\s*</positionwaypoint>';
matches = regexp(xmlData, pattern, 'tokens');

% Extract data from matches
x = cellfun(@(match) str2double(match{1}), matches);
y = cellfun(@(match) str2double(match{2}), matches);
z = cellfun(@(match) str2double(match{3}), matches);

% Display the extracted data
disp('Extracted Position Waypoints:');
disp('x =');
disp(x);
disp('y =');
disp(y);
disp('z =');
disp(z);
