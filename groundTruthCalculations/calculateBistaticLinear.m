% This script calculates the Bistatic range(m)  and Doppler(Hz) From
% Bistatic setup
% Simulation
% Author: Itumeleng Malemela

%The Positions are all in the format Below
%[x,y,z] where x,y,z are the coordinates in meters in the Cartesian Coordinate System

wavelength = 299792458/94e6;

%Transmitter Position
Tx_Pos=[6440;10760;1056];
%Reference Reciever Position 
RefRx_Pos=[0;0;1000];

%Surveillance Reciever Position 
SurvRx_Pos=[0;1;1000];

%Target Positions

%Single Target Scenario 
TargetPos = [[4000;18000;3600;],[-2000;3000;1600]];
TargetWayPoints =[0,60];

%Lane Change Maneuver Scenario 
%TargetPos = [[3806;20680;3200;],[-3116;14010;3200],[-6884;10971;4000],[-14979;3765;4000]];
%TargetWayPoints =[0,20,30,60];

%Landing Maneuver
%TargetPos = [[4000;18000;2000;],[-1000;8000;1700],[-1500;6000;1500],[-1000;4000;1300],[0;4000;1100],[500;6000;1000]];
%TargetWayPoints =[0,30,38,45,50,60];

%360 Maneuver
%TargetPos = [[4000;18000;2000;],[-8528;13851;2000],[-12354;15339;2000],[-7997;17783;2000],[-13311;11938;2000],[-20431;9706;2000]];
%TargetWayPoints =[0,35,50,75,95,120];

%Multit-Target Target 1 
%TargetPos = [[4000;18000;3600;],[-2000;3000;1600]];
%TargetWayPoints =[0,60];

%Multit-Target Target 2
%TargetPos = [[4000;4000;1600;],[2000;20000;3600]];
%TargetWayPoints =[0,60];

UpdateInterval = 1;

Baseline = norm(RefRx_Pos - Tx_Pos);

SimulationTime = TargetWayPoints(end)-TargetWayPoints(1)+1;

bistatic_ranges = [];
bistatic_doppler_shifts = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For every section of waypoints we need
% to calculate the bistatic values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prev_total_range =0;

%Loop through all the waypoints
for waypoint=1:size(TargetWayPoints,2)-1
    
    %Delta-Time : time taken in waypoint section
    delta_time = TargetWayPoints(waypoint+1)-TargetWayPoints(waypoint);

    %Delta-Position : Distance between waypoint section
    delta_postion = TargetPos(:,waypoint+1)-TargetPos(:,waypoint);

    speed_xyz = delta_postion/delta_time;

    %num steps for waypoint
    num_steps = delta_time / UpdateInterval;
    

    Target_Pos_1 = TargetPos(:,waypoint);
    %////////////////////////////////////////////////////////////////////////////////
    % CALCULATE BISTATIC RANGE AND DOPPLER FOR WAYPOINT SECTION
    %///////////////////////////////////////////////////////////////////////////////

    % Estimate the bistatic ranges at each time step
    for step = 1:num_steps
        estimated_position = Target_Pos_1 +   speed_xyz* UpdateInterval * step;
    
        
        % Calculate the distance between target and transmitter
        range_tx_target = norm(estimated_position - Tx_Pos);
        
        % Calculate the distance between target and reference receiver
        range_ref_target = norm(estimated_position - RefRx_Pos);
        
        % Calculate the distance between target and surveillance receiver
        range_surv_target = norm(estimated_position - SurvRx_Pos);
        
        total_range = range_tx_target + range_ref_target;
    
        % Calculate the bistatic ranges
        bistatic_range_ref = range_tx_target + range_ref_target - Baseline;
        bistatic_range_surv = range_tx_target + range_surv_target - Baseline;
        
        % Store the bistatic ranges in the matrix
        bistatic_ranges(end+1) = bistatic_range_ref;
    
        %Calculate the rate of change of the total range with respect to time
        
        d_total_range_dt = (total_range - prev_total_range) / UpdateInterval;
        
        prev_total_range = total_range;
        
        % Calculate the Doppler shift
        doppler_shift = - (1 / wavelength) * (d_total_range_dt);
        
        bistatic_doppler_shifts(end+1)=doppler_shift;
    
         % Calculate the Doppler shift for the previous position
        if step == 1 && waypoint==1
            prev_position = Target_Pos_1 + speed_xyz * UpdateInterval * (step - 1);
            prev_range_tx_target = norm(prev_position - Tx_Pos);
            prev_range_ref_target = norm(prev_position - RefRx_Pos);
            prev_total_range = prev_range_tx_target + prev_range_ref_target;
            
            % Calculate the Doppler shift for the previous position
            prev_doppler_shift = - (1 / wavelength) * ((total_range - prev_total_range) / UpdateInterval);
            
            % Store the previous Doppler shift in the matrix
            bistatic_doppler_shifts(step) = prev_doppler_shift;
            prev_total_range = total_range;
        end
    
    end
end



%figure(1);
%plot(bistatic_doppler_shifts,  'b-');

%figure(2);
%plot(bistatic_ranges,  'b-');

figure(3);
plot(bistatic_ranges,bistatic_doppler_shifts);
title('Bistatic range vs Doppler shift (Ground Truth)');
xlabel('Bistatic range(m)');
ylabel('Bistatic Doppler shift(Hz)');
xlim([0 7e4]);
ylim([-200 200]);

if exist('../true_data.h5', 'file')
    delete('../true_data.h5');
end

%disp(bistatic_ranges);
%disp(bistatic_doppler_shifts);

% Save bistatic ranges and Doppler shifts to an HDF5 file
h5create('true_data.h5', '/bistatic_ranges', size(bistatic_ranges));
h5write('true_data.h5', '/bistatic_ranges', bistatic_ranges);

h5create('true_data.h5', '/doppler_shifts', size(bistatic_doppler_shifts));
h5write('true_data.h5', '/doppler_shifts', bistatic_doppler_shifts);
