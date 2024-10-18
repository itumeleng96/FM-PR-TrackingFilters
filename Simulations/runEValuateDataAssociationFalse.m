clc;
clear;
close all;

addpath('FERS/', ...
        'cfar/', ...
        'meanShiftCluster/', ...
        'multiTargetTracking/', ...
        'DPI_Suppression', ...
        'TrackingFilter-KalmanFilter/', ...
        'TrackingFilter-CSKF/', ...
        'TrackingFilter-ParticleFilter/', ...
        'TrackingFilter-CS-ParticleFilter/', ...
        'TrackingFilter-UKF/', ... 
        'TrackingFilter-CSUKF/', ... 
        'TrackingFilter-RGNF/',...
        'TrackingFilter-CSRGNF/');

%FLIGHT Scenarios
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/flightScenarios/scenario_1_laneChange.fersxml');
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/flightScenarios/scenario_2_landingManeuver.fersxml');
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/flightScenarios/scenario_3_takeoffManeuver.fersxml');
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/flightScenarios/scenario_4_360.fersxml');
system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/flightScenarios/scenario_5_2_targets.fersxml');

%Noise Scenarios
%SCENARIO 3
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/NoiseScenarios/scenario_1_fm_noise.fersxml');

%SCENARIO 2
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/NoiseScenarios/scenario_2_white_noise.fersxml');

%SCENARIO 1 
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/BackupScenarios/scenario_1_singleFile.fersxml'); 


%%%%%
%%% MUST decrease number for track deletion when no update
%%%%%

% h5 Import from FERS simulation
[Ino, Qno, scale_no] = loadfersHDF5('direct.h5');

[Imov, Qmov, scale_mov] = loadfersHDF5('echo.h5');


I_Qmov = Imov + 1i*Qmov;
I_Qmov = I_Qmov.*scale_mov;
I_Qno = Ino + 1i*Qno;
I_Qno = I_Qno.*scale_no;

%I_Qmov=I_Qmov-I_Qno;

fs = 200000;
dopp_bins = 200;
delay = 233e-6;
c=299792458;
range_delay = delay*c;

%DPI Cancellation
proc = struct('cancellationMaxRange_m', range_delay, ...
              'cancellationMaxDoppler_Hz', 4, ...
              'TxToRefRxDistance_m', 12540, ...
              'nSegments', 1, ...
              'nIterations', 20, ...
              'Fs', fs, ...
              'alpha', 0, ...
              'initialAlpha', 0);



s1 = I_Qmov;   %Surv
s2 = I_Qno;    %Ref

initial=1;

current=fs;                                %based on samples in transmitted signal
simulation_time = size(I_Qmov,1)/fs;       %Simulation time: number of data points/sampling frequency


ard = [];
rdm =[];


f3=figure();
f3.Position = [4000 10 1000 800]; 
movegui(f3,'southwest');

%f5=figure();
%f5.Position = [4000 10 1000 800]; 
%movegui(f5,'southeast');

%Create MTT object
confirmationThreshold=4;
deletionThreshold=4;
gatingThreshold=20; %scalar value for ellipsoidal gate

%FilterType 1: Kalman Filter
%FilterType 2: Covariance Scaling Kalman Filter
%FilterType 3: Particle Filter
%FilterType 4: CS Particle Filter
%FilterType 5: UKF  Filter
%FilterType 6: Covariance Scaling UKF  Filter
%FilterType 7: RGNF Filter
%FilterType 8: Covariance scaling RGNF Filter


multiTargetTracker1 = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,2);
multiTargetTracker2 = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,4);
multiTargetTracker3 = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,6);
multiTargetTracker4 = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,8);

%LOG_LIKELIHOODS
doppler_ll_1=[];
range_ll_1=[];

doppler_ll_2=[];
range_ll_2=[];


doppler_ll_3=[];
range_ll_3=[];


doppler_ll_4=[];
range_ll_4=[];

prevCentroids=[];

%True Data for single Target Scenario
rangeTrueData = h5read('./groundTruthCalculations/true_data.h5', '/bistatic_ranges');
dopplerTrueData = h5read('./groundTruthCalculations/true_data.h5', '/doppler_shifts');

rangeTrueData1 = h5read('./groundTruthCalculations/true_data_1.h5', '/bistatic_ranges');
dopplerTrueData1 = h5read('./groundTruthCalculations/true_data_1.h5', '/doppler_shifts');

trueData1 = [rangeTrueData1;dopplerTrueData1];
trueData = [rangeTrueData;dopplerTrueData];

%rangeTrueData = h5read('./measurement_data.h5', '/bistatic_ranges');
%dopplerTrueData = h5read('./measurement_data.h5', '/doppler_shifts');
numSimulations = 1; 
simulation_time  = 10;
% Initialize counters for confirmations

% Define PFA values to iterate over
pfa_values = [10^-5, 10^-6, 10^-7, 10^-8, 10^-9];

% Initialize containers to store results
avgKalmanResults = zeros(size(pfa_values));
avgUKFResults = zeros(size(pfa_values));
avgParticleResults = zeros(size(pfa_values));
avgRGNFResults = zeros(size(pfa_values));

for pfa_index = 1:length(pfa_values)
    pfa = pfa_values(pfa_index);
    fprintf('\nRunning simulation with PFA = %.1e\n', pfa);

    % Reset counters for confirmations
    totalConfirmedFalseKalman = 0;
    totalConfirmedFalseUKF = 0;
    totalConfirmedFalseParticle = 0;
    totalConfirmedFalseRGNF = 0;

    multiTargetTracker1 = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,2);
    multiTargetTracker2 = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,4);
    multiTargetTracker3 = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,6);
    multiTargetTracker4 = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,8);


    for j = 1:numSimulations
        initial = 1;
        current = fs;

        for i = 1:simulation_time
            s1 = I_Qmov(initial:current); % Surveillance
            s2 = I_Qno(initial:current);  % Reference

            s1 = procECA_Optimized(s2, s1, proc);

            % Range-Doppler Map
            [y, ard_] = ardNoPlot(s1, s2, fs, dopp_bins, delay, i, ard);

            % Use PFA in CFAR detection
            [targetClusters, RDM, rdm_] = ca_cfar(y.', pfa, fs, dopp_bins, delay, 40);

            % Extract cluster centroids using MeanShift
            [clusterCentroids, ~, ~, ~] = meanShift(targetClusters, 10, 8);

            % Update trackers
            multiTargetTracker1 = multiTargetTracker1.createNewTracks(clusterCentroids, i);
            multiTargetTracker2 = multiTargetTracker2.createNewTracks(clusterCentroids, i);
            multiTargetTracker3 = multiTargetTracker3.createNewTracks(clusterCentroids, i);
            multiTargetTracker4 = multiTargetTracker4.createNewTracks(clusterCentroids, i);

            % Maintain and predict tracks
            multiTargetTracker1 = multiTargetTracker1.maintainTracks().predictionStage();
            multiTargetTracker2 = multiTargetTracker2.maintainTracks().predictionStage();
            multiTargetTracker3 = multiTargetTracker3.maintainTracks().predictionStage();
            multiTargetTracker4 = multiTargetTracker4.maintainTracks().predictionStage();

            % Plot results
            multiTargetTracker2 = multiTargetTracker2.plotMultiTargetTracking(fs, dopp_bins, delay, i, f3, RDM);

            % Update stage with measurements
            multiTargetTracker1 = multiTargetTracker1.updateStage(clusterCentroids, i);
            multiTargetTracker2 = multiTargetTracker2.updateStage(clusterCentroids, i);
            multiTargetTracker3 = multiTargetTracker3.updateStage(clusterCentroids, i);
            multiTargetTracker4 = multiTargetTracker4.updateStage(clusterCentroids, i);

            % Update indices for the next iteration
            initial = current + 1;
            current = current + fs;
        end

        % Compare confirmed tracks
        totalConfirmedFalseKalman = totalConfirmedFalseKalman + compareTracks(multiTargetTracker1.getConfirmedTracks, {trueData, trueData1}, 0.5, 1, 10, 0.1);
        totalConfirmedFalseUKF = totalConfirmedFalseUKF + compareTracks(multiTargetTracker2.getConfirmedTracks, {trueData, trueData1}, 0.5, 1, 10, 0.1);
        totalConfirmedFalseParticle = totalConfirmedFalseParticle + compareTracks(multiTargetTracker3.getConfirmedTracks, {trueData, trueData1}, 0.5, 1, 10, 0.1);
        totalConfirmedFalseRGNF = totalConfirmedFalseRGNF + compareTracks(multiTargetTracker4.getConfirmedTracks, {trueData, trueData1}, 0.5, 1, 10, 0.1);    
    end

    % Calculate averages for the current PFA value
    avgKalmanResults(pfa_index) = totalConfirmedFalseKalman / numSimulations;
    avgUKFResults(pfa_index) = totalConfirmedFalseUKF / numSimulations;
    avgParticleResults(pfa_index) = totalConfirmedFalseParticle / numSimulations;
    avgRGNFResults(pfa_index) = totalConfirmedFalseRGNF / numSimulations;
end

% Display results
for pfa_index = 1:length(pfa_values)
    fprintf('\nResults for PFA = %.1e:\n', pfa_values(pfa_index));
    fprintf('  Average Confirmed Kalman False Tracks: %.2f\n', avgKalmanResults(pfa_index));
    fprintf('  Average Confirmed UKF False Tracks: %.2f\n', avgUKFResults(pfa_index));
    fprintf('  Average Confirmed Particle False Tracks: %.2f\n', avgParticleResults(pfa_index));
    fprintf('  Average Confirmed RGNF False Tracks: %.2f\n', avgRGNFResults(pfa_index));
end
% Update this function to count confirmed false tracks
function numConfirmedFalseTracks = compareTracks(confirmedTracks, trueTracks, rangeWeight, dopplerWeight, timeWindow, costOfNonAssignment)
    % Inputs:
    % - confirmedTracks: Cell array of confirmed tracks (each 2 x N array)
    % - trueTracks: Cell array of ground truth tracks (each 2 x T array)
    % - rangeWeight: Weight for range similarity in cost calculation
    % - dopplerWeight: Weight for Doppler similarity in cost calculation
    % - timeWindow: Number of time samples to consider for matching
    % - costOfNonAssignment: Scalar penalty for non-assignment

    % Ensure all true tracks are truncated to the same time window
    nConfirmed = length(confirmedTracks);
    nTrue = length(trueTracks);

    % Handle case where there are no confirmed or true tracks
    if nConfirmed == 0 || nTrue == 0
        numConfirmedFalseTracks = nConfirmed; % all confirmed are false if no true tracks
        return;
    end

    % Initialize the cost matrix (rows: confirmed tracks, columns: true tracks)
    costMatrix = zeros(nConfirmed, nTrue);

    % Fill the cost matrix with weighted similarity costs
    for i = 1:nConfirmed
        confirmedTrack = confirmedTracks{i};

        for j = 1:nTrue
            trueTrack = trueTracks{j};

            % Extract and truncate range and Doppler data for both tracks
            confirmedRange = confirmedTrack(1,:);
            confirmedDoppler = confirmedTrack(2,:);
            trueRange = trueTrack(1, 1:timeWindow);
            trueDoppler = trueTrack(2, 1:timeWindow);

            % Calculate normalized cross-correlation
            rangeCost = 1 - max(abs(xcorr(confirmedRange, trueRange, 'none'))) / timeWindow;
            dopplerCost = 1 - max(abs(xcorr(confirmedDoppler, trueDoppler, 'none'))) / timeWindow;

            % Weighted sum of range and Doppler costs
            costMatrix(i, j) = rangeWeight * rangeCost + dopplerWeight * dopplerCost;
        end
    end

    % Perform optimal assignment using the Munkres algorithm
    try
        [assignments, ~, ~] = assignmunkres(costMatrix, costOfNonAssignment);
    catch ME
        error('Error calling assignmunkres: %s', ME.message);
    end

    % Number of confirmed false tracks is the total confirmed tracks minus valid assignments
    numConfirmedFalseTracks = nConfirmed - size(assignments, 1);
end
