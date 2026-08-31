% =============================================================
% runComputationalLoadRawFilter.m
%
% Purpose:
%   Per-scan raw-filter runtime benchmark, bypassing the MTT wrapper.
%   Drives each filter (KF, UKF, PF, RGNF) directly with the same
%   sequence of single-target measurements so per-call cost differences
%   are visible without shared MTT overhead.
%
% Notes vs v2:
%   - Constructs each filter with the exact parameters used by
%     multiTargetTracking/track.m (cases 1, 3, 5, 7), including
%     N = 10000 particles for the PF.
%   - Measurements are the first cluster centroid from each cached scan.
%   - Times filter.predict() and filter.update(z) directly.
%   - Option B aggregation (mean over scans within an MC run, then
%     mean +/- std across 1000 per-run means).
% =============================================================

clc;
clear;
close all;

addpath('01_FERS/', ...
        '04_Detection/', ...
        '05_Clustering/', ...
        '06_Tracking/MTT/', ...
        '03_DPI_Cancellation/', ...
        '06_Tracking/Filters/KF/', ...
        '06_Tracking/Filters/KF/', ...
        '06_Tracking/Filters/PF/', ...
        '06_Tracking/Filters/PF/', ...
        '06_Tracking/Filters/UKF/', ...
        '06_Tracking/Filters/UKF/', ...
        '06_Tracking/Filters/RGNF/', ...
        '06_Tracking/Filters/RGNF/');

% -------- FERS data load (skips fers regen if h5 files exist) --------
fers_bin = '/home/itumeleng/Documents/Academia/MscEng/FERS/build/src/fers';
sys_libs = '/usr/lib/gcc/x86_64-pc-linux-gnu/15:/usr/lib64';
if ~exist('direct.h5', 'file') || ~exist('echo.h5', 'file')
    fers_cmd = ['env LD_LIBRARY_PATH=' sys_libs ':$LD_LIBRARY_PATH ' ...
                fers_bin ' FERS/BackupScenarios/scenario_1_singleFile.fersxml'];
    system(fers_cmd);
end
[Ino,  Qno,  scale_no ] = loadfersHDF5('direct.h5');
[Imov, Qmov, scale_mov] = loadfersHDF5('echo.h5');

I_Qmov = (Imov + 1i*Qmov) .* scale_mov;
I_Qno  = (Ino  + 1i*Qno ) .* scale_no;

fs              = 200000;
simulation_time = 60;
warmup_scans    = 10;
num_simulations = 1000;

dopp_bins   = 200;
delay       = 233e-6;
c           = 299792458;
range_delay = delay * c;

proc = struct('cancellationMaxRange_m',   range_delay, ...
              'cancellationMaxDoppler_Hz', 4, ...
              'TxToRefRxDistance_m',       12540, ...
              'nSegments',                 1, ...
              'nIterations',               20, ...
              'Fs',                        fs, ...
              'alpha',                     0, ...
              'initialAlpha',              0);

% =============================================================
% Step 1 — Pre-compute per-scan measurement (first cluster centroid).
% =============================================================
measurements = zeros(2, simulation_time);   % [range; doppler] per scan
initial = 1; current = fs;
fprintf('Pre-computing measurement pipeline for %d scans...\n', simulation_time);
for i = 1:simulation_time
    s1 = I_Qmov(initial:current);
    s2 = I_Qno(initial:current);
    s1 = procECA_Optimized(s2, s1, proc);
    [y, ~] = ardNoPlot(s1, s2, fs, dopp_bins, delay, i, []);
    [targetClusters, ~, ~] = ca_cfar(y.', 1e-7, fs, dopp_bins, delay, 20);
    [clusterCentroids, ~, ~, ~] = meanShift(targetClusters, 10, 8);
    if ~isempty(clusterCentroids)
        measurements(:, i) = clusterCentroids(1:2, 1);   % first centroid, [range; doppler]
    else
        % Fall back to previous scan's measurement if none found
        if i > 1
            measurements(:, i) = measurements(:, i-1);
        end
    end
    initial = current + 1;
    current = current + fs;
end
initial_state = [measurements(1,1); 0; measurements(2,1); 0];
fprintf('Measurement pipeline cached. Initial state: [%.2f, 0, %.2f, 0]\n\n', ...
        initial_state(1), initial_state(3));

% =============================================================
% Step 2 — Filter parameters (mirrored from track.m cases 1/3/5/7).
% =============================================================
dt = 1;
% KF (case 1)
kf_std_meas = [4.9038, 0.9985];
kf_std_acc  = [0.0048354, 0.0991];
% UKF (case 5)
ukf_std_acc  = [0.0076533, 0.09938];
ukf_std_meas = [0.9707, 0.79739];
ukf_alpha = 1e-4; ukf_kappa = 0; ukf_beta = 2;
% PF (case 3)
N_particles = 10000;
pf_std_acc  = [1.429, 1.9452];
pf_std_meas = [10, 2];
% RGNF (case 7)
rgnf_std_acc  = [0.057027, 0.047789];
rgnf_std_meas = [2.046, 0.98];
rgnf_max_iter = 100; rgnf_lambda = 1;

% =============================================================
% Step 3 — Monte Carlo timing loop, raw filter calls only.
% =============================================================
pred_time_kalman   = zeros(simulation_time, num_simulations);
update_time_kalman = zeros(simulation_time, num_simulations);
pred_time_ukf      = zeros(simulation_time, num_simulations);
update_time_ukf    = zeros(simulation_time, num_simulations);
pred_time_particle   = zeros(simulation_time, num_simulations);
update_time_particle = zeros(simulation_time, num_simulations);
pred_time_rgnf   = zeros(simulation_time, num_simulations);
update_time_rgnf = zeros(simulation_time, num_simulations);

fprintf('Running %d Monte Carlo simulations on raw filters...\n', num_simulations);
progress_step = max(1, floor(num_simulations / 20));

for sim = 1:num_simulations
    kf_obj   = kalmanFilter(dt, kf_std_acc, kf_std_meas(1), kf_std_meas(2), initial_state);
    ukf_obj  = unscentedKalmanFilter(dt, ukf_std_acc, ukf_std_meas(1), ukf_std_meas(2), initial_state', ukf_alpha, ukf_kappa, ukf_beta);
    pf_obj   = particleFilter(dt, pf_std_acc, pf_std_meas, initial_state, N_particles);
    rgnf_obj = RGNF(dt, rgnf_std_acc, rgnf_std_meas(1), rgnf_std_meas(2), initial_state, rgnf_max_iter, rgnf_lambda);

    for i = 1:simulation_time
        z = measurements(:, i);

        % Prediction stage
        tic; [~, kf_obj]   = kf_obj.predict();   pred_time_kalman(i, sim)   = toc;
        tic; [~, ukf_obj]  = ukf_obj.predict();  pred_time_ukf(i, sim)      = toc;
        tic; [~, pf_obj]   = pf_obj.predict();   pred_time_particle(i, sim) = toc;
        tic; [~, rgnf_obj] = rgnf_obj.predict(); pred_time_rgnf(i, sim)     = toc;

        % Update stage
        tic; [~, kf_obj]   = kf_obj.update(z);   update_time_kalman(i, sim)   = toc;
        tic; [~, ukf_obj]  = ukf_obj.update(z);  update_time_ukf(i, sim)      = toc;
        tic; [~, pf_obj]   = pf_obj.update(z);   update_time_particle(i, sim) = toc;
        tic; [~, rgnf_obj] = rgnf_obj.update(z); update_time_rgnf(i, sim)     = toc;
    end

    if mod(sim, progress_step) == 0
        fprintf('  %d / %d MC runs done\n', sim, num_simulations);
    end
end

% =============================================================
% Step 4 — Option B aggregation + print + save.
% =============================================================
[mu_pKF,  sd_pKF]  = agg(pred_time_kalman,     warmup_scans);
[mu_uKF,  sd_uKF]  = agg(update_time_kalman,   warmup_scans);
[mu_pUKF, sd_pUKF] = agg(pred_time_ukf,        warmup_scans);
[mu_uUKF, sd_uUKF] = agg(update_time_ukf,      warmup_scans);
[mu_pPF,  sd_pPF]  = agg(pred_time_particle,   warmup_scans);
[mu_uPF,  sd_uPF]  = agg(update_time_particle, warmup_scans);
[mu_pRG,  sd_pRG]  = agg(pred_time_rgnf,       warmup_scans);
[mu_uRG,  sd_uRG]  = agg(update_time_rgnf,     warmup_scans);

r_uUKF = mu_uUKF / mu_uKF;   r_pUKF = mu_pUKF / mu_pKF;
r_uPF  = mu_uPF  / mu_uKF;   r_pPF  = mu_pPF  / mu_pKF;
r_uRG  = mu_uRG  / mu_uKF;   r_pRG  = mu_pRG  / mu_pKF;

fprintf('\n=========================================================\n');
fprintf(' RAW FILTER per-scan runtime (\\mu s), mean +/- std over %d MC runs\n', num_simulations);
fprintf(' Warm-up: first %d scans of each run discarded\n', warmup_scans);
fprintf('=========================================================\n');
fprintf('%-6s | %-18s | %-8s | %-18s | %-8s\n', ...
        'Filter', 'Update (mu +/- s)', 'vs KF', 'Predict (mu +/- s)', 'vs KF');
fprintf('%-6s | %6.2f +/- %-6.2f  | %6.2f  | %6.2f +/- %-6.2f  | %6.2f\n', ...
        'KF',   mu_uKF,  sd_uKF,  1.00,   mu_pKF,  sd_pKF,  1.00);
fprintf('%-6s | %6.2f +/- %-6.2f  | %6.2f  | %6.2f +/- %-6.2f  | %6.2f\n', ...
        'UKF',  mu_uUKF, sd_uUKF, r_uUKF, mu_pUKF, sd_pUKF, r_pUKF);
fprintf('%-6s | %6.2f +/- %-6.2f  | %6.2f  | %6.2f +/- %-6.2f  | %6.2f\n', ...
        'PF',   mu_uPF,  sd_uPF,  r_uPF,  mu_pPF,  sd_pPF,  r_pPF);
fprintf('%-6s | %6.2f +/- %-6.2f  | %6.2f  | %6.2f +/- %-6.2f  | %6.2f\n', ...
        'RGNF', mu_uRG,  sd_uRG,  r_uRG,  mu_pRG,  sd_pRG,  r_pRG);
fprintf('=========================================================\n');

save('evaluationScripts/ComputationalLoad/computational_load_raw_filter_results.mat', ...
     'pred_time_kalman', 'update_time_kalman', ...
     'pred_time_ukf', 'update_time_ukf', ...
     'pred_time_particle', 'update_time_particle', ...
     'pred_time_rgnf', 'update_time_rgnf', ...
     'num_simulations', 'simulation_time', 'warmup_scans');

fprintf('\nRaw timing matrices saved to evaluationScripts/ComputationalLoad/computational_load_raw_filter_results.mat\n');

% =============================================================
% Local functions
% =============================================================
function [mu_us, sd_us] = agg(t_matrix, warmup)
    keep = t_matrix(warmup+1:end, :);
    per_run_mean = mean(keep, 1);
    mu_us = mean(per_run_mean) * 1e6;
    sd_us = std(per_run_mean)  * 1e6;
end
