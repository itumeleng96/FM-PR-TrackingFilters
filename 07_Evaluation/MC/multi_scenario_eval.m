% =============================================================
% multi_scenario_eval.m
%
% Purpose
% -------
% Rigorous per-scenario evaluation of the four covariance-scaling
% filters (CSKF, CSUKF, CSRGNF, CSPF) against ground truth for the
% four canonical flight profiles used in the CIE paper:
%   levelFlight, landing, takeoff, orbit360.
%
% Metrics reported per (filter, scenario):
%   - RMSE_range, RMSE_doppler   (km, Hz) vs. ground truth
%   - MAE_range,  MAE_doppler    (km, Hz)
%   - Smoothness (std of scan-to-scan estimate first-difference)
%   - Outlier trigger counts (range, doppler)
%   - Mean/min NEFF (CSPF only)
%
% Currently applies Scheme A defaults (see Scheme A block below).
% =============================================================

clc; clear; close all;

addpath('01_FERS/', ...
        '04_Detection/', ...
        '05_Clustering/', ...
        '03_DPI_Cancellation/', ...
        '07_Evaluation/GroundTruth/', ...
        '06_Tracking/Filters/KF/', ...
        '06_Tracking/Filters/UKF/', ...
        '06_Tracking/Filters/RGNF/', ...
        '06_Tracking/Filters/PF/');

% -------- Scenario -> FERS XML mapping --------
scenarios = { ...
    'levelFlight', '01_FERS/BackupScenarios/scenario_1_singleFile.fersxml',   60; ...
    'landing',     '01_FERS/flightScenarios/scenario_2_landingManeuver.fersxml', 120; ...
    'takeoff',     '01_FERS/flightScenarios/scenario_3_takeoffManeuver.fersxml',  60; ...
    'orbit360',    '01_FERS/flightScenarios/scenario_4_360.fersxml',            120; ...
};

fers_bin = '/home/itumeleng/Documents/Academia/MscEng/FERS/build/src/fers';
sys_libs = '/usr/lib/gcc/x86_64-pc-linux-gnu/15:/usr/lib64';

fs          = 200000;
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

% -------- Scheme A parameters (per paper Sec IV.C tuning discussion) --------
%   - alpha  : adaptation rate on trigger  (standard adaptive-filter default)
%   - beta   : fading-memory rate off trigger (~ 1/(1-beta) scan time const)
%   - N_pf   : particle count for CSPF (compute vs. sample-noise trade-off)
% Filter tunings (std_meas, std_acc, ...) mirror runComputationalLoadRawFilter.m
% (which is what the compute-load section uses).  Alpha/beta values are the
% current hard-coded defaults in each filter's update().
dt = 1;
kf_std_meas   = [4.9038, 0.9985];
kf_std_acc    = [0.0048354, 0.0991];
ukf_std_meas  = [0.9707, 0.79739];
ukf_std_acc   = [0.0076533, 0.09938];
ukf_alpha     = 1e-4; ukf_kappa = 0; ukf_beta = 2;
pf_std_meas   = [10, 2];
pf_std_acc    = [1.429, 1.9452];
N_pf          = 5000;
rgnf_std_meas = [2.046, 0.98];
rgnf_std_acc  = [0.057027, 0.047789];
rgnf_max_iter = 100; rgnf_lambda = 1;

% -------- Container for results --------
n_scen  = size(scenarios, 1);
filters = {'CSKF', 'CSUKF', 'CSRGNF', 'CSPF'};
n_filt  = numel(filters);

metrics = struct();
for s = 1:n_scen
    for f = 1:n_filt
        metrics(s, f).scenario = scenarios{s, 1};
        metrics(s, f).filter   = filters{f};
    end
end

% =============================================================
% Main evaluation loop
% =============================================================
for s = 1:n_scen
    scen_name  = scenarios{s, 1};
    scen_xml   = scenarios{s, 2};
    scen_sec   = scenarios{s, 3};

    direct_h5 = sprintf('direct_%s.h5', scen_name);
    echo_h5   = sprintf('echo_%s.h5',   scen_name);

    fprintf('\n================ %s ================\n', scen_name);

    % ---- Ensure FERS output exists ----
    if ~exist(direct_h5, 'file') || ~exist(echo_h5, 'file')
        fprintf('  Running FERS for %s...\n', scen_name);
        fers_cmd = ['env LD_LIBRARY_PATH=' sys_libs ':$LD_LIBRARY_PATH ' ...
                    fers_bin ' ' scen_xml];
        system(fers_cmd);
        % FERS writes to hardcoded names; rename to scenario-scoped files.
        if exist('direct.h5', 'file'); movefile('direct.h5', direct_h5); end
        if exist('echo.h5', 'file');   movefile('echo.h5',   echo_h5);   end
    end

    [Ino,  Qno,  scale_no ] = loadfersHDF5(direct_h5);
    [Imov, Qmov, scale_mov] = loadfersHDF5(echo_h5);

    I_Qmov = (Imov + 1i*Qmov) .* scale_mov;
    I_Qno  = (Ino  + 1i*Qno ) .* scale_no;

    % ---- Ground truth (1 s resolution) ----
    [t_gt, true_range_km, true_doppler_hz] = computeGroundTruth(scen_name);
    fprintf('  Ground truth: %d s, %d samples.\n', scen_sec, numel(t_gt));

    % ---- Per-scan measurement pipeline ----
    N_scans      = scen_sec;
    measurements = zeros(2, N_scans);
    initial = 1; current = fs;
    for i = 1:N_scans
        s1 = I_Qmov(initial:current);
        s2 = I_Qno(initial:current);
        s1 = procECA_Optimized(s2, s1, proc);
        [y, ~] = ardNoPlot(s1, s2, fs, dopp_bins, delay, i, []);
        [targetClusters, ~, ~] = ca_cfar(y.', 1e-7, fs, dopp_bins, delay, 20);
        [clusterCentroids, ~, ~, ~] = meanShift(targetClusters, 3, 5);
        if ~isempty(clusterCentroids)
            measurements(:, i) = clusterCentroids(1:2, 1);
        elseif i > 1
            measurements(:, i) = measurements(:, i-1);
        end
        initial = current + 1;
        current = current + fs;
    end

    % Truth aligned to filter scan index (assume scan k <-> t_gt(k))
    N_align = min(N_scans, numel(t_gt));
    true_r  = true_range_km(1:N_align);
    true_d  = true_doppler_hz(1:N_align);

    initial_state = [measurements(1,1); 0; measurements(2,1); 0];

    % ---- Run each filter and score ----
    for f = 1:n_filt
        name = filters{f};
        try
            switch name
                case 'CSKF'
                    flt = CSKF(dt, kf_std_acc, kf_std_meas(1), kf_std_meas(2), initial_state);
                case 'CSUKF'
                    flt = CSUKF(dt, ukf_std_acc, ukf_std_meas(1), ukf_std_meas(2), initial_state', ukf_alpha, ukf_kappa, ukf_beta);
                case 'CSRGNF'
                    flt = CSRGNF(dt, rgnf_std_acc, rgnf_std_meas(1), rgnf_std_meas(2), initial_state, rgnf_max_iter, rgnf_lambda);
                case 'CSPF'
                    flt = CSPF(dt, pf_std_acc, pf_std_meas, initial_state, N_pf);
            end

            est_r    = zeros(1, N_scans);
            est_d    = zeros(1, N_scans);
            trig_r   = 0;
            trig_d   = 0;
            neff_trace = nan(1, N_scans);
            eps_r_prev = 0; eps_d_prev = 0;

            for k = 1:N_scans
                [~,  flt] = flt.predict();
                [xk, flt] = flt.update(measurements(:, k));
                est_r(k) = xk(1);
                est_d(k) = xk(3);

                % Trigger detection (recompute externally from the object's residual buffers).
                if isprop(flt, 'residuals_x')          % CSKF style
                    if numel(flt.residuals_x) >= 1 && numel(flt.residuals_x) < eps_r_prev
                        trig_r = trig_r + 1;
                    end
                    eps_r_prev = numel(flt.residuals_x);
                    if numel(flt.residuals_y) >= 1 && numel(flt.residuals_y) < eps_d_prev
                        trig_d = trig_d + 1;
                    end
                    eps_d_prev = numel(flt.residuals_y);
                elseif isprop(flt, 'epsRange')         % CSUKF/CSRGNF/CSPF style
                    if numel(flt.epsRange) >= 1 && numel(flt.epsRange) < eps_r_prev
                        trig_r = trig_r + 1;
                    end
                    eps_r_prev = numel(flt.epsRange);
                    if numel(flt.epsDoppler) >= 1 && numel(flt.epsDoppler) < eps_d_prev
                        trig_d = trig_d + 1;
                    end
                    eps_d_prev = numel(flt.epsDoppler);
                end

                if strcmp(name, 'CSPF')
                    neff_trace(k) = 1 / sum(flt.weights .^ 2);
                end
            end

            % ---- Metrics vs. ground truth ----
            N = N_align;
            err_r = est_r(1:N)' - true_r;
            err_d = est_d(1:N)' - true_d;

            rmse_r = sqrt(mean(err_r.^2));
            rmse_d = sqrt(mean(err_d.^2));
            mae_r  = mean(abs(err_r));
            mae_d  = mean(abs(err_d));

            % Smoothness: std of scan-to-scan first-difference of estimate.
            % Small = smooth track; large = jittery.
            smooth_r = std(diff(est_r(1:N)));
            smooth_d = std(diff(est_d(1:N)));

            metrics(s, f).rmse_r   = rmse_r;
            metrics(s, f).rmse_d   = rmse_d;
            metrics(s, f).mae_r    = mae_r;
            metrics(s, f).mae_d    = mae_d;
            metrics(s, f).smooth_r = smooth_r;
            metrics(s, f).smooth_d = smooth_d;
            metrics(s, f).trig_r   = trig_r;
            metrics(s, f).trig_d   = trig_d;
            metrics(s, f).neff_mean = mean(neff_trace, 'omitnan');
            metrics(s, f).neff_min  = min(neff_trace, [], 'omitnan');
            metrics(s, f).status = 'OK';
        catch ME
            metrics(s, f).status = ['FAIL: ' ME.message];
            fprintf('  %-7s FAIL: %s\n', name, ME.message);
        end
    end
end

% =============================================================
% Report
% =============================================================
fprintf('\n\n');
fprintf('=============================================================\n');
fprintf('  Multi-scenario evaluation  (Scheme A defaults)\n');
fprintf('=============================================================\n\n');

fprintf('%-14s %-8s | %8s %8s | %8s %8s | %8s %8s | %5s %5s\n', ...
        'Scenario','Filter','RMSE_r','RMSE_d','MAE_r','MAE_d','sm_r','sm_d','tr_r','tr_d');
fprintf('%s\n', repmat('-', 1, 100));
for s = 1:n_scen
    for f = 1:n_filt
        m = metrics(s, f);
        if strcmp(m.status, 'OK')
            fprintf('%-14s %-8s | %8.3f %8.3f | %8.3f %8.3f | %8.3f %8.3f | %5d %5d\n', ...
                    m.scenario, m.filter, m.rmse_r, m.rmse_d, m.mae_r, m.mae_d, ...
                    m.smooth_r, m.smooth_d, m.trig_r, m.trig_d);
        else
            fprintf('%-14s %-8s | %s\n', m.scenario, m.filter, m.status);
        end
    end
    if s < n_scen; fprintf('%s\n', repmat('-', 1, 100)); end
end

% ---- Per-filter aggregate (mean across scenarios) ----
fprintf('\n');
fprintf('Per-filter mean across scenarios (lower is better):\n');
fprintf('%-8s | %8s %8s | %8s %8s\n', 'Filter','RMSE_r','RMSE_d','sm_r','sm_d');
fprintf('%s\n', repmat('-', 1, 55));
for f = 1:n_filt
    rr = []; rd = []; sr = []; sd = [];
    for s = 1:n_scen
        if strcmp(metrics(s,f).status, 'OK')
            rr(end+1) = metrics(s,f).rmse_r; %#ok<*SAGROW>
            rd(end+1) = metrics(s,f).rmse_d;
            sr(end+1) = metrics(s,f).smooth_r;
            sd(end+1) = metrics(s,f).smooth_d;
        end
    end
    fprintf('%-8s | %8.3f %8.3f | %8.3f %8.3f\n', ...
            filters{f}, mean(rr), mean(rd), mean(sr), mean(sd));
end

% Save the full metrics struct for later inspection / follow-up runs.
save('multi_scenario_metrics.mat', 'metrics', 'scenarios', 'filters');
fprintf('\nMetrics saved to multi_scenario_metrics.mat\n');
