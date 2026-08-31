% =============================================================
% runComputationalLoadScaling.m
%
% Purpose:
%   Definitive MTT-inclusive per-scan cost measurement across
%   N = {1, 3, 5, 10} target scenarios, for each of KF / UKF / PF / RGNF.
%   Produces the scaling data for Paragraph 3 of Section IV.C.
%
% Correct MTT flow (per runTrackingFilter.m:140-153):
%     mtt = mtt.createNewTracks(centroids, i);
%     mtt = mtt.maintainTracks();
%     mtt = mtt.predictionStage();
%     mtt = mtt.updateStage(centroids, i);
%
% Methodology:
%   - 1000 Monte Carlo runs per scenario.
%   - First 10 scans of every run discarded as warm-up (MTT confirmation
%     ramp + MATLAB JIT compilation transient).
%   - Only predictionStage() and updateStage() are timed. createNewTracks
%     and maintainTracks are shared across all filters and would only
%     add constant per-scan noise to every timing.
%   - Option B aggregation: within a run, compute per-scan mean; then
%     report mean +/- std across the 1000 per-run means.
%   - Track n_confirmed_per_scan alongside timing so we can report the
%     actual mean N_t observed under each scenario.
%
% Scenarios:
%   Nominal N | Fersxml                                        | h5 files
%   1         | scenario_1_singleFile.fersxml   (existing)     | direct.h5     / echo.h5
%   3         | scenario_N3_targets.fersxml     (redesigned)   | direct_N3.h5  / echo_N3.h5
%   5         | scenario_N5_targets.fersxml     (redesigned)   | direct_N5.h5  / echo_N5.h5
%   10        | scenario_N10_targets.fersxml    (redesigned)   | direct_N10.h5 / echo_N10.h5
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

fs              = 200000;
simulation_time = 60;
warmup_scans    = 10;
num_simulations = 1000;

dopp_bins  = 200;
delay      = 233e-6;
c          = 299792458;
range_delay = delay * c;

proc = struct('cancellationMaxRange_m',   range_delay, ...
              'cancellationMaxDoppler_Hz', 4, ...
              'TxToRefRxDistance_m',       12540, ...
              'nSegments',                 1, ...
              'nIterations',               20, ...
              'Fs',                        fs, ...
              'alpha',                     0, ...
              'initialAlpha',              0);

confirmationThreshold = 4;
deletionThreshold     = 4;
gatingThreshold       = 20;

fers_bin = '/home/itumeleng/Documents/Academia/MscEng/FERS/build/src/fers';
sys_libs = '/usr/lib/gcc/x86_64-pc-linux-gnu/15:/usr/lib64';

% -------- Scenarios --------
scenarios = struct( ...
    'nominal',   {1, 3, 5, 10}, ...
    'label',     {'N=1', 'N=3', 'N=5', 'N=10'}, ...
    'fersxml',   {'01_FERS/BackupScenarios/scenario_1_singleFile.fersxml', ...
                  '01_FERS/BackupScenarios/scenario_N3_targets.fersxml', ...
                  '01_FERS/BackupScenarios/scenario_N5_targets.fersxml', ...
                  '01_FERS/BackupScenarios/scenario_N10_targets.fersxml'}, ...
    'direct_h5', {'direct.h5',    'direct_N3.h5', 'direct_N5.h5', 'direct_N10.h5'}, ...
    'echo_h5',   {'echo.h5',      'echo_N3.h5',   'echo_N5.h5',   'echo_N10.h5'});

% -------- Result container --------
results = struct();

for s = 1:numel(scenarios)
    sc = scenarios(s);
    fprintf('\n============================================================\n');
    fprintf(' Scenario %s (nominal N_target = %d)\n', sc.label, sc.nominal);
    fprintf('============================================================\n');

    % Ensure h5 files exist
    if ~exist(sc.direct_h5, 'file') || ~exist(sc.echo_h5, 'file')
        fprintf('Running FERS on %s ...\n', sc.fersxml);
        system(['env LD_LIBRARY_PATH=' sys_libs ':$LD_LIBRARY_PATH ' fers_bin ' ' sc.fersxml]);
    end

    [Ino,  Qno,  scale_no ] = loadfersHDF5(sc.direct_h5);
    [Imov, Qmov, scale_mov] = loadfersHDF5(sc.echo_h5);
    I_Qmov = (Imov + 1i*Qmov) .* scale_mov;
    I_Qno  = (Ino  + 1i*Qno ) .* scale_no;

    % ------- Step 1: Pre-compute per-scan cluster centroids -------
    clusterCentroidsCache = cell(simulation_time, 1);
    initial = 1; current = fs;
    fprintf('Pre-computing measurement pipeline for %d scans...\n', simulation_time);
    for i = 1:simulation_time
        s1 = I_Qmov(initial:current);
        s2 = I_Qno(initial:current);
        s1 = procECA_Optimized(s2, s1, proc);
        [y, ~] = ardNoPlot(s1, s2, fs, dopp_bins, delay, i, []);
        [targetClusters, ~, ~] = ca_cfar(y.', 1e-5, fs, dopp_bins, delay, 20);
        [clusterCentroids, ~, ~, ~] = meanShift(targetClusters, 10, 8);
        clusterCentroidsCache{i} = clusterCentroids;
        initial = current + 1;
        current = current + fs;
    end
    fprintf('Measurement pipeline cached.\n');

    % ------- Step 2: MC timing loop with correct MTT flow -------
    pred_time_kalman   = zeros(simulation_time, num_simulations);
    update_time_kalman = zeros(simulation_time, num_simulations);
    pred_time_ukf      = zeros(simulation_time, num_simulations);
    update_time_ukf    = zeros(simulation_time, num_simulations);
    pred_time_pf       = zeros(simulation_time, num_simulations);
    update_time_pf     = zeros(simulation_time, num_simulations);
    pred_time_rgnf     = zeros(simulation_time, num_simulations);
    update_time_rgnf   = zeros(simulation_time, num_simulations);

    % Track confirmed count per scan per filter (mean across MC)
    n_conf_kalman = zeros(simulation_time, num_simulations);
    n_conf_ukf    = zeros(simulation_time, num_simulations);
    n_conf_pf     = zeros(simulation_time, num_simulations);
    n_conf_rgnf   = zeros(simulation_time, num_simulations);

    fprintf('Running %d Monte Carlo simulations...\n', num_simulations);
    progress_step = max(1, floor(num_simulations / 20));

    for sim = 1:num_simulations
        mttKF   = multiTargetTracker(confirmationThreshold, deletionThreshold, gatingThreshold, 1);
        mttPF   = multiTargetTracker(confirmationThreshold, deletionThreshold, gatingThreshold, 3);
        mttUKF  = multiTargetTracker(confirmationThreshold, deletionThreshold, gatingThreshold, 5);
        mttRGNF = multiTargetTracker(confirmationThreshold, deletionThreshold, gatingThreshold, 7);

        for i = 1:simulation_time
            centroids = clusterCentroidsCache{i};

            % Setup calls (not timed) - identical for every filter
            mttKF   = mttKF.createNewTracks(centroids, i);
            mttKF   = mttKF.maintainTracks();
            mttPF   = mttPF.createNewTracks(centroids, i);
            mttPF   = mttPF.maintainTracks();
            mttUKF  = mttUKF.createNewTracks(centroids, i);
            mttUKF  = mttUKF.maintainTracks();
            mttRGNF = mttRGNF.createNewTracks(centroids, i);
            mttRGNF = mttRGNF.maintainTracks();

            % Prediction stage (timed)
            tic; mttPF   = mttPF.predictionStage();   pred_time_pf(i, sim)     = toc;
            tic; mttKF   = mttKF.predictionStage();   pred_time_kalman(i, sim) = toc;
            tic; mttUKF  = mttUKF.predictionStage();  pred_time_ukf(i, sim)    = toc;
            tic; mttRGNF = mttRGNF.predictionStage(); pred_time_rgnf(i, sim)   = toc;

            % Update stage (timed)
            tic; mttPF   = mttPF.updateStage(centroids, i);   update_time_pf(i, sim)     = toc;
            tic; mttKF   = mttKF.updateStage(centroids, i);   update_time_kalman(i, sim) = toc;
            tic; mttUKF  = mttUKF.updateStage(centroids, i);  update_time_ukf(i, sim)    = toc;
            tic; mttRGNF = mttRGNF.updateStage(centroids, i); update_time_rgnf(i, sim)   = toc;

            % Confirmed-track counts (post-update, per filter, per scan)
            n_conf_kalman(i, sim) = countConfirmed(mttKF);
            n_conf_ukf(i, sim)    = countConfirmed(mttUKF);
            n_conf_pf(i, sim)     = countConfirmed(mttPF);
            n_conf_rgnf(i, sim)   = countConfirmed(mttRGNF);
        end

        if mod(sim, progress_step) == 0
            fprintf('  %d / %d MC runs done\n', sim, num_simulations);
        end
    end

    % ------- Step 3: Option B aggregation -------
    r = struct();
    [r.mu_pKF,  r.sd_pKF]  = agg(pred_time_kalman,   warmup_scans);
    [r.mu_uKF,  r.sd_uKF]  = agg(update_time_kalman, warmup_scans);
    [r.mu_pUKF, r.sd_pUKF] = agg(pred_time_ukf,      warmup_scans);
    [r.mu_uUKF, r.sd_uUKF] = agg(update_time_ukf,    warmup_scans);
    [r.mu_pPF,  r.sd_pPF]  = agg(pred_time_pf,       warmup_scans);
    [r.mu_uPF,  r.sd_uPF]  = agg(update_time_pf,     warmup_scans);
    [r.mu_pRG,  r.sd_pRG]  = agg(pred_time_rgnf,     warmup_scans);
    [r.mu_uRG,  r.sd_uRG]  = agg(update_time_rgnf,   warmup_scans);

    % Actual mean confirmed count per filter (aggregate across warm-up-excluded scans and runs)
    r.n_conf_kalman = mean(mean(n_conf_kalman(warmup_scans+1:end, :), 1));
    r.n_conf_ukf    = mean(mean(n_conf_ukf   (warmup_scans+1:end, :), 1));
    r.n_conf_pf     = mean(mean(n_conf_pf    (warmup_scans+1:end, :), 1));
    r.n_conf_rgnf   = mean(mean(n_conf_rgnf  (warmup_scans+1:end, :), 1));

    r.label   = sc.label;
    r.nominal = sc.nominal;
    results(s).r = r;

    fprintf('\n----- %s per-scan cost (mean +/- std across %d MC runs) -----\n', ...
            sc.label, num_simulations);
    fprintf('%-6s | Update (us)      | Predict (us)     | mean N_conf\n', 'Filter');
    fprintf('%-6s | %8.2f +/- %-6.2f | %8.2f +/- %-6.2f | %5.2f\n', ...
            'KF',   r.mu_uKF,  r.sd_uKF,  r.mu_pKF,  r.sd_pKF,  r.n_conf_kalman);
    fprintf('%-6s | %8.2f +/- %-6.2f | %8.2f +/- %-6.2f | %5.2f\n', ...
            'UKF',  r.mu_uUKF, r.sd_uUKF, r.mu_pUKF, r.sd_pUKF, r.n_conf_ukf);
    fprintf('%-6s | %8.2f +/- %-6.2f | %8.2f +/- %-6.2f | %5.2f\n', ...
            'PF',   r.mu_uPF,  r.sd_uPF,  r.mu_pPF,  r.sd_pPF,  r.n_conf_pf);
    fprintf('%-6s | %8.2f +/- %-6.2f | %8.2f +/- %-6.2f | %5.2f\n', ...
            'RGNF', r.mu_uRG,  r.sd_uRG,  r.mu_pRG,  r.sd_pRG,  r.n_conf_rgnf);
end

% =============================================================
% Step 4: Scaling summary Table X.
% Total per-scan cost = update + predict, in microseconds.
% =============================================================
fprintf('\n============================================================\n');
fprintf(' TABLE X: MTT-inclusive per-scan cost (Update + Predict, us)\n');
fprintf(' vs number of nominal targets N. Mean over %d MC runs.\n', num_simulations);
fprintf('============================================================\n');
fprintf('%-6s | %-14s | %-14s | %-14s | %-14s\n', ...
        'Filter', 'N=1', 'N=3', 'N=5', 'N=10');
for f = {'KF', 'UKF', 'PF', 'RGNF'}
    fname = f{1};
    row = sprintf('%-6s |', fname);
    for s = 1:numel(scenarios)
        r = results(s).r;
        switch fname
            case 'KF';   tot = r.mu_uKF  + r.mu_pKF;
            case 'UKF';  tot = r.mu_uUKF + r.mu_pUKF;
            case 'PF';   tot = r.mu_uPF  + r.mu_pPF;
            case 'RGNF'; tot = r.mu_uRG  + r.mu_pRG;
        end
        row = [row, sprintf(' %13.2f  |', tot)]; %#ok<AGROW>
    end
    fprintf('%s\n', row);
end

fprintf('\nMean confirmed N_t observed per scenario:\n');
fprintf('%-6s | %-6s | %-6s | %-6s | %-6s\n', 'Filter', 'N=1', 'N=3', 'N=5', 'N=10');
for f = {'KF', 'UKF', 'PF', 'RGNF'}
    fname = f{1};
    row = sprintf('%-6s |', fname);
    for s = 1:numel(scenarios)
        r = results(s).r;
        switch fname
            case 'KF';   v = r.n_conf_kalman;
            case 'UKF';  v = r.n_conf_ukf;
            case 'PF';   v = r.n_conf_pf;
            case 'RGNF'; v = r.n_conf_rgnf;
        end
        row = [row, sprintf(' %5.2f  |', v)]; %#ok<AGROW>
    end
    fprintf('%s\n', row);
end

% Scaling ratio N=10 / N=1 (linear scaling predicts approx equal to mean N_conf ratio)
fprintf('\nScaling ratio (N=10 total cost) / (N=1 total cost):\n');
for f = {'KF', 'UKF', 'PF', 'RGNF'}
    fname = f{1};
    r1  = results(1).r;
    r10 = results(4).r;
    switch fname
        case 'KF';   tot1 = r1.mu_uKF  + r1.mu_pKF;  tot10 = r10.mu_uKF  + r10.mu_pKF;
        case 'UKF';  tot1 = r1.mu_uUKF + r1.mu_pUKF; tot10 = r10.mu_uUKF + r10.mu_pUKF;
        case 'PF';   tot1 = r1.mu_uPF  + r1.mu_pPF;  tot10 = r10.mu_uPF  + r10.mu_pPF;
        case 'RGNF'; tot1 = r1.mu_uRG  + r1.mu_pRG;  tot10 = r10.mu_uRG  + r10.mu_pRG;
    end
    fprintf('  %-4s : cost ratio = %5.2fx, N_conf ratio = %5.2fx\n', ...
            fname, tot10/tot1, ...
            get_nconf(results(4).r, fname) / max(get_nconf(results(1).r, fname), eps));
end

save('evaluationScripts/ComputationalLoad/computational_load_scaling_results.mat', 'results', ...
     'num_simulations', 'simulation_time', 'warmup_scans');
fprintf('\nSaved evaluationScripts/ComputationalLoad/computational_load_scaling_results.mat\n');

% =============================================================
% Local functions
% =============================================================
function [mu_us, sd_us] = agg(t_matrix, warmup)
    keep = t_matrix(warmup+1:end, :);
    per_run_mean = mean(keep, 1);
    mu_us = mean(per_run_mean) * 1e6;
    sd_us = std(per_run_mean)  * 1e6;
end

function n = countConfirmed(mtt)
    n = 0;
    for k = 1:numel(mtt.tracks)
        if ~mtt.tracks(k).deleted && mtt.tracks(k).confirmed == 1
            n = n + 1;
        end
    end
end

function v = get_nconf(r, fname)
    switch fname
        case 'KF';   v = r.n_conf_kalman;
        case 'UKF';  v = r.n_conf_ukf;
        case 'PF';   v = r.n_conf_pf;
        case 'RGNF'; v = r.n_conf_rgnf;
    end
end
