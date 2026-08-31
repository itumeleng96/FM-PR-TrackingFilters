% =============================================================
% scheme_b_grid_sweep.m
%
% Purpose
% -------
% Scheme B tuning: grid sweep over the covariance-scaling parameters
% (cs_alpha, cs_beta) AND a process-noise scale factor (Q_scale)
% for each of the four filters (CSKF, CSUKF, CSRGNF, CSPF), evaluated
% across all four canonical scenarios (levelFlight, landing, takeoff,
% orbit360).
%
% Optimisation criterion
% ----------------------
% Per (filter, combo): compute per-scenario RMSE_r (km) and RMSE_d (Hz).
% Normalise each metric column-wise across scenarios (so no single
% scenario dominates), then take mean over scenarios and mean of the
% two normalised metrics. Pick the combo with the smallest score per
% filter -> per-filter "best" tuning.
%
% Runtime
% -------
% Measurements are cached to <scenario>_measurements.mat after the
% first pass (so FERS + CFAR + mean-shift run once per scenario).
% Subsequent tracker runs are cheap.
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

% -------- Baseline filter tunings (mirror runComputationalLoadRawFilter.m) --
dt = 1;
kf_std_meas   = [4.9038, 0.9985];
kf_std_acc0   = [0.0048354, 0.0991];
ukf_std_meas  = [0.9707, 0.79739];
ukf_std_acc0  = [0.0076533, 0.09938];
ukf_alpha     = 1e-4; ukf_kappa = 0; ukf_beta = 2;
pf_std_meas   = [10, 2];
pf_std_acc0   = [1.429, 1.9452];
N_pf          = 3000;                         % reduced from 5000 to speed sweep
rgnf_std_meas = [2.046, 0.98];
rgnf_std_acc0 = [0.057027, 0.047789];
rgnf_max_iter = 100; rgnf_lambda = 1;

% -------- Scheme B grid --------
alpha_grid   = [0.3, 0.6, 0.9];        % adaptation rate on trigger
beta_grid    = [0.90, 0.98];           % fading-memory rate off trigger
qscale_grid  = [0.5, 1.0, 2.0, 5.0];   % scale factor on baseline std_acc
combos       = combvec(alpha_grid, beta_grid, qscale_grid).';
n_combo      = size(combos, 1);

fprintf('Scheme B: %d combos per filter (alpha=%d x beta=%d x Qscale=%d).\n', ...
        n_combo, numel(alpha_grid), numel(beta_grid), numel(qscale_grid));

% =============================================================
% Stage 1: precompute measurements per scenario (cached)
% =============================================================
n_scen = size(scenarios, 1);
scen_data = cell(n_scen, 1);

for s = 1:n_scen
    scen_name = scenarios{s, 1};
    scen_xml  = scenarios{s, 2};
    scen_sec  = scenarios{s, 3};
    cache     = sprintf('meas_%s.mat', scen_name);

    if exist(cache, 'file')
        fprintf('  Loading cached measurements for %s...\n', scen_name);
        S = load(cache);
        scen_data{s} = S.data;
        continue;
    end

    fprintf('  Extracting measurements for %s...\n', scen_name);
    direct_h5 = sprintf('direct_%s.h5', scen_name);
    echo_h5   = sprintf('echo_%s.h5',   scen_name);

    if ~exist(direct_h5, 'file') || ~exist(echo_h5, 'file')
        fprintf('    Running FERS...\n');
        fers_cmd = ['env LD_LIBRARY_PATH=' sys_libs ':$LD_LIBRARY_PATH ' ...
                    fers_bin ' ' scen_xml];
        system(fers_cmd);
        if exist('direct.h5', 'file'); movefile('direct.h5', direct_h5); end
        if exist('echo.h5', 'file');   movefile('echo.h5',   echo_h5);   end
    end

    [Ino,  Qno,  scale_no ] = loadfersHDF5(direct_h5);
    [Imov, Qmov, scale_mov] = loadfersHDF5(echo_h5);

    I_Qmov = (Imov + 1i*Qmov) .* scale_mov;
    I_Qno  = (Ino  + 1i*Qno ) .* scale_no;

    [t_gt, true_range_km, true_doppler_hz] = computeGroundTruth(scen_name);

    N_scans = scen_sec;
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

    N_align = min(N_scans, numel(t_gt));
    data = struct( ...
        'name',    scen_name, ...
        'N',       N_scans, ...
        'N_align', N_align, ...
        'meas',    measurements, ...
        'true_r',  true_range_km(1:N_align), ...
        'true_d',  true_doppler_hz(1:N_align));
    scen_data{s} = data;

    save(cache, 'data');
    fprintf('    Cached to %s.\n', cache);
end

% =============================================================
% Stage 2: per-filter grid sweep
% =============================================================
filters = {'CSKF', 'CSUKF', 'CSRGNF', 'CSPF'};
n_filt  = numel(filters);

% results{f} : n_combo x n_scen x 2 (rmse_r, rmse_d)
results = cell(n_filt, 1);
for f = 1:n_filt
    results{f} = nan(n_combo, n_scen, 2);
end

for f = 1:n_filt
    name = filters{f};
    fprintf('\n---- Sweeping %s (%d combos x %d scenarios) ----\n', ...
            name, n_combo, n_scen);
    for c = 1:n_combo
        a  = combos(c, 1);
        b  = combos(c, 2);
        qs = combos(c, 3);
        for s = 1:n_scen
            d = scen_data{s};
            initial_state = [d.meas(1,1); 0; d.meas(2,1); 0];
            try
                switch name
                    case 'CSKF'
                        flt = CSKF(dt, kf_std_acc0 * qs, kf_std_meas(1), kf_std_meas(2), initial_state);
                    case 'CSUKF'
                        flt = CSUKF(dt, ukf_std_acc0 * qs, ukf_std_meas(1), ukf_std_meas(2), initial_state', ukf_alpha, ukf_kappa, ukf_beta);
                    case 'CSRGNF'
                        flt = CSRGNF(dt, rgnf_std_acc0 * qs, rgnf_std_meas(1), rgnf_std_meas(2), initial_state, rgnf_max_iter, rgnf_lambda);
                    case 'CSPF'
                        flt = CSPF(dt, pf_std_acc0 * qs, pf_std_meas, initial_state, N_pf);
                end
                flt.cs_alpha = a;
                flt.cs_beta  = b;

                est_r = zeros(1, d.N);
                est_d = zeros(1, d.N);
                for k = 1:d.N
                    [~,  flt] = flt.predict();
                    [xk, flt] = flt.update(d.meas(:, k));
                    est_r(k) = xk(1);
                    est_d(k) = xk(3);
                end

                N = d.N_align;
                err_r = est_r(1:N)' - d.true_r;
                err_d = est_d(1:N)' - d.true_d;
                results{f}(c, s, 1) = sqrt(mean(err_r.^2));
                results{f}(c, s, 2) = sqrt(mean(err_d.^2));
            catch ME
                fprintf('   %s combo %d scen %s FAIL: %s\n', name, c, d.name, ME.message);
            end
        end
        if mod(c, 6) == 0
            fprintf('    ...combo %d/%d\n', c, n_combo);
        end
    end
end

% =============================================================
% Stage 3: pick best combo per filter (normalised score)
% =============================================================
best = struct();
for f = 1:n_filt
    R = results{f};                    % n_combo x n_scen x 2

    % Normalise each (scenario, metric) column to [0, 1] using min/max
    % across combos (so bad-in-worst-case combos are penalised, and no
    % single scenario dominates the ranking).
    Rn = R;
    for s = 1:n_scen
        for m = 1:2
            col = R(:, s, m);
            lo  = min(col, [], 'omitnan');
            hi  = max(col, [], 'omitnan');
            if hi > lo
                Rn(:, s, m) = (col - lo) ./ (hi - lo);
            else
                Rn(:, s, m) = 0;
            end
        end
    end

    % Aggregate: mean over scenarios AND both metrics.
    score = mean(mean(Rn, 3, 'omitnan'), 2, 'omitnan');
    [best_score, best_c] = min(score);

    best(f).filter    = filters{f};
    best(f).combo_idx = best_c;
    best(f).cs_alpha  = combos(best_c, 1);
    best(f).cs_beta   = combos(best_c, 2);
    best(f).q_scale   = combos(best_c, 3);
    best(f).score     = best_score;
    best(f).rmse_r    = squeeze(R(best_c, :, 1));
    best(f).rmse_d    = squeeze(R(best_c, :, 2));
end

% =============================================================
% Report
% =============================================================
fprintf('\n\n=============================================================\n');
fprintf('  Scheme B grid-sweep results\n');
fprintf('=============================================================\n\n');

fprintf('Best tuning per filter (normalised-RMSE score, lower better):\n');
fprintf('%-8s | %8s %8s %8s | %8s\n', 'Filter','cs_alpha','cs_beta','Q_scale','score');
fprintf('%s\n', repmat('-', 1, 55));
for f = 1:n_filt
    fprintf('%-8s | %8.2f %8.2f %8.2f | %8.4f\n', ...
            best(f).filter, best(f).cs_alpha, best(f).cs_beta, ...
            best(f).q_scale, best(f).score);
end

fprintf('\nPer-scenario RMSE at best tuning:\n');
fprintf('%-8s | ', 'Filter');
for s = 1:n_scen
    fprintf('%-22s', scenarios{s,1});
end
fprintf('\n%s\n', repmat('-', 1, 8 + 3 + 22*n_scen));
for f = 1:n_filt
    fprintf('%-8s | ', best(f).filter);
    for s = 1:n_scen
        fprintf('%7.2f km / %6.1f Hz  ', best(f).rmse_r(s), best(f).rmse_d(s));
    end
    fprintf('\n');
end

save('scheme_b_results.mat', 'results', 'best', 'combos', 'scenarios', 'filters');
fprintf('\nSaved to scheme_b_results.mat\n');
