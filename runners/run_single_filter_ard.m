function out = run_single_filter_ard(filter_name, scenario_name, snapshot_scans)
% RUN_SINGLE_FILTER_ARD  Run one filter on one scenario, produce ARD plots.
%
%   out = run_single_filter_ard('CSPF', 'landing')
%   out = run_single_filter_ard('CSKF', 'orbit360', [10 40 80 110])
%
% Produces:
%   Fig 1: Range-time waterfall (ARD max-projected over doppler) with
%          truth, measurement, and filter estimate overlaid.
%   Fig 2: Doppler-time waterfall (ARD max-projected over range) with
%          overlays.
%   Fig 3: ARD snapshot montage at snapshot_scans (default: 4 evenly
%          spaced scans through the run). Filter estimate marked with
%          a red '+' on each frame.
%
% Also prints RMSE and log-likelihood (same metrics as run_single_filter).
%
% Requires FERS/CFAR/meanShift toolchain + cached FERS h5 files
% (direct_<scenario>.h5, echo_<scenario>.h5) from earlier runs.

    addpath('FERS/', 'cfar/', 'meanShiftCluster/', 'DPI_Suppression/', ...
            'groundTruthCalculations/', ...
            'TrackingFilter-CSKF/',   'TrackingFilter-CSUKF/', ...
            'TrackingFilter-CSRGNF/', 'TrackingFilter-CS-ParticleFilter/');

    scen_secs = struct('levelFlight', 60, 'landing', 120, 'takeoff', 60, 'orbit360', 120);
    if ~isfield(scen_secs, scenario_name)
        error('Unknown scenario "%s".', scenario_name);
    end
    N_scans = scen_secs.(scenario_name);

    direct_h5 = sprintf('direct_%s.h5', scenario_name);
    echo_h5   = sprintf('echo_%s.h5',   scenario_name);
    if ~exist(direct_h5, 'file') || ~exist(echo_h5, 'file')
        error('Missing %s or %s. Run scheme_b_grid_sweep.m first.', direct_h5, echo_h5);
    end

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

    % ---- Filter tuning (Scheme B2) ----
    dt = 1;
    baselines = struct();
    baselines.CSKF   = struct('std_meas', [4.9038, 0.9985],   'std_acc', [0.0048354, 0.0991]);
    baselines.CSUKF  = struct('std_meas', [0.9707, 0.79739],  'std_acc', [0.0076533, 0.09938], ...
                              'alpha', 1e-4, 'kappa', 0, 'beta', 2);
    baselines.CSRGNF = struct('std_meas', [2.046, 0.98],      'std_acc', [0.057027, 0.047789], ...
                              'max_iter', 100, 'lambda', 1);
    baselines.CSPF   = struct('std_meas', [10, 2],            'std_acc', [1.429, 1.9452], 'N', 3000);
    tuning = struct();
    tuning.CSKF   = struct('a', 0.90, 'b', 0.90, 'qr', 2.0,  'qd', 1.0);
    tuning.CSUKF  = struct('a', 0.60, 'b', 0.90, 'qr', 0.5,  'qd', 5.0);
    tuning.CSRGNF = struct('a', 0.90, 'b', 0.90, 'qr', 1.0,  'qd', 10.0);
    tuning.CSPF   = struct('a', 0.60, 'b', 0.90, 'qr', 1.0,  'qd', 2.0);
    if ~isfield(baselines, filter_name)
        error('Unknown filter "%s".', filter_name);
    end
    B  = baselines.(filter_name);
    tp = tuning.(filter_name);

    % ---- Ground truth ----
    [~, true_r, true_d] = computeGroundTruth(scenario_name);
    N = min(N_scans, numel(true_r));

    % ---- Load FERS complex baseband ----
    fprintf('Loading FERS data for %s...\n', scenario_name);
    [Ino,  Qno,  scale_no ] = loadfersHDF5(direct_h5);
    [Imov, Qmov, scale_mov] = loadfersHDF5(echo_h5);
    I_Qmov = (Imov + 1i*Qmov) .* scale_mov;
    I_Qno  = (Ino  + 1i*Qno ) .* scale_no;

    % ---- ARD axes ----
    N_samp   = fs;                           % 1-second window
    Ndelay   = floor(delay * fs);            % range bins - 1
    Ndop     = ceil(N_samp * dopp_bins / fs);% one-sided doppler bin count
    range_axis_km = (0:Ndelay).' * (c/fs) / 1000;  % (Ndelay+1) x 1
    dopp_axis_hz  = (-Ndop:Ndop).';                % (2*Ndop+1) x 1

    % ---- Preallocate ARD stack ----
    ard_stack = zeros(numel(range_axis_km), numel(dopp_axis_hz), N_scans, 'single');
    measurements = zeros(2, N_scans);

    initial = 1; current = fs;
    fprintf('Building ARD stack (%d scans)...\n', N_scans);
    for i = 1:N_scans
        s1 = I_Qmov(initial:current);
        s2 = I_Qno(initial:current);
        s1 = procECA_Optimized(s2, s1, proc);
        [y, ~] = ardNoPlot(s1, s2, fs, dopp_bins, delay, i, []);
        ard_stack(:, :, i) = single(y);

        [targetClusters, ~, ~] = ca_cfar(y.', 1e-7, fs, dopp_bins, delay, 20);
        [clusterCentroids, ~, ~, ~] = meanShift(targetClusters, 3, 5);
        if ~isempty(clusterCentroids)
            measurements(:, i) = clusterCentroids(1:2, 1);
        elseif i > 1
            measurements(:, i) = measurements(:, i-1);
        end
        initial = current + 1;
        current = current + fs;
        if mod(i, 20) == 0
            fprintf('  ...scan %d/%d\n', i, N_scans);
        end
    end

    % ---- Run filter with log-likelihood accumulation ----
    initial_state = [measurements(1,1); 0; measurements(2,1); 0];
    sa = [B.std_acc(1) * tp.qr, B.std_acc(2) * tp.qd];
    switch filter_name
        case 'CSKF'
            flt = CSKF(dt, sa, B.std_meas(1), B.std_meas(2), initial_state);
        case 'CSUKF'
            flt = CSUKF(dt, sa, B.std_meas(1), B.std_meas(2), initial_state', B.alpha, B.kappa, B.beta);
        case 'CSRGNF'
            flt = CSRGNF(dt, sa, B.std_meas(1), B.std_meas(2), initial_state, B.max_iter, B.lambda);
        case 'CSPF'
            flt = CSPF(dt, sa, B.std_meas, initial_state, B.N);
    end
    flt.cs_alpha = tp.a;
    flt.cs_beta  = tp.b;

    est_r = zeros(1, N_scans);
    est_d = zeros(1, N_scans);
    loglik = 0;
    log2pi = log(2 * pi);

    for k = 1:N_scans
        [x_pred, flt] = flt.predict();
        zpred = [x_pred(1); x_pred(3)];
        nu    = measurements(:, k) - zpred;
        Sk    = flt.S;
        detS  = det(Sk);
        if detS > 0 && all(isfinite(nu))
            loglik = loglik - 0.5 * (2*log2pi + log(detS) + nu.' * (Sk \ nu));
        end
        [xk, flt] = flt.update(measurements(:, k));
        est_r(k)  = xk(1);
        est_d(k)  = xk(3);
    end

    err_r = est_r(1:N)' - true_r(1:N);
    err_d = est_d(1:N)' - true_d(1:N);
    out = struct( ...
        'filter',      filter_name, ...
        'scenario',    scenario_name, ...
        'N',           N, ...
        'rmse_r',      sqrt(mean(err_r.^2)), ...
        'rmse_d',      sqrt(mean(err_d.^2)), ...
        'loglik',      loglik, ...
        'loglik_mean', loglik / N, ...
        'est_r',       est_r, ...
        'est_d',       est_d, ...
        'meas',        measurements);

    fprintf('\n==== %s on %s (N=%d) ====\n', filter_name, scenario_name, N);
    fprintf('  RMSE range   : %8.3f km\n', out.rmse_r);
    fprintf('  RMSE doppler : %8.3f Hz\n', out.rmse_d);
    fprintf('  Log-likelihood         : %10.2f\n', out.loglik);
    fprintf('  Log-likelihood / scan  : %10.4f\n', out.loglik_mean);

    % ---- Fig 1: range-time waterfall (max over doppler) ----
    range_time = squeeze(max(ard_stack, [], 2));     % (Nrange x Nscan)
    scans = 1:N_scans;
    figure('Name', sprintf('%s | %s | range waterfall', filter_name, scenario_name), ...
           'Position', [80 80 900 500]);
    imagesc(scans, range_axis_km, 10*log10(range_time + eps)); hold on;
    axis xy; colormap(jet); colorbar;
    xlabel('scan'); ylabel('bistatic range (km)');
    title(sprintf('%s on %s | Range-time (max over doppler)  RMSE_r=%.2f km', ...
                  filter_name, scenario_name, out.rmse_r));
    plot(scans(1:N), true_r(1:N),           'w-', 'LineWidth', 2, 'DisplayName', 'truth');
    plot(scans,      measurements(1, :),     'k.', 'MarkerSize', 8, 'DisplayName', 'meas');
    plot(scans,      est_r,                  'r-', 'LineWidth', 1.5, 'DisplayName', 'estimate');
    legend('truth','meas','estimate','Location','best','TextColor','w','Color',[0 0 0 0.4]);

    % ---- Fig 2: doppler-time waterfall (max over range) ----
    dopp_time = squeeze(max(ard_stack, [], 1));       % (Ndoppler x Nscan)
    figure('Name', sprintf('%s | %s | doppler waterfall', filter_name, scenario_name), ...
           'Position', [80 80 900 500]);
    imagesc(scans, dopp_axis_hz, 10*log10(dopp_time + eps)); hold on;
    axis xy; colormap(jet); colorbar;
    xlabel('scan'); ylabel('bistatic doppler (Hz)');
    title(sprintf('%s on %s | Doppler-time (max over range)  RMSE_d=%.2f Hz', ...
                  filter_name, scenario_name, out.rmse_d));
    plot(scans(1:N), true_d(1:N),           'w-', 'LineWidth', 2);
    plot(scans,      measurements(2, :),     'k.', 'MarkerSize', 8);
    plot(scans,      est_d,                  'r-', 'LineWidth', 1.5);
    legend('truth','meas','estimate','Location','best','TextColor','w','Color',[0 0 0 0.4]);

    % ---- Fig 3: ARD snapshot montage ----
    if nargin < 3 || isempty(snapshot_scans)
        snapshot_scans = round(linspace(2, N_scans-1, 4));
    end
    n_snap = numel(snapshot_scans);
    figure('Name', sprintf('%s | %s | ARD snapshots', filter_name, scenario_name), ...
           'Position', [80 80 350*n_snap 400]);
    for ii = 1:n_snap
        k = snapshot_scans(ii);
        subplot(1, n_snap, ii);
        imagesc(dopp_axis_hz, range_axis_km, 10*log10(ard_stack(:, :, k) + eps));
        axis xy; colormap(jet);
        xlabel('doppler (Hz)'); ylabel('range (km)');
        title(sprintf('scan %d', k));
        hold on;
        plot(est_d(k),           est_r(k),           'r+', 'MarkerSize', 12, 'LineWidth', 2);
        plot(measurements(2, k), measurements(1, k), 'wo', 'MarkerSize', 10, 'LineWidth', 1.5);
        if k <= N
            plot(true_d(k), true_r(k), 'g*', 'MarkerSize', 10, 'LineWidth', 1.5);
        end
        if ii == 1
            legend('estimate (+)', 'meas (o)', 'truth (*)', 'Location', 'southoutside', 'TextColor','w');
        end
    end
    sgtitle(sprintf('%s on %s | ARD snapshots', filter_name, scenario_name), ...
            'FontWeight','bold');
end
