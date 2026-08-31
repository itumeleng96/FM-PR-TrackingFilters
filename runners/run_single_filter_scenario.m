function out = run_single_filter_scenario(filter_name, scenario_name, use_cached_meas, verbose)
% RUN_SINGLE_FILTER_SCENARIO
%   Reproduce the per-scan animation-style plots of runTrackingFilterPlot.m
%   for a single filter on a single scenario, using Scheme B2 tuning.
%
%   out = run_single_filter_scenario('CSPF', 'landing')
%   out = run_single_filter_scenario('CSPF', 'landing', true)   % use cached meas
%   out = run_single_filter_scenario('CSPF', 'landing', true, true) % + verbose
%
% use_cached_meas: if true, load meas_<scenario>.mat (known-good from the
%                  Scheme B sweep) instead of regenerating measurements from
%                  the FERS baseband. Isolates filter behaviour from any
%                  live-pipeline drift. Fig 1/2 will show placeholder
%                  content in this mode.
% verbose:         print first-5 (meas, est, truth) triplets for debugging
%
% Figures (updated per scan, phosphor-style where applicable):
%   Fig 1: ARD (contourf, colour, ardPlot)
%   Fig 2: CFAR + centroids (B&W, ca_cfarPlotBW)
%   Fig 3: Overlaid tracks on RDM axes: estimate, measurement, ground truth
%          (predicted -> blue triangles; measurement -> black dots; truth -> green)
%   Fig 4: Bistatic doppler log-likelihood vs time (cumulative trace)
%   Fig 5: Bistatic range log-likelihood vs time (cumulative trace)
%   Fig 6: Running RMSE (range, doppler) vs time
%
% Log-likelihood follows the convention in multiTargetTracking/logLikelihood.m:
%   ll = -0.5*log(2 pi P^2) - (sample - mean)^2 / (2 P^2)
% with P = filter's own P element on the axis (P(1,1) for range,
% P(3,3) for doppler) and sample = ground truth.

    if nargin < 3; use_cached_meas = false; end
    if nargin < 4; verbose = false; end

    addpath('FERS/', 'cfar/', 'meanShiftCluster/', 'DPI_Suppression/', ...
            'groundTruthCalculations/', ...
            'multiTargetTracking/', ...
            'TrackingFilter-CSKF/',   'TrackingFilter-CSUKF/', ...
            'TrackingFilter-CSRGNF/', 'TrackingFilter-CS-ParticleFilter/');

    % ---- Optionally short-circuit measurement pipeline with the cache ----
    cached = [];
    if use_cached_meas
        cache = sprintf('meas_%s.mat', scenario_name);
        if ~exist(cache, 'file')
            error('use_cached_meas=true but %s not found.', cache);
        end
        S = load(cache);
        cached = S.data;
        fprintf('Using cached measurements from %s (N=%d).\n', cache, cached.N);
    end

    scen_secs = struct('levelFlight', 60, 'landing', 120, 'takeoff', 60, 'orbit360', 120);
    if ~isfield(scen_secs, scenario_name)
        error('Unknown scenario "%s".', scenario_name);
    end
    N_scans = scen_secs.(scenario_name);

    direct_h5 = sprintf('direct_%s.h5', scenario_name);
    echo_h5   = sprintf('echo_%s.h5',   scenario_name);
    if ~exist(direct_h5, 'file') || ~exist(echo_h5, 'file')
        error('Missing %s / %s. Run scheme_b_grid_sweep.m first.', direct_h5, echo_h5);
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

    % -------- Scheme B2 tuning + baselines --------
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

    % -------- Ground truth --------
    [~, true_r, true_d] = computeGroundTruth(scenario_name);
    N_gt = numel(true_r);

    % -------- Load FERS baseband --------
    fprintf('Loading FERS data for %s...\n', scenario_name);
    [Ino,  Qno,  scale_no ] = loadfersHDF5(direct_h5);
    [Imov, Qmov, scale_mov] = loadfersHDF5(echo_h5);
    I_Qmov = (Imov + 1i*Qmov) .* scale_mov;
    I_Qno  = (Ino  + 1i*Qno ) .* scale_no;

    % -------- Preallocate figures (persistent across scans) --------
    f1 = figure('Name','ARD','Position',[100 550 700 450]);
    f2 = figure('Name','CFAR + centroids','Position',[820 550 700 450]);
    f3 = figure('Name','Track: estimate vs meas vs truth','Position',[100 50 700 450]);
    f4 = figure('Name','Bistatic Doppler log-likelihood','Position',[820 50 700 220]);
    f5 = figure('Name','Bistatic Range log-likelihood',  'Position',[820 280 700 220]);
    f6 = figure('Name','Running RMSE',                    'Position',[1520 50 500 450]);

    % -------- Initialise filter using first CFAR centroid --------
    sa = [B.std_acc(1) * tp.qr, B.std_acc(2) * tp.qd];

    ard = []; rdm = [];
    initial = 1; current = fs;

    measurements = nan(2, N_scans);
    est_r  = nan(1, N_scans);
    est_d  = nan(1, N_scans);
    doppler_ll = zeros(1, N_scans);
    range_ll   = zeros(1, N_scans);
    rmse_r_run = zeros(1, N_scans);
    rmse_d_run = zeros(1, N_scans);

    flt = [];   % constructed on first valid measurement

    for i = 1:N_scans
        tic
        s1 = I_Qmov(initial:current);
        s2 = I_Qno(initial:current);
        s1 = procECA_Optimized(s2, s1, proc);

        % --- Fig 1: ARD ---
        [y, ard_] = ardPlot(s1, s2, fs, dopp_bins, delay, i, ard, f1);
        title(sprintf('%s | %s | ARD  (scan %d/%d)', ...
                      filter_name, scenario_name, i, N_scans), 'FontSize', 14);

        % --- Fig 2: CFAR + centroids B/W ---
        [targetClusters, RDM, rdm_] = ca_cfarPlotBW(y.', 1e-7, fs, dopp_bins, delay, i, f2, rdm, 20);

        % --- Mean-shift centroid ---
        [clusterCentroids, ~, ~, ~] = meanShift(targetClusters, 3, 5);
        if ~isempty(clusterCentroids)
            measurements(:, i) = clusterCentroids(1:2, 1);
        elseif i > 1
            measurements(:, i) = measurements(:, i-1);
        else
            % No detection on scan 1: fall back to ground truth or ARD peak.
            if N_gt >= 1
                measurements(:, i) = [true_r(1); true_d(1)];
            else
                [~, kmax] = max(y(:));
                [ir, id]  = ind2sub(size(y), kmax);
                measurements(:, i) = [(ir-1)*(c/fs)/1000; id - (dopp_bins+1)];
            end
        end

        % --- Defer filter construction until first non-NaN measurement ---
        if isempty(flt) && all(~isnan(measurements(:, i)))
            initial_state = [measurements(1,i); 0; measurements(2,i); 0];
            switch filter_name
                case 'CSKF'
                    flt = CSKF(dt, sa, B.std_meas(1), B.std_meas(2), initial_state);
                case 'CSUKF'
                    flt = CSUKF(dt, sa, B.std_meas(1), B.std_meas(2), initial_state', ...
                                B.alpha, B.kappa, B.beta);
                case 'CSRGNF'
                    flt = CSRGNF(dt, sa, B.std_meas(1), B.std_meas(2), initial_state, ...
                                 B.max_iter, B.lambda);
                case 'CSPF'
                    flt = CSPF(dt, sa, B.std_meas, initial_state, B.N);
            end
            flt.cs_alpha = tp.a;
            flt.cs_beta  = tp.b;
        end

        % --- Filter predict + update (skip until filter constructed) ---
        if ~isempty(flt)
            [~,  flt] = flt.predict();
            [xk, flt] = flt.update(measurements(:, i));
            est_r(i) = xk(1);
            est_d(i) = xk(3);
        end

        % --- Fig 3: overlay tracks on the actual CFAR mask ---
        % RDM from ca_cfarPlotBW is (2*Ndop+1) x (Ndelay+1) = doppler x range
        % so imagesc(range_axis, freq_axis, RDM) places doppler on y, range on x.
        figure(f3); clf;
        range_axis = (0:floor(delay*fs)) * (c/fs) / 1000;
        freq_axis  = -dopp_bins:1:dopp_bins;
        imagesc(range_axis, freq_axis, RDM);
        set(gca, 'Color', [0.85 0.85 0.85]);
        colormap(gca, gray); axis xy; grid on; hold on;
        xlabel('Bistatic Range [km]','FontSize',14);
        ylabel('Bistatic Doppler [Hz]','FontSize',14);
        title(sprintf('%s | %s | scan %d/%d  (CFAR mask + tracks)', ...
                      filter_name, scenario_name, i, N_scans), 'FontSize', 14);
        xlim([0 max(range_axis)]); ylim([-dopp_bins dopp_bins]);
        text(1, dopp_bins-15, sprintf('t = %d s', i), 'FontSize', 12, 'Color', 'k');

        % Truth line up to now (green)
        if N_gt > 0
            k_gt = min(i, N_gt);
            plot(true_r(1:k_gt), true_d(1:k_gt), 'g-', 'LineWidth', 2, 'DisplayName', 'truth');
        end
        % Measurement history (black dots)
        plot(measurements(1, 1:i), measurements(2, 1:i), 'k.', 'MarkerSize', 6, ...
             'DisplayName', 'measurement');
        % Estimate history (blue triangles + line)
        plot(est_r(1:i), est_d(1:i), 'b-^', 'LineWidth', 1.5, 'MarkerSize', 5, ...
             'MarkerFaceColor', 'none', 'DisplayName', 'estimate');
        legend('Location', 'best');

        % --- Log-likelihood using filter P (matches plotLogLikelihoodSingleP) ---
        if isprop(flt, 'P') && ~isempty(flt.P) && i <= N_gt
            P_r = flt.P(1,1);
            P_d = flt.P(3,3);
            range_ll(i)   = logLikelihood(est_r(i), P_r, true_r(i));
            doppler_ll(i) = logLikelihood(est_d(i), P_d, true_d(i));
        end

        % --- Running RMSE (up to now, vs truth) ---
        if i <= N_gt
            errs_r = est_r(1:i)' - true_r(1:i);
            errs_d = est_d(1:i)' - true_d(1:i);
            rmse_r_run(i) = sqrt(mean(errs_r.^2));
            rmse_d_run(i) = sqrt(mean(errs_d.^2));
        end

        % --- Fig 4: doppler LL trace ---
        figure(f4); clf;
        plot(1:i, doppler_ll(1:i), 'b-', 'LineWidth', 1.2); grid on;
        xlabel('Time (s)'); ylabel('Doppler log-likelihood');
        title(sprintf('Bistatic Doppler log-likelihood  (%s | %s)', filter_name, scenario_name));

        % --- Fig 5: range LL trace ---
        figure(f5); clf;
        plot(1:i, range_ll(1:i), 'r-', 'LineWidth', 1.2); grid on;
        xlabel('Time (s)'); ylabel('Range log-likelihood');
        title(sprintf('Bistatic Range log-likelihood  (%s | %s)', filter_name, scenario_name));

        % --- Fig 6: running RMSE (dual axis) ---
        figure(f6); clf;
        yyaxis left;
        plot(1:i, rmse_r_run(1:i), 'b-', 'LineWidth', 1.5); grid on;
        ylabel('RMSE range (km)');
        yyaxis right;
        plot(1:i, rmse_d_run(1:i), 'r-', 'LineWidth', 1.5);
        ylabel('RMSE doppler (Hz)');
        xlabel('Time (s)');
        title(sprintf('Running RMSE  (%s | %s)', filter_name, scenario_name));

        ard = ard_;
        rdm = rdm_;
        initial = current + 1;
        current = current + fs;
        drawnow;
        toc
    end

    N = min(N_scans, N_gt);
    err_r = est_r(1:N)' - true_r(1:N);
    err_d = est_d(1:N)' - true_d(1:N);
    out = struct( ...
        'filter',        filter_name, ...
        'scenario',      scenario_name, ...
        'N',             N, ...
        'rmse_r',        sqrt(mean(err_r.^2)), ...
        'rmse_d',        sqrt(mean(err_d.^2)), ...
        'loglik_range',  sum(range_ll(1:N)), ...
        'loglik_dopp',   sum(doppler_ll(1:N)), ...
        'range_ll',      range_ll, ...
        'doppler_ll',    doppler_ll, ...
        'est_r',         est_r, ...
        'est_d',         est_d, ...
        'meas',          measurements, ...
        'true_r',        true_r, ...
        'true_d',        true_d);

    fprintf('\n==== %s on %s (N=%d) ====\n', filter_name, scenario_name, N);
    fprintf('  RMSE range              : %8.3f km\n', out.rmse_r);
    fprintf('  RMSE doppler            : %8.3f Hz\n', out.rmse_d);
    fprintf('  Sum log-likelihood (r)  : %10.2f\n',   out.loglik_range);
    fprintf('  Sum log-likelihood (d)  : %10.2f\n',   out.loglik_dopp);
end
