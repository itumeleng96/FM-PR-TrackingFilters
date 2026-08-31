function out = run_single_filter(filter_name, scenario_name, plot_flag)
% RUN_SINGLE_FILTER  Run one filter on one scenario, report RMSE + log-likelihood.
%
%   out = run_single_filter('CSPF',   'landing')
%   out = run_single_filter('CSKF',   'orbit360', true)   % + plot
%
% filter_name  : 'CSKF' | 'CSUKF' | 'CSRGNF' | 'CSPF'
% scenario_name: 'levelFlight' | 'landing' | 'takeoff' | 'orbit360'
% plot_flag    : optional bool (default false)
%
% Returns struct with fields:
%   rmse_r, rmse_d   - RMSE vs ground truth (km, Hz)
%   loglik           - cumulative measurement log-likelihood
%                      sum_k log N(z_k | H x_k^-, S_k)
%   loglik_mean      - loglik / N (per-scan average)
%   est_r, est_d     - filter estimates
%   true_r, true_d   - ground truth
%
% Uses Scheme B2 winning tuning (from scheme_b2_results.mat).
% Requires cached measurements meas_<scenario>.mat (from scheme_b_grid_sweep.m).

    if nargin < 3; plot_flag = false; end

    addpath('01_FERS/', '04_Detection/', '05_Clustering/', '03_DPI_Cancellation/', ...
            '07_Evaluation/GroundTruth/', ...
            '06_Tracking/Filters/KF/',   '06_Tracking/Filters/UKF/', ...
            '06_Tracking/Filters/RGNF/', '06_Tracking/Filters/PF/');

    % ---- Cached measurements ----
    cache = sprintf('meas_%s.mat', scenario_name);
    if ~exist(cache, 'file')
        error('Missing %s. Run scheme_b_grid_sweep.m first.', cache);
    end
    S = load(cache); d = S.data;
    N = d.N_align;

    % ---- Baselines (mirror sweep script) ----
    dt = 1;
    baselines = struct();
    baselines.CSKF   = struct('std_meas', [4.9038, 0.9985],   'std_acc', [0.0048354, 0.0991]);
    baselines.CSUKF  = struct('std_meas', [0.9707, 0.79739],  'std_acc', [0.0076533, 0.09938], ...
                              'alpha', 1e-4, 'kappa', 0, 'beta', 2);
    baselines.CSRGNF = struct('std_meas', [2.046, 0.98],      'std_acc', [0.057027, 0.047789], ...
                              'max_iter', 100, 'lambda', 1);
    baselines.CSPF   = struct('std_meas', [10, 2],            'std_acc', [1.429, 1.9452], 'N', 3000);

    % ---- Scheme B2 tuning ----
    tuning = struct();
    tuning.CSKF   = struct('a', 0.90, 'b', 0.90, 'qr', 2.0,  'qd', 1.0);
    tuning.CSUKF  = struct('a', 0.60, 'b', 0.90, 'qr', 0.5,  'qd', 5.0);
    tuning.CSRGNF = struct('a', 0.90, 'b', 0.90, 'qr', 1.0,  'qd', 10.0);
    tuning.CSPF   = struct('a', 0.60, 'b', 0.90, 'qr', 1.0,  'qd', 2.0);

    if ~isfield(baselines, filter_name)
        error('Unknown filter "%s". Use CSKF | CSUKF | CSRGNF | CSPF.', filter_name);
    end
    B  = baselines.(filter_name);
    tp = tuning.(filter_name);

    initial_state = [d.meas(1,1); 0; d.meas(2,1); 0];
    sa = [B.std_acc(1) * tp.qr, B.std_acc(2) * tp.qd];

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

    % ---- Loop: run filter, accumulate log-likelihood ----
    % Measurement log-likelihood (Anderson & Moore 1979, Sec 10.5):
    %   log p(z_k | z_{1:k-1}) = -0.5 [ 2 log(2 pi) + log|S_k| + nu_k' S_k^{-1} nu_k ]
    % where nu_k = z_k - H x_k^- and S_k = H P_k^- H' + R (or PF equivalent).
    est_r  = zeros(1, N);
    est_d  = zeros(1, N);
    loglik = 0;
    log2pi = log(2 * pi);

    for k = 1:d.N
        [x_pred, flt] = flt.predict();

        % Predicted measurement (H*x for KF-family, x([1;3]) for PF).
        switch filter_name
            case 'CSUKF'
                zpred = [x_pred(1); x_pred(3)];
            otherwise
                zpred = [x_pred(1); x_pred(3)];
        end
        nu = d.meas(:, k) - zpred;
        Sk = flt.S;

        if k <= N
            % Numerically safe log-det for 2x2 diagonal-ish S.
            detS = det(Sk);
            if detS > 0 && all(isfinite(nu))
                loglik = loglik - 0.5 * (2*log2pi + log(detS) + nu.' * (Sk \ nu));
            end
        end

        [xk, flt] = flt.update(d.meas(:, k));
        if k <= N
            est_r(k) = xk(1);
            est_d(k) = xk(3);
        end
    end

    err_r = est_r(:) - d.true_r(:);
    err_d = est_d(:) - d.true_d(:);
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
        'true_r',      d.true_r, ...
        'true_d',      d.true_d, ...
        'meas',        d.meas(:, 1:N));

    fprintf('\n==== %s on %s (N=%d) ====\n', filter_name, scenario_name, N);
    fprintf('  RMSE range   : %8.3f km\n',  out.rmse_r);
    fprintf('  RMSE doppler : %8.3f Hz\n',  out.rmse_d);
    fprintf('  Log-likelihood         : %10.2f\n', out.loglik);
    fprintf('  Log-likelihood / scan  : %10.4f\n', out.loglik_mean);

    if plot_flag
        t = (1:N).';
        figure('Position', [100 100 900 500]);
        subplot(2,1,1); hold on; grid on;
        plot(t, d.meas(1,1:N), 'k.', 'MarkerSize', 5, 'DisplayName', 'meas');
        plot(t, d.true_r,      'g-', 'LineWidth', 2, 'DisplayName', 'truth');
        plot(t, est_r,         'b-', 'LineWidth', 1.5, 'DisplayName', 'estimate');
        title(sprintf('%s on %s | range RMSE=%.2f km, loglik/scan=%.3f', ...
                      filter_name, scenario_name, out.rmse_r, out.loglik_mean));
        xlabel('scan'); ylabel('bistatic range (km)'); legend('Location','best');

        subplot(2,1,2); hold on; grid on;
        plot(t, d.meas(2,1:N), 'k.', 'MarkerSize', 5);
        plot(t, d.true_d,      'g-', 'LineWidth', 2);
        plot(t, est_d,         'b-', 'LineWidth', 1.5);
        title(sprintf('doppler RMSE=%.2f Hz', out.rmse_d));
        xlabel('scan'); ylabel('bistatic doppler (Hz)');
    end
end
