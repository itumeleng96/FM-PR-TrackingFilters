% =============================================================
% run_tuned_all_scenarios.m
%
% Purpose
% -------
% Run each of the four filters (CSKF, CSUKF, CSRGNF, CSPF) with
% their Scheme B2 winning tuning across all four scenarios and
% produce diagnostic plots.
%
% Outputs
% -------
%   figures/tuned_<scenario>.png : per-scenario 2xN_filt panel
%       (row 1 = range estimate vs GT vs measurement,
%        row 2 = doppler estimate vs GT vs measurement)
%   Summary table printed to stdout.
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

if ~exist('figures', 'dir'); mkdir('figures'); end

scenarios = { ...
    'levelFlight', 60; ...
    'landing',     120; ...
    'takeoff',     60; ...
    'orbit360',    120; ...
};

dt = 1;

% -------- Baseline (unscaled) filter params (mirror sweep script) --------
kf_std_meas   = [4.9038, 0.9985];
kf_std_acc0   = [0.0048354, 0.0991];
ukf_std_meas  = [0.9707, 0.79739];
ukf_std_acc0  = [0.0076533, 0.09938];
ukf_alpha     = 1e-4; ukf_kappa = 0; ukf_beta = 2;
pf_std_meas   = [10, 2];
pf_std_acc0   = [1.429, 1.9452];
N_pf          = 3000;
rgnf_std_meas = [2.046, 0.98];
rgnf_std_acc0 = [0.057027, 0.047789];
rgnf_max_iter = 100; rgnf_lambda = 1;

% -------- Scheme B2 winning tuning (from scheme_b2_results.mat) --------
% cs_alpha, cs_beta, Q_scale_r, Q_scale_d
tuning = struct( ...
    'CSKF',   struct('a', 0.90, 'b', 0.90, 'qr', 2.0,  'qd', 1.0), ...
    'CSUKF',  struct('a', 0.60, 'b', 0.90, 'qr', 0.5,  'qd', 5.0), ...
    'CSRGNF', struct('a', 0.90, 'b', 0.90, 'qr', 1.0,  'qd', 10.0), ...
    'CSPF',   struct('a', 0.60, 'b', 0.90, 'qr', 1.0,  'qd', 2.0));

filters = {'CSKF', 'CSUKF', 'CSRGNF', 'CSPF'};
n_filt  = numel(filters);
n_scen  = size(scenarios, 1);

summary = struct();

for s = 1:n_scen
    scen_name = scenarios{s, 1};
    scen_sec  = scenarios{s, 2};
    cache     = sprintf('meas_%s.mat', scen_name);
    if ~exist(cache, 'file')
        error('Missing %s. Run scheme_b_grid_sweep.m first.', cache);
    end
    S = load(cache); d = S.data;
    N = d.N_align;
    t = (1:N).';

    fprintf('\n==== %s (N=%d) ====\n', scen_name, N);

    % Prepare figure: 2 rows (range/doppler), n_filt cols
    fig = figure('Position', [100 100 350*n_filt 500], 'Visible', 'off');

    for f = 1:n_filt
        name = filters{f};
        tp   = tuning.(name);

        initial_state = [d.meas(1,1); 0; d.meas(2,1); 0];
        switch name
            case 'CSKF'
                sa = [kf_std_acc0(1)*tp.qr, kf_std_acc0(2)*tp.qd];
                flt = CSKF(dt, sa, kf_std_meas(1), kf_std_meas(2), initial_state);
            case 'CSUKF'
                sa = [ukf_std_acc0(1)*tp.qr, ukf_std_acc0(2)*tp.qd];
                flt = CSUKF(dt, sa, ukf_std_meas(1), ukf_std_meas(2), initial_state', ukf_alpha, ukf_kappa, ukf_beta);
            case 'CSRGNF'
                sa = [rgnf_std_acc0(1)*tp.qr, rgnf_std_acc0(2)*tp.qd];
                flt = CSRGNF(dt, sa, rgnf_std_meas(1), rgnf_std_meas(2), initial_state, rgnf_max_iter, rgnf_lambda);
            case 'CSPF'
                sa = [pf_std_acc0(1)*tp.qr, pf_std_acc0(2)*tp.qd];
                flt = CSPF(dt, sa, pf_std_meas, initial_state, N_pf);
        end
        flt.cs_alpha = tp.a;
        flt.cs_beta  = tp.b;

        est_r = zeros(1, d.N);
        est_d = zeros(1, d.N);
        for k = 1:d.N
            [~,  flt] = flt.predict();
            [xk, flt] = flt.update(d.meas(:, k));
            est_r(k) = xk(1);
            est_d(k) = xk(3);
        end

        err_r = est_r(1:N)' - d.true_r;
        err_d = est_d(1:N)' - d.true_d;
        rmse_r = sqrt(mean(err_r.^2));
        rmse_d = sqrt(mean(err_d.^2));
        summary.(name).(scen_name) = struct('rmse_r', rmse_r, 'rmse_d', rmse_d);

        fprintf('  %-8s RMSE_r=%7.2f km   RMSE_d=%7.2f Hz\n', name, rmse_r, rmse_d);

        % ---- Range panel ----
        subplot(2, n_filt, f); hold on; grid on;
        plot(t, d.meas(1,1:N), 'k.', 'MarkerSize', 4, 'DisplayName', 'meas');
        plot(t, d.true_r,      'g-',  'LineWidth', 2, 'DisplayName', 'truth');
        plot(t, est_r(1:N),    'b-',  'LineWidth', 1.5, 'DisplayName', 'estimate');
        title(sprintf('%s | range  RMSE=%.2f km', name, rmse_r));
        xlabel('scan'); ylabel('bistatic range (km)');
        if f == 1; legend('Location','best'); end

        % ---- Doppler panel ----
        subplot(2, n_filt, n_filt + f); hold on; grid on;
        plot(t, d.meas(2,1:N), 'k.', 'MarkerSize', 4);
        plot(t, d.true_d,      'g-',  'LineWidth', 2);
        plot(t, est_d(1:N),    'b-',  'LineWidth', 1.5);
        title(sprintf('%s | doppler  RMSE=%.2f Hz', name, rmse_d));
        xlabel('scan'); ylabel('bistatic doppler (Hz)');
    end

    sgtitle(sprintf('Scenario: %s (tuned Scheme B2)', scen_name), ...
            'FontWeight', 'bold', 'FontSize', 12);
    out_png = sprintf('figures/tuned_%s.png', scen_name);
    exportgraphics(fig, out_png, 'Resolution', 130);
    close(fig);
    fprintf('  Saved %s\n', out_png);
end

% -------- Print summary table --------
fprintf('\n\n=============================================================\n');
fprintf('  Summary: RMSE per (filter, scenario) with tuned Scheme B2\n');
fprintf('=============================================================\n');
fprintf('%-8s | ', 'Filter');
for s = 1:n_scen
    fprintf('%-22s', scenarios{s,1});
end
fprintf('\n%s\n', repmat('-', 1, 8 + 3 + 22*n_scen));
for f = 1:n_filt
    name = filters{f};
    fprintf('%-8s | ', name);
    for s = 1:n_scen
        scn = scenarios{s,1};
        m   = summary.(name).(scn);
        fprintf('%6.2f km / %6.1f Hz   ', m.rmse_r, m.rmse_d);
    end
    fprintf('\n');
end

save('tuned_all_scenarios.mat', 'summary', 'tuning');
fprintf('\nSaved to tuned_all_scenarios.mat\n');
