% =============================================================
% scheme_b2_grid_sweep.m
%
% Purpose
% -------
% Extends Scheme B by decoupling Q-scaling into per-axis factors
% (Q_scale_r for the range chain, Q_scale_d for the doppler chain).
% This lets filters that maneuver hard in doppler (landing, orbit360)
% inflate only the axis that needs it, without smearing range.
%
% Reuses cached measurements (meas_<scenario>.mat) from
% scheme_b_grid_sweep.m so this stage is tracker-runs only.
%
% Grid
% ----
%   cs_alpha   in {0.6, 0.9}
%   cs_beta    in {0.90, 0.95}
%   Q_scale_r  in {0.5, 1, 2, 5}
%   Q_scale_d  in {0.5, 1, 2, 5, 10}
% Total: 2 x 2 x 4 x 5 = 80 combos per filter.
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

scenarios = { ...
    'levelFlight', 60; ...
    'landing',     120; ...
    'takeoff',     60; ...
    'orbit360',    120; ...
};

dt = 1;
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

% -------- Grid --------
alpha_grid    = [0.6, 0.9];
beta_grid     = [0.90, 0.95];
qs_r_grid     = [0.5, 1, 2, 5];
qs_d_grid     = [0.5, 1, 2, 5, 10];
combos        = combvec(alpha_grid, beta_grid, qs_r_grid, qs_d_grid).';
n_combo       = size(combos, 1);
fprintf('Scheme B2: %d combos per filter (alpha=%d x beta=%d x Qr=%d x Qd=%d).\n', ...
        n_combo, numel(alpha_grid), numel(beta_grid), numel(qs_r_grid), numel(qs_d_grid));

% -------- Load cached measurements --------
n_scen = size(scenarios, 1);
scen_data = cell(n_scen, 1);
for s = 1:n_scen
    scen_name = scenarios{s, 1};
    cache     = sprintf('meas_%s.mat', scen_name);
    if ~exist(cache, 'file')
        error('Missing measurement cache %s. Run scheme_b_grid_sweep.m first.', cache);
    end
    S = load(cache);
    scen_data{s} = S.data;
end
fprintf('Loaded cached measurements for %d scenarios.\n', n_scen);

% -------- Per-filter sweep --------
filters = {'CSKF', 'CSUKF', 'CSRGNF', 'CSPF'};
n_filt  = numel(filters);

results = cell(n_filt, 1);
for f = 1:n_filt
    results{f} = nan(n_combo, n_scen, 2);
end

for f = 1:n_filt
    name = filters{f};
    fprintf('\n---- Sweeping %s (%d combos x %d scenarios) ----\n', ...
            name, n_combo, n_scen);
    for c = 1:n_combo
        a    = combos(c, 1);
        b    = combos(c, 2);
        qs_r = combos(c, 3);
        qs_d = combos(c, 4);
        for s = 1:n_scen
            d = scen_data{s};
            initial_state = [d.meas(1,1); 0; d.meas(2,1); 0];
            try
                switch name
                    case 'CSKF'
                        std_acc = [kf_std_acc0(1)*qs_r, kf_std_acc0(2)*qs_d];
                        flt = CSKF(dt, std_acc, kf_std_meas(1), kf_std_meas(2), initial_state);
                    case 'CSUKF'
                        std_acc = [ukf_std_acc0(1)*qs_r, ukf_std_acc0(2)*qs_d];
                        flt = CSUKF(dt, std_acc, ukf_std_meas(1), ukf_std_meas(2), initial_state', ukf_alpha, ukf_kappa, ukf_beta);
                    case 'CSRGNF'
                        std_acc = [rgnf_std_acc0(1)*qs_r, rgnf_std_acc0(2)*qs_d];
                        flt = CSRGNF(dt, std_acc, rgnf_std_meas(1), rgnf_std_meas(2), initial_state, rgnf_max_iter, rgnf_lambda);
                    case 'CSPF'
                        std_acc = [pf_std_acc0(1)*qs_r, pf_std_acc0(2)*qs_d];
                        flt = CSPF(dt, std_acc, pf_std_meas, initial_state, N_pf);
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
        if mod(c, 20) == 0
            fprintf('    ...combo %d/%d\n', c, n_combo);
        end
    end
end

% -------- Best combo per filter (normalised score) --------
best = struct();
for f = 1:n_filt
    R  = results{f};
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
    score = mean(mean(Rn, 3, 'omitnan'), 2, 'omitnan');
    [best_score, best_c] = min(score);

    best(f).filter    = filters{f};
    best(f).combo_idx = best_c;
    best(f).cs_alpha  = combos(best_c, 1);
    best(f).cs_beta   = combos(best_c, 2);
    best(f).q_scale_r = combos(best_c, 3);
    best(f).q_scale_d = combos(best_c, 4);
    best(f).score     = best_score;
    best(f).rmse_r    = squeeze(R(best_c, :, 1));
    best(f).rmse_d    = squeeze(R(best_c, :, 2));
end

% -------- Report --------
fprintf('\n\n=============================================================\n');
fprintf('  Scheme B2 grid-sweep (Q per-axis) results\n');
fprintf('=============================================================\n\n');

fprintf('Best tuning per filter (normalised-RMSE score, lower better):\n');
fprintf('%-8s | %8s %8s %8s %8s | %8s\n', ...
        'Filter','cs_alpha','cs_beta','Q_scale_r','Q_scale_d','score');
fprintf('%s\n', repmat('-', 1, 72));
for f = 1:n_filt
    fprintf('%-8s | %8.2f %8.2f %9.2f %9.2f | %8.4f\n', ...
            best(f).filter, best(f).cs_alpha, best(f).cs_beta, ...
            best(f).q_scale_r, best(f).q_scale_d, best(f).score);
end

fprintf('\nPer-scenario RMSE at best tuning (km / Hz):\n');
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

save('scheme_b2_results.mat', 'results', 'best', 'combos', 'scenarios', 'filters');
fprintf('\nSaved to scheme_b2_results.mat\n');
