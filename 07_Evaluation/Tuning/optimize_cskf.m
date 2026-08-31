% =============================================================
% optimize_cskf.m
%
% Targeted grid search for CSKF tuning that minimises a composite
% objective combining RMSE AND smoothness (std of scan-to-scan
% first-difference of the estimate). Weighted for FM PR emphasis:
% doppler matters more than range; smoothness matters as much as
% RMSE (per axis).
%
% Composite score (lower is better):
%   score = w_rmse_r * RMSE_r/d_r + w_rmse_d * RMSE_d/d_d + ...
%           w_sm_r   * smooth_r/d_r + w_sm_d   * smooth_d/d_d
%
% Units: bins-of-measurement-resolution so range and doppler are
% directly comparable. d_r = 1.5 km, d_d = 1.0 Hz.
%
% Averaged across the paper's three scenarios (levelFlight, landing,
% orbit360). Takeoff excluded (not in paper text).
% =============================================================

clc; clear;

addpath('07_Evaluation/GroundTruth/', '06_Tracking/Filters/KF/');

scenarios = {'levelFlight','landing','orbit360'};

dt = 1;
kf_std_meas   = [4.9038, 0.9985];
base_std_acc  = [0.0048354, 0.0991];

% ---- Grid (targeted at smoothness-relevant params) ----
Qr_grid    = [0.5, 1.0, 2.0, 5.0];        % range process-noise scale
Qd_grid    = [0.5, 1.0, 2.0, 5.0];        % doppler process-noise scale
alpha_grid = [0.7, 0.9];                  % R adapt rate on trigger
beta_grid  = [0.90, 0.95, 0.98];          % R fading-memory rate
qamax_grid = [2.0, 5.0];                  % Q adaptation cap
qbeta_grid = [0.95, 0.98];                % Q fading-memory rate

combos = combvec(Qr_grid, Qd_grid, alpha_grid, beta_grid, qamax_grid, qbeta_grid).';
n_combo = size(combos, 1);
fprintf('CSKF grid: %d combos (Qr=%d x Qd=%d x a=%d x b=%d x qamax=%d x qbeta=%d).\n', ...
        n_combo, numel(Qr_grid), numel(Qd_grid), numel(alpha_grid), ...
        numel(beta_grid), numel(qamax_grid), numel(qbeta_grid));

% ---- Load caches ----
n_scen = numel(scenarios);
scen_data = cell(n_scen, 1);
for s = 1:n_scen
    cache = sprintf('meas_%s.mat', scenarios{s});
    S = load(cache); scen_data{s} = S.data;
end

% ---- Composite objective weights ----
w_rmse_r = 1.0;   w_rmse_d = 1.0;
w_sm_r   = 1.0;   w_sm_d   = 1.0;
d_r = 1.5;        d_d = 1.0;   % FM PR bin sizes

% ---- Sweep ----
% results: n_combo x n_scen x 4  (rmse_r, rmse_d, smooth_r, smooth_d)
results = nan(n_combo, n_scen, 4);

for c = 1:n_combo
    Qr    = combos(c, 1);   Qd    = combos(c, 2);
    a     = combos(c, 3);   b     = combos(c, 4);
    qamax = combos(c, 5);   qbeta = combos(c, 6);
    sa    = [base_std_acc(1)*Qr, base_std_acc(2)*Qd];

    for s = 1:n_scen
        d = scen_data{s};
        try
            init_state = [d.meas(1,1); 0; d.meas(2,1); 0];
            flt = CSKF(dt, sa, kf_std_meas(1), kf_std_meas(2), init_state);
            flt.cs_alpha       = a;
            flt.cs_beta        = b;
            flt.cs_adapt_q     = true;
            flt.cs_q_alpha_max = qamax;
            flt.cs_q_beta      = qbeta;

            est_r = zeros(1, d.N); est_d = zeros(1, d.N);
            for k = 1:d.N
                [~,  flt] = flt.predict();
                [xk, flt] = flt.update(d.meas(:,k));
                est_r(k) = xk(1); est_d(k) = xk(3);
            end
            N = d.N_align;
            er = est_r(1:N)' - d.true_r;
            ed = est_d(1:N)' - d.true_d;
            results(c, s, 1) = sqrt(mean(er.^2));
            results(c, s, 2) = sqrt(mean(ed.^2));
            results(c, s, 3) = std(diff(est_r(1:N)));
            results(c, s, 4) = std(diff(est_d(1:N)));
        catch
        end
    end
    if mod(c, 50) == 0
        fprintf('  combo %d/%d\n', c, n_combo);
    end
end

% ---- Composite score per combo (averaged across scenarios) ----
score = zeros(n_combo, 1);
for c = 1:n_combo
    s_r  = mean(results(c, :, 1), 'omitnan') / d_r;
    s_d  = mean(results(c, :, 2), 'omitnan') / d_d;
    sm_r = mean(results(c, :, 3), 'omitnan') / d_r;
    sm_d = mean(results(c, :, 4), 'omitnan') / d_d;
    score(c) = w_rmse_r*s_r + w_rmse_d*s_d + w_sm_r*sm_r + w_sm_d*sm_d;
end

% ---- Top-10 combos ----
[~, order] = sort(score);
fprintf('\n\nTop 10 CSKF tunings (composite score, lower better):\n');
fprintf('%s\n', repmat('=', 1, 110));
fprintf('%-5s | %5s %5s %5s %5s %6s %6s | %-30s | %s\n', ...
        'rank','Qr','Qd','alph','beta','qamax','qbeta', ...
        'per-scenario RMSE (km/Hz)', 'score');
fprintf('%s\n', repmat('-', 1, 110));
for i = 1:min(10, n_combo)
    c = order(i);
    fprintf('%-5d | %5.2f %5.2f %5.2f %5.2f %6.2f %6.2f | ', ...
            i, combos(c,1), combos(c,2), combos(c,3), combos(c,4), combos(c,5), combos(c,6));
    for s = 1:n_scen
        fprintf('%5.1f/%5.1f  ', results(c,s,1), results(c,s,2));
    end
    fprintf('| %6.3f\n', score(c));
end
fprintf('%s\n', repmat('=', 1, 110));

% ---- Best combo full breakdown ----
best_c = order(1);
fprintf('\nBEST combo (rank 1):  Qr=%.2f  Qd=%.2f  alpha=%.2f  beta=%.2f  qamax=%.2f  qbeta=%.2f\n', ...
        combos(best_c,1),combos(best_c,2),combos(best_c,3),combos(best_c,4),combos(best_c,5),combos(best_c,6));
fprintf('%-14s | %8s %8s %8s %8s\n', ...
        'Scenario','RMSE_r','RMSE_d','smooth_r','smooth_d');
fprintf('%s\n', repmat('-', 1, 55));
for s = 1:n_scen
    fprintf('%-14s | %6.2fkm %6.1fHz %6.2fkm %6.1fHz\n', ...
            scenarios{s}, results(best_c,s,1), results(best_c,s,2), ...
            results(best_c,s,3), results(best_c,s,4));
end

save('cskf_optimization.mat', 'results','combos','score','order','scenarios');
fprintf('\nSaved cskf_optimization.mat\n');
