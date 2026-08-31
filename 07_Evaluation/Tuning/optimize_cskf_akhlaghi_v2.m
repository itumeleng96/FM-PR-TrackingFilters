% =============================================================
% optimize_cskf_akhlaghi_v2.m
%
% Full grid tuning of CSKF + Akhlaghi + GNN across all knobs and
% all three paper scenarios. Reports:
%   * top-15 global combos (all scenarios balanced)
%   * per-scenario best combo (optimum tuning for each scenario)
%   * final verification of global best on each scenario
%
% Parameters tuned:
%   alpha      Akhlaghi forgetting factor
%   gamma      chi^2 GNN gate threshold
%   R_scale    measurement-noise nominal scale
%   Qr, Qd     process-noise nominal scales (per block)
% =============================================================

clc; clear;

addpath('07_Evaluation/GroundTruth/', '06_Tracking/Filters/KF/');

scenarios = {'levelFlight','landing','orbit360'};

dt = 1;
kf_std_meas_base  = [4.9038, 0.9985];
base_std_acc      = [0.0048354, 0.0991];

% ---- Grid ----
alpha_grid   = [0.9, 0.95, 0.98];
gamma_grid   = [16, 25];
R_grid       = [0.5, 1.0, 2.0];
Qr_grid      = [0.5, 1.0, 2.0, 5.0];
Qd_grid      = [0.25, 0.5, 1.0, 2.0];

combos = combvec(alpha_grid, gamma_grid, R_grid, Qr_grid, Qd_grid).';
n_combo = size(combos, 1);
fprintf('CSKF+Akhlaghi+GNN v2 sweep: %d combos over %d scenarios.\n\n', ...
        n_combo, numel(scenarios));

% ---- Load cached clusters ----
n_scen = numel(scenarios);
cache = cell(n_scen, 1);
for s = 1:n_scen
    cname = sprintf('clusters_%s.mat', scenarios{s});
    if ~isfile(cname)
        error('Missing %s. Run cache_all_clusters.m first.', cname);
    end
    cache{s} = load(cname);
end

% ---- Sweep ----
% metrics: [rmse_r, rmse_d, smooth_r, smooth_d, mean_nis, mean_nees, n_coast]
metrics = nan(n_combo, n_scen, 7);
t0 = tic;

for c = 1:n_combo
    a   = combos(c, 1);
    g   = combos(c, 2);
    Rs  = combos(c, 3);
    Qr  = combos(c, 4);
    Qd  = combos(c, 5);
    sa  = [base_std_acc(1)*Qr, base_std_acc(2)*Qd];
    rst = kf_std_meas_base * Rs;

    for s = 1:n_scen
        d       = cache{s};
        N       = d.N_scans;
        N_gt    = numel(d.true_r);
        Nc      = min(N, N_gt);
        est_r   = nan(1, N);
        est_d   = nan(1, N);
        nis_arr = nan(1, N);
        nees_arr= nan(1, N);
        n_coast = 0;

        if isempty(d.scans{1})
            first_z = [d.true_r(1); d.true_d(1)];
        else
            first_z = d.scans{1}(:, 1);
        end
        init_state = [first_z(1); 0; first_z(2); 0];

        try
            flt = CSKF(dt, sa, rst(1), rst(2), init_state);
            flt.cs_alpha   = a;
            flt.cs_beta    = a;
            flt.cs_adapt_q = true;

            [~, flt] = flt.predict();
            [xk, flt] = flt.update(first_z);
            est_r(1) = xk(1); est_d(1) = xk(3);

            for k = 2:N
                [x_pred, flt] = flt.predict();
                Sk = flt.S;
                [z_gated, ~, ~] = gated_association(d.scans{k}, x_pred, flt.H, Sk, g);
                if any(isnan(z_gated))
                    n_coast = n_coast + 1;
                    est_r(k) = x_pred(1); est_d(k) = x_pred(3);
                else
                    nu = z_gated - flt.H * x_pred;
                    if det(Sk) > 0
                        nis_arr(k) = nu.' * (Sk \ nu);
                    end
                    [xk, flt] = flt.update(z_gated);
                    est_r(k) = xk(1); est_d(k) = xk(3);
                end
                if k <= N_gt
                    HPH = flt.H * flt.P * flt.H.';
                    ze  = [d.true_r(k); d.true_d(k)] - flt.H * flt.X;
                    if det(HPH) > 0
                        nees_arr(k) = ze.' * (HPH \ ze);
                    end
                end
            end

            er = est_r(1:Nc)' - d.true_r(1:Nc);
            ed = est_d(1:Nc)' - d.true_d(1:Nc);
            metrics(c, s, 1) = sqrt(mean(er.^2));
            metrics(c, s, 2) = sqrt(mean(ed.^2));
            metrics(c, s, 3) = std(diff(est_r(1:Nc)));
            metrics(c, s, 4) = std(diff(est_d(1:Nc)));
            metrics(c, s, 5) = mean(nis_arr, 'omitnan');
            metrics(c, s, 6) = mean(nees_arr, 'omitnan');
            metrics(c, s, 7) = n_coast;
        catch
        end
    end

    if mod(c, 30) == 0
        fprintf('  combo %3d/%d   elapsed %.1fs\n', c, n_combo, toc(t0));
    end
end
fprintf('\nSweep complete (%.1fs).\n\n', toc(t0));

% ---- Composite score ----
d_r = 1.5;  d_d = 1.0;
w_smooth = 2.0;   w_rmse = 1.0;

score = nan(n_combo, 1);
per_scen_score = nan(n_combo, n_scen);
diverged = false(n_combo, 1);
for c = 1:n_combo
    ok = true;
    for s = 1:n_scen
        if isnan(metrics(c,s,1)) || metrics(c,s,1) > 8 || metrics(c,s,2) > 15
            ok = false; break;
        end
    end
    diverged(c) = ~ok;
    for s = 1:n_scen
        if isnan(metrics(c,s,1))
            continue;
        end
        per_scen_score(c, s) = w_rmse *(metrics(c,s,1)/d_r + metrics(c,s,2)/d_d) + ...
                               w_smooth*(metrics(c,s,3)/d_r + metrics(c,s,4)/d_d);
    end
    if ok
        sr = mean(metrics(c,:,1)) / d_r;
        sd = mean(metrics(c,:,2)) / d_d;
        mr = mean(metrics(c,:,3)) / d_r;
        md = mean(metrics(c,:,4)) / d_d;
        score(c) = w_rmse*(sr+sd) + w_smooth*(mr+md);
    end
end

fprintf('Divergence: %d / %d combos\n', sum(diverged), n_combo);
fprintf('Valid:      %d / %d combos\n\n', sum(~diverged), n_combo);

% ---- Top-15 global ----
[~, order] = sort(score);
fprintf('Top 15 GLOBAL tunings (all scenarios balanced):\n');
fprintf('%s\n', repmat('=', 1, 150));
fprintf('%-5s | %6s %6s %5s %5s %5s | RMSE_r/RMSE_d per scenario                | avg smR  smD  NIS  NEES  coast | score\n', ...
        'rank','alpha','gamma','Rs','Qr','Qd');
fprintf('%s\n', repmat('-', 1, 150));
for i = 1:min(15, sum(~diverged))
    c = order(i);
    fprintf('%-5d | %6.2f %6.2f %5.2f %5.2f %5.2f | ', ...
            i, combos(c,1), combos(c,2), combos(c,3), combos(c,4), combos(c,5));
    for s = 1:n_scen
        fprintf('%5.2f/%5.2f  ', metrics(c,s,1), metrics(c,s,2));
    end
    fprintf('| %5.2f %5.2f %5.2f %5.2f %5d | %6.3f\n', ...
            mean(metrics(c,:,3)), mean(metrics(c,:,4)), ...
            mean(metrics(c,:,5),'omitnan'), mean(metrics(c,:,6),'omitnan'), ...
            round(mean(metrics(c,:,7))), score(c));
end

% ---- Per-scenario best ----
fprintf('\nPer-scenario BEST tuning:\n');
fprintf('%s\n', repmat('=', 1, 120));
for s = 1:n_scen
    [~, ord_s] = sort(per_scen_score(:, s));
    best = ord_s(1);
    fprintf('%-14s : alpha=%.2f  gamma=%.2f  Rs=%.2f  Qr=%.2f  Qd=%.2f\n', ...
            scenarios{s}, combos(best,1), combos(best,2), combos(best,3), combos(best,4), combos(best,5));
    fprintf('%14s   RMSE=%.2fkm / %.2fHz    Smooth=%.2fkm / %.2fHz    NIS=%.2f  NEES=%.2f  Coasts=%d\n', ...
            '', metrics(best,s,1), metrics(best,s,2), ...
            metrics(best,s,3), metrics(best,s,4), ...
            metrics(best,s,5), metrics(best,s,6), round(metrics(best,s,7)));
end

% ---- Verification: global best on each scenario ----
best_c = order(1);
fprintf('\n%s\n', repmat('=', 1, 120));
fprintf('GLOBAL BEST verification\n');
fprintf('  alpha=%.2f  gamma=%.2f  Rs=%.2f  Qr=%.2f  Qd=%.2f    score=%.3f\n\n', ...
        combos(best_c,1), combos(best_c,2), combos(best_c,3), ...
        combos(best_c,4), combos(best_c,5), score(best_c));
fprintf('%-14s | %8s %8s %8s %8s %6s %6s %6s\n', ...
        'Scenario','RMSE_r','RMSE_d','sm_r','sm_d','NIS','NEES','coast');
fprintf('%s\n', repmat('-', 1, 75));
for s = 1:n_scen
    fprintf('%-14s | %6.2fkm %6.2fHz %6.2fkm %6.2fHz %6.2f %6.2f %6d\n', ...
            scenarios{s}, metrics(best_c,s,1), metrics(best_c,s,2), ...
            metrics(best_c,s,3), metrics(best_c,s,4), ...
            metrics(best_c,s,5), metrics(best_c,s,6), round(metrics(best_c,s,7)));
end

save('cskf_akhlaghi_optimization_v2.mat', 'metrics', 'combos', 'score', ...
     'per_scen_score', 'order', 'diverged', 'scenarios', ...
     'alpha_grid', 'gamma_grid', 'R_grid', 'Qr_grid', 'Qd_grid');
fprintf('\nSaved cskf_akhlaghi_optimization_v2.mat\n');
