% =============================================================
% optimize_cskf_akhlaghi.m
%
% Grid search for CSKF + Akhlaghi (Eqs. 11 & 15) + GNN gating.
% Uses cached ALL-cluster data from cache_all_clusters.m so the
% sweep can vary the GNN threshold without re-running ECA/CFAR.
%
% Objective (lower is better):
%   score = 2 * (smooth_r/1.5 + smooth_d/1.0)      <-- smoothness weight
%         +      (RMSE_r/1.5   + RMSE_d/1.0)       <-- accuracy weight
%   averaged over scenarios (levelFlight, landing, orbit360).
% Units: bins-of-measurement (1.5 km range bin, 1.0 Hz doppler bin).
%
% Also enforces a soft floor: skip combos whose per-scenario RMSE_r > 8 km
% or RMSE_d > 15 Hz (blown up filters).
% =============================================================

clc; clear;

addpath('groundTruthCalculations/', 'TrackingFilter-CSKF/');

scenarios = {'levelFlight','landing','orbit360'};

dt = 1;
kf_std_meas   = [4.9038, 0.9985];
base_std_acc  = [0.0048354, 0.0991];

% ---- Grid ----
alpha_grid  = [0.7, 0.8, 0.9, 0.95, 0.98];  % Akhlaghi forgetting factor
gamma_grid  = [9.21, 13.82, 16, 25];         % chi^2 gate (2 dof: 99%, 99.9%, 99.97%, ~5-sigma)
Qr_grid     = [0.5, 1.0, 2.0];               % range Q_nom scale
Qd_grid     = [0.5, 1.0, 2.0];               % doppler Q_nom scale

combos = combvec(alpha_grid, gamma_grid, Qr_grid, Qd_grid).';
n_combo = size(combos, 1);
fprintf('CSKF + Akhlaghi + GNN sweep: %d combos over %d scenarios.\n\n', ...
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
    fprintf('  loaded %s (%d scans)\n', cname, cache{s}.N_scans);
end

% ---- Sweep ----
% metrics: [rmse_r, rmse_d, smooth_r, smooth_d, mean_nis, mean_nees, n_coast]
metrics = nan(n_combo, n_scen, 7);

t0 = tic;
for c = 1:n_combo
    a  = combos(c, 1);
    g  = combos(c, 2);
    Qr = combos(c, 3);
    Qd = combos(c, 4);
    sa = [base_std_acc(1)*Qr, base_std_acc(2)*Qd];

    for s = 1:n_scen
        d       = cache{s};
        N       = d.N_scans;
        N_gt    = numel(d.true_r);
        Nc      = min(N, N_gt);
        scans   = d.scans;
        est_r   = nan(1, N);
        est_d   = nan(1, N);
        nis_arr = nan(1, N);
        nees_arr= nan(1, N);
        n_coast = 0;

        % Bootstrap on strongest centroid of scan 1
        if isempty(scans{1})
            first_z = [d.true_r(1); d.true_d(1)];
        else
            first_z = scans{1}(:, 1);
        end
        init_state = [first_z(1); 0; first_z(2); 0];

        try
            flt = CSKF(dt, sa, kf_std_meas(1), kf_std_meas(2), init_state);
            flt.cs_alpha   = a;
            flt.cs_beta    = a;              % Akhlaghi: single alpha for both
            flt.cs_adapt_q = true;

            [~, flt] = flt.predict();
            [xk, flt] = flt.update(first_z);
            est_r(1) = xk(1); est_d(1) = xk(3);

            for k = 2:N
                [x_pred, flt] = flt.predict();
                Sk = flt.S;
                [z_gated, ~, ~] = gated_association(scans{k}, x_pred, flt.H, Sk, g);
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
            % filter blew up — leave NaN
        end
    end
    if mod(c, 20) == 0
        fprintf('  combo %3d/%d   elapsed %.1fs\n', c, n_combo, toc(t0));
    end
end
fprintf('\nSweep complete (%.1fs).\n\n', toc(t0));

% ---- Composite score with divergence floor ----
d_r = 1.5;  d_d = 1.0;
w_smooth = 2.0;   w_rmse = 1.0;

score = nan(n_combo, 1);
diverged = false(n_combo, 1);
for c = 1:n_combo
    ok = true;
    for s = 1:n_scen
        if isnan(metrics(c,s,1)) || metrics(c,s,1) > 8 || metrics(c,s,2) > 15
            ok = false; break;
        end
    end
    diverged(c) = ~ok;
    if ok
        sr = mean(metrics(c,:,1)) / d_r;
        sd = mean(metrics(c,:,2)) / d_d;
        mr = mean(metrics(c,:,3)) / d_r;
        md = mean(metrics(c,:,4)) / d_d;
        score(c) = w_rmse*(sr+sd) + w_smooth*(mr+md);
    end
end

fprintf('Divergence: %d / %d combos (RMSE > 8km or 15Hz on some scenario)\n', ...
        sum(diverged), n_combo);
fprintf('Valid:      %d / %d combos\n\n', sum(~diverged), n_combo);

[~, order] = sort(score);
fprintf('Top 15 CSKF+Akhlaghi+GNN tunings (composite score, lower better):\n');
fprintf('%s\n', repmat('=', 1, 130));
fprintf('%-5s | %6s %6s %5s %5s | RMSE_r  RMSE_d   sm_r   sm_d  NIS   NEES  coast | score\n', ...
        'rank','alpha','gamma','Qr','Qd');
fprintf('%s\n', repmat('-', 1, 130));
for i = 1:min(15, sum(~diverged))
    c = order(i);
    fprintf('%-5d | %6.2f %6.2f %5.2f %5.2f | ', ...
            i, combos(c,1), combos(c,2), combos(c,3), combos(c,4));
    for s = 1:n_scen
        fprintf('%5.2f/%5.2f  ', metrics(c,s,1), metrics(c,s,2));
    end
    fprintf('| smR=%.2f smD=%.2f NIS=%.2f NEES=%.2f coast=%d | %.3f\n', ...
            mean(metrics(c,:,3)), mean(metrics(c,:,4)), ...
            mean(metrics(c,:,5),'omitnan'), mean(metrics(c,:,6),'omitnan'), ...
            round(mean(metrics(c,:,7))), score(c));
end

% ---- Best combo full breakdown ----
best_c = order(1);
fprintf('\nBEST combo: alpha=%.2f  gamma=%.2f  Qr=%.2f  Qd=%.2f  score=%.3f\n', ...
        combos(best_c,1), combos(best_c,2), combos(best_c,3), combos(best_c,4), score(best_c));
fprintf('%-14s | %8s %8s %8s %8s %6s %6s %6s\n', ...
        'Scenario','RMSE_r','RMSE_d','sm_r','sm_d','NIS','NEES','coast');
fprintf('%s\n', repmat('-', 1, 75));
for s = 1:n_scen
    fprintf('%-14s | %6.2fkm %6.2fHz %6.2fkm %6.2fHz %6.2f %6.2f %6d\n', ...
            scenarios{s}, metrics(best_c,s,1), metrics(best_c,s,2), ...
            metrics(best_c,s,3), metrics(best_c,s,4), ...
            metrics(best_c,s,5), metrics(best_c,s,6), round(metrics(best_c,s,7)));
end

save('cskf_akhlaghi_optimization.mat', 'metrics', 'combos', 'score', 'order', ...
     'diverged', 'scenarios', 'alpha_grid', 'gamma_grid', 'Qr_grid', 'Qd_grid');
fprintf('\nSaved cskf_akhlaghi_optimization.mat\n');
