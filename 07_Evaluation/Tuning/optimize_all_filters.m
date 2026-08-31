% =============================================================
% optimize_all_filters.m
%
% Full grid tuning of all 4 filters (CSKF, CSUKF, CSRGNF, CSPF) with
% Akhlaghi (2018) continuous R + Q + GNN association. Uses cached
% cluster centroids from cache_all_clusters.m so the sweep is fast.
%
% Reports per-filter global-best tuning and per-scenario metrics.
% Same objective as before: composite score of RMSE + smoothness.
% =============================================================

clc; clear;

addpath('07_Evaluation/GroundTruth/', ...
        '06_Tracking/Filters/KF/','06_Tracking/Filters/UKF/', ...
        '06_Tracking/Filters/RGNF/','06_Tracking/Filters/PF/');

scenarios = {'landing','takeoff','orbit360'};   % matches paper §IV.B scenarios
n_scen    = numel(scenarios);
cache     = cell(n_scen, 1);
for s = 1:n_scen
    cname = sprintf('clusters_%s.mat', scenarios{s});
    if ~isfile(cname)
        error('Missing %s. Run cache_all_clusters.m first.', cname);
    end
    cache{s} = load(cname);
end

dt = 1;

% Baseline nominal params per filter
base = struct();
base.KF     = struct('std_meas',[4.9038, 0.9985],       'std_acc',[0.0048354, 0.0991]);
base.UKF    = struct('std_meas',[0.9707, 0.79739],      'std_acc',[0.0076533, 0.09938], ...
                     'alpha',0.5,'kappa',0,'beta',2);   % widened sigma-spread (was 1e-4 -> UKF degenerated to EKF)
base.RGNF   = struct('std_meas',[2.046, 0.98],          'std_acc',[0.057027, 0.047789], ...
                     'max_iter',100,'lambda',1);
base.PF     = struct('std_meas',[10, 2],                'std_acc',[1.0, 1.0], ...
                     'N',1000);

% ---- Grid (shared across filters where applicable) ----
% Extended Qr/Qd upper end so UKF/PF have room to track fast take-off.
alpha_grid = [0.9, 0.95, 0.98];
gamma_grid = [16, 25];
R_grid     = [0.5, 1.0, 2.0];
Qr_grid    = [0.5, 1.0, 2.0, 5.0, 10.0];
Qd_grid    = [0.25, 0.5, 1.0, 2.0, 5.0];

combos = combvec(alpha_grid, gamma_grid, R_grid, Qr_grid, Qd_grid).';
n_combo = size(combos, 1);

filters = {'KF','UKF','RGNF','PF'};

results = struct();
fprintf('Total combos per filter: %d, scenarios: %d, filters: %d\n', ...
        n_combo, n_scen, numel(filters));

for f = 1:numel(filters)
    fname = filters{f};
    fprintf('\n============ Tuning %s ============\n', fname);
    metrics = nan(n_combo, n_scen, 7);
    t0 = tic;

    for c = 1:n_combo
        a  = combos(c, 1);
        g  = combos(c, 2);
        Rs = combos(c, 3);
        Qr = combos(c, 4);
        Qd = combos(c, 5);

        for s = 1:n_scen
            d   = cache{s};
            N   = d.N_scans; N_gt = numel(d.true_r); Nc = min(N,N_gt);
            est_r = nan(1,N); est_d = nan(1,N);
            nis_arr = nan(1,N); nees_arr = nan(1,N); n_coast = 0;

            if isempty(d.scans{1}), first_z = [d.true_r(1); d.true_d(1)];
            else, first_z = d.scans{1}(:,1); end

            try
                flt = local_make_filter(fname, base, dt, Rs, Qr, Qd, first_z);
                flt.cs_alpha   = a;
                if isprop(flt,'cs_adapt_q'), flt.cs_adapt_q = true; end

                [~, flt] = flt.predict();
                [xk, flt] = flt.update(first_z);
                est_r(1) = xk(1); est_d(1) = xk(3);

                for k = 2:N
                    [x_pred, flt] = flt.predict();
                    Sk = flt.S;
                    [z_gated, ~, ~] = gated_association(d.scans{k}, x_pred, ...
                                                         flt.H, Sk, g);
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
                        x_col = flt.X(:);
                        ze    = [d.true_r(k); d.true_d(k)] - flt.H * x_col;
                        if det(HPH) > 0
                            nees_arr(k) = ze.' * (HPH \ ze);
                        end
                    end
                end

                er = est_r(1:Nc)' - d.true_r(1:Nc);
                ed = est_d(1:Nc)' - d.true_d(1:Nc);
                metrics(c,s,1) = sqrt(mean(er.^2));
                metrics(c,s,2) = sqrt(mean(ed.^2));
                metrics(c,s,3) = std(diff(est_r(1:Nc)));
                metrics(c,s,4) = std(diff(est_d(1:Nc)));
                metrics(c,s,5) = mean(nis_arr,'omitnan');
                metrics(c,s,6) = mean(nees_arr,'omitnan');
                metrics(c,s,7) = n_coast;
            catch
            end
        end
        if mod(c, 30) == 0
            fprintf('  combo %3d/%d  elapsed %.1fs\n', c, n_combo, toc(t0));
        end
    end

    % ---- Composite score ----
    d_r = 1.5; d_d = 1.0;
    w_smooth = 2.0; w_rmse = 1.0;
    score = nan(n_combo, 1); diverged = false(n_combo, 1);
    for c = 1:n_combo
        ok = true;
        for s = 1:n_scen
            if isnan(metrics(c,s,1)) || metrics(c,s,1) > 8 || metrics(c,s,2) > 15
                ok = false; break;
            end
        end
        diverged(c) = ~ok;
        if ok
            sr = mean(metrics(c,:,1))/d_r; sd = mean(metrics(c,:,2))/d_d;
            mr = mean(metrics(c,:,3))/d_r; md = mean(metrics(c,:,4))/d_d;
            score(c) = w_rmse*(sr+sd) + w_smooth*(mr+md);
        end
    end
    [~, order] = sort(score);
    best_c = order(1);

    fprintf('\n%s Global-Best Tuning:\n', fname);
    fprintf('  alpha=%.2f  gamma=%.2f  Rs=%.2f  Qr=%.2f  Qd=%.2f    score=%.3f\n', ...
            combos(best_c,1),combos(best_c,2),combos(best_c,3), ...
            combos(best_c,4),combos(best_c,5), score(best_c));
    fprintf('%-14s | %8s %8s %8s %8s %6s %6s %6s\n', ...
            'Scenario','RMSE_r','RMSE_d','sm_r','sm_d','NIS','NEES','coast');
    fprintf('%s\n', repmat('-',1,75));
    for s = 1:n_scen
        fprintf('%-14s | %6.2fkm %6.2fHz %6.2fkm %6.2fHz %6.2f %6.2f %6d\n', ...
                scenarios{s}, metrics(best_c,s,1), metrics(best_c,s,2), ...
                metrics(best_c,s,3), metrics(best_c,s,4), ...
                metrics(best_c,s,5), metrics(best_c,s,6), round(metrics(best_c,s,7)));
    end
    results.(fname).metrics  = metrics;
    results.(fname).score    = score;
    results.(fname).order    = order;
    results.(fname).best     = combos(best_c,:);
    results.(fname).combos   = combos;
end

save('all_filters_optimization.mat','results','scenarios', ...
     'alpha_grid','gamma_grid','R_grid','Qr_grid','Qd_grid');
fprintf('\nSaved all_filters_optimization.mat\n');


% ---- Local helper to construct each filter with correct signature ----
function flt = local_make_filter(fname, base, dt, Rs, Qr, Qd, first_z)
    init = [first_z(1); 0; first_z(2); 0];
    switch fname
        case 'KF'
            sa  = [base.KF.std_acc(1)*Qr, base.KF.std_acc(2)*Qd];
            rst = base.KF.std_meas * Rs;
            flt = CSKF(dt, sa, rst(1), rst(2), init);
        case 'UKF'
            sa  = [base.UKF.std_acc(1)*Qr, base.UKF.std_acc(2)*Qd];
            rst = base.UKF.std_meas * Rs;
            flt = CSUKF(dt, sa, rst(1), rst(2), init', ...
                        base.UKF.alpha, base.UKF.kappa, base.UKF.beta);
        case 'RGNF'
            sa  = [base.RGNF.std_acc(1)*Qr, base.RGNF.std_acc(2)*Qd];
            rst = base.RGNF.std_meas * Rs;
            flt = CSRGNF(dt, sa, rst(1), rst(2), init, ...
                         base.RGNF.max_iter, base.RGNF.lambda);
        case 'PF'
            sa  = [base.PF.std_acc(1)*Qr, base.PF.std_acc(2)*Qd];
            rst = base.PF.std_meas * Rs;
            flt = CSPF(dt, sa, rst, init, base.PF.N);   % 4x1 [r; 0; d; 0]
    end
end
