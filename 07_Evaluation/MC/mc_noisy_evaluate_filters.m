function mc_noisy_evaluate_filters(N_MC, sigma_r_j, sigma_d_j, warmup)
% MC_NOISY_EVALUATE_FILTERS  Injected-noise Monte Carlo evaluation.
%
%   mc_noisy_evaluate_filters(N_MC, sigma_r_j, sigma_d_j, warmup)
%
% For each MC run r = 1..N_MC, the pre-cached mean-shift centroids
% from clusters_<scenario>.mat are perturbed by additive Gaussian noise
%   N(0, diag(sigma_r_j^2, sigma_d_j^2))
% independently per (scan, centroid). The four tuned filters
% (KF, UKF, PF, RGNF) are then run with GNN gating on the perturbed
% centroids and metrics collected.
%
% Metrics per (r, scen, filt), averaged over the steady-state window
% (scans warmup+1 : end) to exclude filter initialisation transient:
%   LL_r      : mean per-scan log-likelihood, range component
%   LL_d      : mean per-scan log-likelihood, doppler component
%   RMSE_r    : range RMSE vs ground truth  (km)
%   RMSE_d    : doppler RMSE vs ground truth (Hz)
%   mean_NIS  : mean d^T S^-1 d
%
% Defaults: N_MC=100, sigma_r_j=0.5 km, sigma_d_j=0.2 Hz, warmup=10 scans.
% Saves mc_results.mat (overwriting).

    if nargin < 1 || isempty(N_MC),      N_MC = 100;    end
    if nargin < 2 || isempty(sigma_r_j), sigma_r_j = 0.5; end   % km
    if nargin < 3 || isempty(sigma_d_j), sigma_d_j = 0.2; end   % Hz
    if nargin < 4 || isempty(warmup),    warmup = 10;   end     % scans

    addpath('06_Tracking/Filters/KF/','06_Tracking/Filters/UKF/', ...
            '06_Tracking/Filters/RGNF/','06_Tracking/Filters/PF/', ...
            '07_Evaluation/GroundTruth/');

    scenarios = {'landing','orbit360','takeoff'};
    filters   = {'KF','UKF','PF','RGNF'};     % paper Table filter_perf column order
    metrics   = {'LL_r','LL_d','RMSE_r','RMSE_d','mean_NIS'};
    dt        = 1;

    % Tuned params (from optimize_all_filters.m, take-off included in sweep).
    T = struct();
    T.KF   = struct('std_meas',[4.9038,0.9985],   'std_acc',[0.0048354,0.0991], ...
                    'alpha',0.98,'gamma',25,'Rs',0.5,'Qr',2.0,'Qd',0.25);
    T.UKF  = struct('std_meas',[0.9707,0.79739],  'std_acc',[0.0076533,0.09938], ...
                    'alpha',0.98,'gamma',25,'Rs',2.0,'Qr',2.0,'Qd',0.25, ...
                    'ukf_a',0.5,'ukf_k',0,'ukf_b',2);
    T.RGNF = struct('std_meas',[2.046,0.98],      'std_acc',[0.057027,0.047789], ...
                    'alpha',0.98,'gamma',25,'Rs',1.0,'Qr',0.5,'Qd',0.25, ...
                    'max_iter',100,'lambda',1);
    T.PF   = struct('std_meas',[10,2],            'std_acc',[1.0,1.0], ...
                    'alpha',0.90,'gamma',16,'Rs',2.0,'Qr',0.5,'Qd',5.0, ...
                    'N',10000);

    % Load baseline centroids once per scenario.
    base = struct();
    for si = 1:numel(scenarios)
        scen  = scenarios{si};
        cfile = sprintf('clusters_%s.mat', scen);
        if ~isfile(cfile)
            error('Missing %s. Run cache_all_clusters.m or copy from seeds/seed_001/.', cfile);
        end
        base.(scen) = load(cfile);
    end

    raw    = nan(N_MC, numel(scenarios), numel(filters), numel(metrics));
    log2pi = log(2*pi);

    for r = 1:N_MC
        rng(r, 'twister');   % reproducible per-run noise

        for si = 1:numel(scenarios)
            scen = scenarios{si};
            d    = base.(scen);
            N    = d.N_scans; N_gt = numel(d.true_r); Nc = min(N, N_gt);

            % Perturb centroids once per MC run (all filters see the same
            % perturbed measurements for fair comparison).
            scans_p = cell(1, N);
            for k = 1:N
                c = d.scans{k};
                if isempty(c), scans_p{k} = []; continue; end
                jitter = [sigma_r_j; sigma_d_j] .* randn(size(c));
                scans_p{k} = c + jitter;
            end

            for fi = 1:numel(filters)
                fname = filters{fi};
                tp    = T.(fname);
                sa    = [tp.std_acc(1)*tp.Qr, tp.std_acc(2)*tp.Qd];
                rst   = tp.std_meas * tp.Rs;

                if isempty(scans_p{1})
                    z0 = [d.true_r(1); d.true_d(1)];
                else
                    z0 = scans_p{1}(:,1);
                end
                init = [z0(1); 0; z0(2); 0];

                try
                    switch fname
                        case 'KF';   flt = CSKF(dt, sa, rst(1), rst(2), init);
                        case 'UKF';  flt = CSUKF(dt, sa, rst(1), rst(2), init', ...
                                                 tp.ukf_a, tp.ukf_k, tp.ukf_b);
                        case 'RGNF'; flt = CSRGNF(dt, sa, rst(1), rst(2), init, ...
                                                  tp.max_iter, tp.lambda);
                        case 'PF';   flt = CSPF(dt, sa, rst, init, tp.N);
                    end
                    flt.cs_alpha = tp.alpha;
                    if isprop(flt,'cs_adapt_q'), flt.cs_adapt_q = true; end

                    est_r = nan(1,N); est_d = nan(1,N);
                    ll_r  = nan(1,N); ll_d  = nan(1,N);
                    nis_a = nan(1,N);

                    [~, flt]  = flt.predict();
                    [xk, flt] = flt.update(z0);
                    est_r(1) = xk(1); est_d(1) = xk(3);

                    for k = 2:N
                        [x_pred, flt] = flt.predict();
                        Sk = flt.S;
                        [z_g, ~, ~] = gated_association(scans_p{k}, x_pred, ...
                                                        flt.H, Sk, tp.gamma);
                        if any(isnan(z_g))
                            est_r(k) = x_pred(1); est_d(k) = x_pred(3);
                        else
                            nu = z_g - flt.H * x_pred;
                            if det(Sk) > 0
                                sr = max(Sk(1,1), 1e-9);
                                sd = max(Sk(2,2), 1e-9);
                                ll_r(k)  = -0.5*(log2pi + log(sr) + nu(1)^2/sr);
                                ll_d(k)  = -0.5*(log2pi + log(sd) + nu(2)^2/sd);
                                nis_a(k) = nu.' * (Sk \ nu);
                            end
                            [xk, flt] = flt.update(z_g);
                            est_r(k) = xk(1); est_d(k) = xk(3);
                        end
                    end

                    er = est_r(1:Nc)' - d.true_r(1:Nc);
                    ed = est_d(1:Nc)' - d.true_d(1:Nc);

                    % Steady-state window: exclude warmup scans to remove
                    % filter-initialisation transient bias.
                    ss_ll  = (warmup+1):numel(ll_r);
                    ss_nis = (warmup+1):numel(nis_a);
                    ss_er  = (warmup+1):numel(er);

                    raw(r, si, fi, 1) = mean(ll_r(ss_ll),  'omitnan');
                    raw(r, si, fi, 2) = mean(ll_d(ss_ll),  'omitnan');
                    raw(r, si, fi, 3) = sqrt(mean(er(ss_er).^2, 'omitnan'));
                    raw(r, si, fi, 4) = sqrt(mean(ed(ss_er).^2, 'omitnan'));
                    raw(r, si, fi, 5) = mean(nis_a(ss_nis), 'omitnan');
                catch ME
                    warning('run=%d scen=%s filt=%s failed: %s', r, scen, fname, ME.message);
                end
            end
        end
        if mod(r, 10) == 0 || r == 1
            fprintf('  MC run %3d/%d\n', r, N_MC);
        end
    end

    mu    = squeeze(mean(raw, 1, 'omitnan'));      % scen x filt x metric
    sigma = squeeze(std (raw, 0, 1, 'omitnan'));   % scen x filt x metric

    save('mc_results.mat', 'raw', 'mu', 'sigma', ...
         'scenarios', 'filters', 'metrics', 'N_MC', 'T', ...
         'sigma_r_j', 'sigma_d_j', 'warmup');
    fprintf('Saved mc_results.mat  (N_MC=%d, sigma_j=[%.2f km, %.2f Hz], warmup=%d scans)\n', ...
            N_MC, sigma_r_j, sigma_d_j, warmup);

    print_table(mu, sigma, scenarios, filters, metrics);
end

function print_table(mu, sigma, scenarios, filters, metrics)
    fprintf('\n=== Table filter_perf (mean +/- std over MC runs) ===\n');
    for si = 1:numel(scenarios)
        fprintf('\n[%s]\n', scenarios{si});
        fprintf('%-10s', 'Metric');
        for fi = 1:numel(filters), fprintf(' %-20s', filters{fi}); end
        fprintf('\n');
        for mi = 1:numel(metrics)
            fprintf('%-10s', metrics{mi});
            for fi = 1:numel(filters)
                fprintf(' %9.3f +/- %6.3f  ', mu(si,fi,mi), sigma(si,fi,mi));
            end
            fprintf('\n');
        end
    end
end
