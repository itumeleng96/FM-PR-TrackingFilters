function mc_evaluate_filters(N_MC, warmup)
% MC_EVALUATE_FILTERS  Populate Table filter_perf with Monte Carlo statistics.
%
%   mc_evaluate_filters(N_MC, warmup)
%
% FERS-seeded MC: each seed carries a distinct FERS thermal-noise realisation,
% so ARD/CFAR/mean-shift centroids differ per seed. Filter-quality only
% (single-track, always-associated via Mahalanobis GNN gate); MTT
% track-management metrics are in mtt_mc.m.
%
% For each seed s = 1..N_MC, scenario in {landing, orbit360, takeoff},
% filter in {KF, UKF, RGNF, PF}, loads cached centroids from
% seeds/seed_SSS/clusters_<scen>.mat and runs the tuned filter with GNN
% gating. Metrics per (seed, scen, filt), evaluated over steady-state
% window (scans warmup+1 : end) to exclude filter-initialisation transient:
%   LL_r      : mean per-scan log-likelihood, range component
%   LL_d      : mean per-scan log-likelihood, doppler component
%   RMSE_r    : range RMSE vs ground truth  (km)
%   RMSE_d    : doppler RMSE vs ground truth (Hz)
%   mean_NIS  : mean nu^T S^-1 nu (target 2 for 2D Gaussian innovation)
%
% Aggregates mean+/-std across seeds; saves mc_results.mat.
%
% Defaults: N_MC=100, warmup=10 scans.

    if nargin < 1 || isempty(N_MC),  N_MC = 100; end
    if nargin < 2 || isempty(warmup), warmup = 10; end

    addpath('06_Tracking/Filters/KF/','06_Tracking/Filters/UKF/', ...
            '06_Tracking/Filters/RGNF/','06_Tracking/Filters/PF/', ...
            '07_Evaluation/GroundTruth/');

    scenarios = {'landing','orbit360','takeoff'};
    filters   = {'KF','UKF','RGNF','PF'};   % matches paper Table filter_perf column order
    metrics   = {'LL_r','LL_d','RMSE_r','RMSE_d','mean_NIS'};
    dt = 1;

    % Tuned params (from optimize_all_filters.m best).
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

    % Storage: raw(seed, scen, filter, metric)
    raw = nan(N_MC, numel(scenarios), numel(filters), numel(metrics));

    log2pi = log(2*pi);

    for s = 1:N_MC
        seed_pad = sprintf('%03d', s);
        seed_dir = fullfile('seeds', ['seed_' seed_pad]);
        for si = 1:numel(scenarios)
            scen = scenarios{si};
            cfile = fullfile(seed_dir, sprintf('clusters_%s.mat', scen));
            if ~isfile(cfile)
                warning('Missing %s; skipping.', cfile);
                continue;
            end
            d = load(cfile);
            N = d.N_scans; N_gt = numel(d.true_r); Nc = min(N, N_gt);

            for fi = 1:numel(filters)
                fname = filters{fi};
                tp    = T.(fname);
                sa    = [tp.std_acc(1)*tp.Qr, tp.std_acc(2)*tp.Qd];
                rst   = tp.std_meas * tp.Rs;

                if isempty(d.scans{1})
                    z0 = [d.true_r(1); d.true_d(1)];
                else
                    z0 = d.scans{1}(:,1);
                end
                init = [z0(1); 0; z0(2); 0];

                try
                    switch fname
                        case 'KF'
                            flt = CSKF(dt, sa, rst(1), rst(2), init);
                        case 'UKF'
                            flt = CSUKF(dt, sa, rst(1), rst(2), init', ...
                                        tp.ukf_a, tp.ukf_k, tp.ukf_b);
                        case 'RGNF'
                            flt = CSRGNF(dt, sa, rst(1), rst(2), init, ...
                                         tp.max_iter, tp.lambda);
                        case 'PF'
                            flt = CSPF(dt, sa, rst, init, tp.N);
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
                        [z_g, ~, ~] = gated_association(d.scans{k}, x_pred, ...
                                                        flt.H, Sk, tp.gamma);
                        if any(isnan(z_g))
                            % Coast: no LL contribution, no NIS sample.
                            est_r(k) = x_pred(1); est_d(k) = x_pred(3);
                        else
                            nu = z_g - flt.H * x_pred;
                            if det(Sk) > 0
                                % Per-channel LL using diagonal of Sk.
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

                    % Steady-state window: exclude warmup scans.
                    ss_ll  = (warmup+1):numel(ll_r);
                    ss_nis = (warmup+1):numel(nis_a);
                    ss_er  = (warmup+1):numel(er);

                    raw(s, si, fi, 1) = mean(ll_r(ss_ll),  'omitnan');
                    raw(s, si, fi, 2) = mean(ll_d(ss_ll),  'omitnan');
                    raw(s, si, fi, 3) = sqrt(mean(er(ss_er).^2, 'omitnan'));
                    raw(s, si, fi, 4) = sqrt(mean(ed(ss_er).^2, 'omitnan'));
                    raw(s, si, fi, 5) = mean(nis_a(ss_nis), 'omitnan');
                catch ME
                    warning('seed=%s scen=%s filt=%s failed: %s', ...
                            seed_pad, scen, fname, ME.message);
                end
            end
            fprintf('  seed %s | %-9s | done\n', seed_pad, scen);
        end
    end

    % Aggregate mean+/-std across seeds.
    mu    = squeeze(mean(raw, 1, 'omitnan'));      % scen x filt x metric
    sigma = squeeze(std (raw, 0, 1, 'omitnan'));   % scen x filt x metric

    save('mc_results.mat', 'raw', 'mu', 'sigma', ...
         'scenarios', 'filters', 'metrics', 'N_MC', 'T', 'warmup');
    fprintf('Saved mc_results.mat  (N_MC=%d, warmup=%d scans)\n', N_MC, warmup);

    print_table(mu, sigma, scenarios, filters, metrics);
end

function print_table(mu, sigma, scenarios, filters, metrics)
    fprintf('\n=== Table filter_perf (mean +/- std over seeds) ===\n');
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
