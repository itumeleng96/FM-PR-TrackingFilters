function mc_kf_ukf_orbit(N_MC, warmup)
% MC_KF_UKF_ORBIT  Targeted re-run: KF + UKF, orbit scenario only, 100 seeds.
%
% Applies retuned UKF (Qd=1.0 vs previous 0.25) to fix orbit divergence.
% KF unchanged (baseline is best-available).
% Splices new KF/UKF orbit rows into existing mc_results.mat.

    if nargin < 1 || isempty(N_MC),   N_MC   = 100; end
    if nargin < 2 || isempty(warmup), warmup = 10;  end

    addpath('TrackingFilter-CSKF/','TrackingFilter-CSUKF/','groundTruthCalculations/');

    filters = {'KF','UKF'};
    dt = 1; log2pi = log(2*pi);

    T.KF  = struct('std_meas',[4.9038,0.9985],  'std_acc',[0.0048354,0.0991], ...
                   'alpha',0.98,'gamma',49,'Rs',0.5,'Qr',2.0,'Qd',0.5);
    T.UKF = struct('std_meas',[0.9707,0.79739], 'std_acc',[0.0076533,0.09938], ...
                   'alpha',0.98,'gamma',25,'Rs',2.0,'Qr',2.0,'Qd',1.0, ...
                   'ukf_a',0.5,'ukf_k',0,'ukf_b',2);

    % raw(seed, filter, metric): LL_r, LL_d, RMSE_r, RMSE_d, mean_NIS
    raw = nan(N_MC, numel(filters), 5);

    for s = 1:N_MC
        cfile = fullfile('seeds', sprintf('seed_%03d', s), 'clusters_orbit360.mat');
        if ~isfile(cfile); warning('Missing %s; skipping.', cfile); continue; end
        d = load(cfile);
        N = d.N_scans; Nc = min(N, numel(d.true_r));

        if isempty(d.scans{1}); z0 = [d.true_r(1); d.true_d(1)];
        else; z0 = d.scans{1}(:,1); end
        init = [z0(1); 0; z0(2); 0];

        for fi = 1:numel(filters)
            f = filters{fi}; tp = T.(f);
            sa  = [tp.std_acc(1)*tp.Qr, tp.std_acc(2)*tp.Qd];
            rst = tp.std_meas * tp.Rs;

            try
                switch f
                    case 'KF';  flt = CSKF(dt, sa, rst(1), rst(2), init);
                    case 'UKF'; flt = CSUKF(dt, sa, rst(1), rst(2), init', ...
                                            tp.ukf_a, tp.ukf_k, tp.ukf_b);
                end
                flt.cs_alpha = tp.alpha;
                if isprop(flt,'cs_adapt_q'); flt.cs_adapt_q = true; end

                est_r = nan(1,N); est_d = nan(1,N);
                ll_r = nan(1,N); ll_d = nan(1,N); nis_a = nan(1,N);

                [~, flt]  = flt.predict();
                [xk, flt] = flt.update(z0);
                est_r(1) = xk(1); est_d(1) = xk(3);

                for k = 2:N
                    [x_pred, flt] = flt.predict();
                    Sk = flt.S;
                    [z_g, ~, ~] = gated_association(d.scans{k}, x_pred, flt.H, Sk, tp.gamma);
                    if any(isnan(z_g))
                        est_r(k) = x_pred(1); est_d(k) = x_pred(3);
                    else
                        nu = z_g - flt.H * x_pred;
                        if det(Sk) > 0
                            sr = max(Sk(1,1),1e-9); sd = max(Sk(2,2),1e-9);
                            ll_r(k) = -0.5*(log2pi + log(sr) + nu(1)^2/sr);
                            ll_d(k) = -0.5*(log2pi + log(sd) + nu(2)^2/sd);
                            nis_a(k) = nu.' * (Sk \ nu);
                        end
                        [xk, flt] = flt.update(z_g);
                        est_r(k) = xk(1); est_d(k) = xk(3);
                    end
                end

                er = est_r(1:Nc)' - d.true_r(1:Nc);
                ed = est_d(1:Nc)' - d.true_d(1:Nc);
                ss = (warmup+1):numel(ll_r); sse = (warmup+1):numel(er);

                raw(s, fi, 1) = mean(ll_r(ss), 'omitnan');
                raw(s, fi, 2) = mean(ll_d(ss), 'omitnan');
                raw(s, fi, 3) = sqrt(mean(er(sse).^2, 'omitnan'));
                raw(s, fi, 4) = sqrt(mean(ed(sse).^2, 'omitnan'));
                raw(s, fi, 5) = mean(nis_a(ss),'omitnan');
            catch ME
                warning('seed=%03d filt=%s failed: %s', s, f, ME.message);
            end
        end
        if mod(s,10)==0 || s==1
            fprintf('  seed %03d/%d done\n', s, N_MC);
        end
    end

    mu    = squeeze(mean(raw, 1, 'omitnan'));
    sigma = squeeze(std (raw, 0, 1, 'omitnan'));

    save('mc_kf_ukf_orbit.mat', 'raw', 'mu', 'sigma', 'filters', 'N_MC', 'T', 'warmup');

    fprintf('\n=== Orbit (KF & UKF retuned; UKF Qd=1.0) — 100 FERS seeded MC ===\n');
    metrics = {'LL_r','LL_d','RMSE_r','RMSE_d','NIS'};
    fprintf('%-8s', 'Metric');
    for fi=1:numel(filters), fprintf(' %-20s', filters{fi}); end
    fprintf('\n');
    for mi = 1:5
        fprintf('%-8s', metrics{mi});
        for fi=1:numel(filters), fprintf(' %9.3f +/- %6.3f  ', mu(fi,mi), sigma(fi,mi)); end
        fprintf('\n');
    end

    % Splice into existing mc_results.mat if present.
    if isfile('mc_results.mat')
        M = load('mc_results.mat');
        orbit_idx = find(strcmp(M.scenarios,'orbit360'), 1);
        kf_idx    = find(strcmp(M.filters,'KF'),  1);
        ukf_idx   = find(strcmp(M.filters,'UKF'), 1);
        if ~isempty(orbit_idx) && ~isempty(kf_idx) && ~isempty(ukf_idx)
            M.raw(:, orbit_idx, kf_idx,  :) = raw(:, 1, :);
            M.raw(:, orbit_idx, ukf_idx, :) = raw(:, 2, :);
            M.mu    = squeeze(mean(M.raw, 1, 'omitnan'));
            M.sigma = squeeze(std (M.raw, 0, 1, 'omitnan'));
            save('mc_results.mat', '-struct', 'M');
            fprintf('\nSpliced KF+UKF orbit rows into mc_results.mat\n');
        else
            warning('Could not locate orbit / KF / UKF indices in mc_results.mat');
        end
    end
end
