function tune_kf_ukf_orbit()
% TUNE_KF_UKF_ORBIT  Small sweep to find a KF/UKF tuning that survives the
% 360-orbit manoeuvre. 3 seeds x orbit only. Print LL/RMSE/NIS per config.

    addpath('TrackingFilter-CSKF/','TrackingFilter-CSUKF/','groundTruthCalculations/');

    seeds  = [1, 25, 50];
    dt     = 1;
    warmup = 10;
    log2pi = log(2*pi);

    % Baseline tunings (from mc_evaluate_filters.m).
    base.KF  = struct('std_meas',[4.9038,0.9985],'std_acc',[0.0048354,0.0991], ...
                      'alpha',0.98,'gamma',25,'Rs',0.5,'Qr',2.0,'Qd',0.25);
    base.UKF = struct('std_meas',[0.9707,0.79739],'std_acc',[0.0076533,0.09938], ...
                      'alpha',0.98,'gamma',25,'Rs',2.0,'Qr',2.0,'Qd',0.25, ...
                      'ukf_a',0.5,'ukf_k',0,'ukf_b',2);

    % Candidate tunings for orbit. Strategy: raise Qd (Doppler-axis process
    % noise) to keep up with continuous heading changes; widen gamma so
    % turning measurements are not gated out. Keep alpha unchanged.
    cand.KF = {
        struct('label','KF-base'     , 'Qr',2.0,'Qd',0.25 ,'gamma',25);
        struct('label','KF-Q2x'      , 'Qr',2.0,'Qd',0.5  ,'gamma',25);
        struct('label','KF-Q4x'      , 'Qr',2.0,'Qd',1.0  ,'gamma',25);
        struct('label','KF-Q4x-gate' , 'Qr',2.0,'Qd',1.0  ,'gamma',36);
        struct('label','KF-Q8x-gate' , 'Qr',2.0,'Qd',2.0  ,'gamma',36);
    };
    cand.UKF = {
        struct('label','UKF-base'    , 'Qr',2.0,'Qd',0.25 ,'gamma',25);
        struct('label','UKF-Q2x'     , 'Qr',2.0,'Qd',0.5  ,'gamma',25);
        struct('label','UKF-Q4x'     , 'Qr',2.0,'Qd',1.0  ,'gamma',25);
        struct('label','UKF-Q4x-gate', 'Qr',2.0,'Qd',1.0  ,'gamma',36);
        struct('label','UKF-Q8x-gate', 'Qr',2.0,'Qd',2.0  ,'gamma',36);
    };

    fprintf('%-16s | %-4s | %8s %8s %8s %8s %8s\n', ...
        'tuning','seed','LL_r','LL_d','RMSE_r','RMSE_d','NIS');
    fprintf('%s\n', repmat('-', 1, 80));

    for fname = {'KF','UKF'}
        f = fname{1};
        tp0 = base.(f);
        for ci = 1:numel(cand.(f))
            c = cand.(f){ci};
            for r = seeds
                cfile = fullfile('seeds', sprintf('seed_%03d', r), 'clusters_orbit360.mat');
                if ~isfile(cfile); continue; end
                d = load(cfile);
                N = d.N_scans; Nc = min(N, numel(d.true_r));

                if isempty(d.scans{1}); z0 = [d.true_r(1); d.true_d(1)];
                else; z0 = d.scans{1}(:,1); end
                init = [z0(1); 0; z0(2); 0];

                sa  = [tp0.std_acc(1)*c.Qr, tp0.std_acc(2)*c.Qd];
                rst = tp0.std_meas * tp0.Rs;

                switch f
                    case 'KF';  flt = CSKF(dt, sa, rst(1), rst(2), init);
                    case 'UKF'; flt = CSUKF(dt, sa, rst(1), rst(2), init', ...
                                           tp0.ukf_a, tp0.ukf_k, tp0.ukf_b);
                end
                flt.cs_alpha = tp0.alpha;
                if isprop(flt,'cs_adapt_q'); flt.cs_adapt_q = true; end

                est_r = nan(1,N); est_d = nan(1,N);
                ll_r = nan(1,N); ll_d = nan(1,N); nis_a = nan(1,N);

                [~, flt]  = flt.predict();
                [xk, flt] = flt.update(z0);
                est_r(1) = xk(1); est_d(1) = xk(3);

                for k = 2:N
                    [x_pred, flt] = flt.predict();
                    Sk = flt.S;
                    [z_g, ~, ~] = gated_association(d.scans{k}, x_pred, flt.H, Sk, c.gamma);
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
                ss = (warmup+1):numel(ll_r);
                sse = (warmup+1):numel(er);

                fprintf('%-16s | %-4d | %8.3f %8.3f %8.3f %8.3f %8.3f\n', ...
                    c.label, r, ...
                    mean(ll_r(ss),'omitnan'), mean(ll_d(ss),'omitnan'), ...
                    sqrt(mean(er(sse).^2,'omitnan')), sqrt(mean(ed(sse).^2,'omitnan')), ...
                    mean(nis_a(ss),'omitnan'));
            end
            fprintf('\n');
        end
    end
end
