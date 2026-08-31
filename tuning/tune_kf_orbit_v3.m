function tune_kf_orbit_v3()
% TUNE_KF_ORBIT_V3  Third KF sweep — vary the knobs we haven't tried:
% forgetting factor alpha, per-axis Rs (range vs Doppler independently),
% high-Rs + wide-gate combos, and initial P0.

    addpath('TrackingFilter-CSKF/','groundTruthCalculations/');

    seeds  = [1, 25, 50];
    dt     = 1;
    warmup = 10;
    log2pi = log(2*pi);

    base = struct('std_meas',[4.9038,0.9985],'std_acc',[0.0048354,0.0991], ...
                  'alpha',0.98,'gamma',25,'Rs_r',0.5,'Rs_d',0.5, ...
                  'Qr',2.0,'Qd',0.25);

    % Try: slower forgetting, faster forgetting, decoupled Rs on range vs Doppler,
    % wider gate + higher R (opposite of v2), and no adaptation.
    cand = {
        struct('label','KF-base'        ,'alpha',0.98,'gamma',25,'Rs_r',0.5,'Rs_d',0.5,'Qr',2.0,'Qd',0.25);
        struct('label','KF-alpha0.90'   ,'alpha',0.90,'gamma',25,'Rs_r',0.5,'Rs_d',0.5,'Qr',2.0,'Qd',0.25);
        struct('label','KF-alpha0.80'   ,'alpha',0.80,'gamma',25,'Rs_r',0.5,'Rs_d',0.5,'Qr',2.0,'Qd',0.25);
        struct('label','KF-alpha0.995'  ,'alpha',0.995,'gamma',25,'Rs_r',0.5,'Rs_d',0.5,'Qr',2.0,'Qd',0.25);
        struct('label','KF-RsD2'        ,'alpha',0.98,'gamma',25,'Rs_r',0.5,'Rs_d',2.0,'Qr',2.0,'Qd',0.25);
        struct('label','KF-RsD4'        ,'alpha',0.98,'gamma',25,'Rs_r',0.5,'Rs_d',4.0,'Qr',2.0,'Qd',0.25);
        struct('label','KF-RsR2'        ,'alpha',0.98,'gamma',25,'Rs_r',2.0,'Rs_d',0.5,'Qr',2.0,'Qd',0.25);
        struct('label','KF-Rsboth2-g49' ,'alpha',0.98,'gamma',49,'Rs_r',2.0,'Rs_d',2.0,'Qr',2.0,'Qd',0.25);
        struct('label','KF-Rsboth4-g49' ,'alpha',0.98,'gamma',49,'Rs_r',4.0,'Rs_d',4.0,'Qr',2.0,'Qd',0.25);
        struct('label','KF-noadapt'     ,'alpha',1.00,'gamma',25,'Rs_r',0.5,'Rs_d',0.5,'Qr',2.0,'Qd',0.25);
        struct('label','KF-noadapt-Rs2' ,'alpha',1.00,'gamma',25,'Rs_r',2.0,'Rs_d',2.0,'Qr',2.0,'Qd',0.25);
        struct('label','KF-Qd2-g49'     ,'alpha',0.98,'gamma',49,'Rs_r',0.5,'Rs_d',0.5,'Qr',2.0,'Qd',0.5);
    };

    fprintf('%-20s | %-4s | %8s %8s %8s %8s %8s\n', ...
        'tuning','seed','LL_r','LL_d','RMSE_r','RMSE_d','NIS');
    fprintf('%s\n', repmat('-', 1, 86));

    for ci = 1:numel(cand)
        c = cand{ci};
        for r = seeds
            cfile = fullfile('seeds', sprintf('seed_%03d', r), 'clusters_orbit360.mat');
            if ~isfile(cfile); continue; end
            d = load(cfile);
            N = d.N_scans; Nc = min(N, numel(d.true_r));

            if isempty(d.scans{1}); z0 = [d.true_r(1); d.true_d(1)];
            else; z0 = d.scans{1}(:,1); end
            init = [z0(1); 0; z0(2); 0];

            sa  = [base.std_acc(1)*c.Qr, base.std_acc(2)*c.Qd];
            rst_r = base.std_meas(1) * c.Rs_r;
            rst_d = base.std_meas(2) * c.Rs_d;

            flt = CSKF(dt, sa, rst_r, rst_d, init);
            flt.cs_alpha = c.alpha;
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
            ss = (warmup+1):numel(ll_r); sse = (warmup+1):numel(er);

            fprintf('%-20s | %-4d | %8.3f %8.3f %8.3f %8.3f %8.3f\n', ...
                c.label, r, ...
                mean(ll_r(ss),'omitnan'), mean(ll_d(ss),'omitnan'), ...
                sqrt(mean(er(sse).^2,'omitnan')), sqrt(mean(ed(sse).^2,'omitnan')), ...
                mean(nis_a(ss),'omitnan'));
        end
        fprintf('\n');
    end
end
