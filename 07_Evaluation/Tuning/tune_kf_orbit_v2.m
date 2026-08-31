function tune_kf_orbit_v2()
% TUNE_KF_ORBIT_V2  Second KF sweep — go the OTHER way: shrink Qd, tighten
% the gate, and vary Rs to see if a stiffer filter kills the seed-25/50
% divergence without exploding.

    addpath('06_Tracking/Filters/KF/','07_Evaluation/GroundTruth/');

    seeds  = [1, 25, 50];
    dt     = 1;
    warmup = 10;
    log2pi = log(2*pi);

    base = struct('std_meas',[4.9038,0.9985],'std_acc',[0.0048354,0.0991], ...
                  'alpha',0.98,'gamma',25,'Rs',0.5,'Qr',2.0,'Qd',0.25);

    % Candidate directions: shrink Qd, tighten gate, vary Rs.
    cand = {
        struct('label','KF-base'      , 'Qr',2.0,'Qd',0.25 ,'gamma',25,'Rs',0.5);
        struct('label','KF-Qd/2'      , 'Qr',2.0,'Qd',0.125,'gamma',25,'Rs',0.5);
        struct('label','KF-Qd/4'      , 'Qr',2.0,'Qd',0.06 ,'gamma',25,'Rs',0.5);
        struct('label','KF-gate16'    , 'Qr',2.0,'Qd',0.25 ,'gamma',16,'Rs',0.5);
        struct('label','KF-gate9'     , 'Qr',2.0,'Qd',0.25 ,'gamma', 9,'Rs',0.5);
        struct('label','KF-Rs2'       , 'Qr',2.0,'Qd',0.25 ,'gamma',25,'Rs',2.0);
        struct('label','KF-Rs2-gate16', 'Qr',2.0,'Qd',0.25 ,'gamma',16,'Rs',2.0);
        struct('label','KF-Qd/4-g16'  , 'Qr',2.0,'Qd',0.06 ,'gamma',16,'Rs',0.5);
        struct('label','KF-Qr/2'      , 'Qr',1.0,'Qd',0.25 ,'gamma',25,'Rs',0.5);
        struct('label','KF-Qr/2-Qd/2' , 'Qr',1.0,'Qd',0.125,'gamma',25,'Rs',0.5);
    };

    fprintf('%-18s | %-4s | %8s %8s %8s %8s %8s\n', ...
        'tuning','seed','LL_r','LL_d','RMSE_r','RMSE_d','NIS');
    fprintf('%s\n', repmat('-', 1, 82));

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
            rst = base.std_meas * c.Rs;

            flt = CSKF(dt, sa, rst(1), rst(2), init);
            flt.cs_alpha = base.alpha;
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

            fprintf('%-18s | %-4d | %8.3f %8.3f %8.3f %8.3f %8.3f\n', ...
                c.label, r, ...
                mean(ll_r(ss),'omitnan'), mean(ll_d(ss),'omitnan'), ...
                sqrt(mean(er(sse).^2,'omitnan')), sqrt(mean(ed(sse).^2,'omitnan')), ...
                mean(nis_a(ss),'omitnan'));
        end
        fprintf('\n');
    end
end
