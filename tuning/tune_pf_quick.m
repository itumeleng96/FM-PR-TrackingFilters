function tune_pf_quick()
% TUNE_PF_QUICK  Compare 3 PF tunings (baseline + 2 candidates) on a small
% number of seeds per scenario. No MC — this is a diagnostic sweep to pick
% an explainable tuning before committing to the full 100-seed run.

    addpath('TrackingFilter-CS-ParticleFilter/','groundTruthCalculations/');

    seeds     = [1, 25, 50];              % 3 seeds spread across the batch
    scenarios = {'landing','orbit360','takeoff'};
    dt        = 1;
    warmup    = 10;
    log2pi    = log(2*pi);

    % Candidate tunings (Rs, Qr, Qd, alpha, gamma, N particles)
    tunings = struct( ...
        'Baseline', struct('Rs',2.0,'Qr',0.5,'Qd',5.0,'alpha',0.90,'gamma',16,'N',10000), ...
        'CandA',    struct('Rs',0.5,'Qr',0.5,'Qd',0.25,'alpha',0.98,'gamma',25,'N',10000), ...
        'CandB',    struct('Rs',1.0,'Qr',0.5,'Qd',0.5,'alpha',0.98,'gamma',25,'N',10000));
    tune_names = fieldnames(tunings);

    fprintf('%-9s | %-9s | %-8s | %8s | %8s | %8s | %8s\n', ...
        'scenario', 'seed', 'tuning', 'LL_r', 'LL_d', 'RMSE_r', 'NIS');
    fprintf('%s\n', repmat('-', 1, 72));

    base_meas = [10, 2];    % PF nominal std_meas
    base_acc  = [1.0, 1.0]; % PF nominal std_acc

    for si = 1:numel(scenarios)
        scen = scenarios{si};
        for r = seeds
            cfile = fullfile('seeds', sprintf('seed_%03d', r), ...
                             sprintf('clusters_%s.mat', scen));
            if ~isfile(cfile), continue; end
            d = load(cfile);
            N = d.N_scans; Nc = min(N, numel(d.true_r));

            if isempty(d.scans{1}), z0 = [d.true_r(1); d.true_d(1)];
            else, z0 = d.scans{1}(:,1); end
            init = [z0(1); 0; z0(2); 0];

            for ti = 1:numel(tune_names)
                tn = tune_names{ti};
                tp = tunings.(tn);
                sa = [base_acc(1)*tp.Qr, base_acc(2)*tp.Qd];
                rst = base_meas * tp.Rs;

                try
                    flt = CSPF(dt, sa, rst, init, tp.N);
                    flt.cs_alpha = tp.alpha;
                    if isprop(flt,'cs_adapt_q'), flt.cs_adapt_q = true; end

                    est_r = nan(1,N); est_d = nan(1,N);
                    ll_r = nan(1,N); ll_d = nan(1,N); nis_a = nan(1,N);

                    [~, flt]  = flt.predict();
                    [xk, flt] = flt.update(z0);
                    est_r(1) = xk(1); est_d(1) = xk(3);

                    for k = 2:N
                        [x_pred, flt] = flt.predict();
                        Sk = flt.S;
                        [z_g, ~, ~] = gated_association(d.scans{k}, x_pred, ...
                                                        flt.H, Sk, tp.gamma);
                        if any(isnan(z_g))
                            est_r(k) = x_pred(1); est_d(k) = x_pred(3);
                        else
                            nu = z_g - flt.H * x_pred;
                            if det(Sk) > 0
                                sr = max(Sk(1,1), 1e-9);
                                sd = max(Sk(2,2), 1e-9);
                                ll_r(k) = -0.5*(log2pi + log(sr) + nu(1)^2/sr);
                                ll_d(k) = -0.5*(log2pi + log(sd) + nu(2)^2/sd);
                                nis_a(k) = nu.' * (Sk \ nu);
                            end
                            [xk, flt] = flt.update(z_g);
                            est_r(k) = xk(1); est_d(k) = xk(3);
                        end
                    end

                    er = est_r(1:Nc)' - d.true_r(1:Nc);
                    ss = (warmup+1):numel(ll_r);
                    ss_e = (warmup+1):numel(er);

                    LLr = mean(ll_r(ss), 'omitnan');
                    LLd = mean(ll_d(ss), 'omitnan');
                    RMSEr = sqrt(mean(er(ss_e).^2, 'omitnan'));
                    NIS = mean(nis_a(ss), 'omitnan');

                    fprintf('%-9s | %-9d | %-8s | %8.3f | %8.3f | %8.3f | %8.3f\n', ...
                        scen, r, tn, LLr, LLd, RMSEr, NIS);
                catch ME
                    fprintf('%-9s | %-9d | %-8s | FAIL: %s\n', scen, r, tn, ME.message);
                end
            end
            fprintf('\n');
        end
    end
end
