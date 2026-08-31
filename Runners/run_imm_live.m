function out = run_imm_live(scenario_name)
% RUN_IMM_LIVE
%   Two-model IMM (CV low-Q + CV high-Q) on a scenario, live ARD update
%   per scan with model-probability trace. Same visual layout as
%   run_cskf_live.m.

    addpath('01_FERS/', '04_Detection/', '05_Clustering/', '03_DPI_Cancellation/', ...
            '07_Evaluation/GroundTruth/', '06_Tracking/Filters/IMM/');

    scen_secs = struct('levelFlight',60,'landing',120,'takeoff',60,'orbit360',120);
    if ~isfield(scen_secs, scenario_name)
        error('Unknown scenario "%s".', scenario_name);
    end
    N_scans   = scen_secs.(scenario_name);
    direct_h5 = sprintf('direct_%s.h5', scenario_name);
    echo_h5   = sprintf('echo_%s.h5',   scenario_name);

    fs = 200000; dopp_bins = 200; delay = 233e-6; c = 299792458;
    range_delay = delay*c;
    proc = struct('cancellationMaxRange_m',range_delay,'cancellationMaxDoppler_Hz',4, ...
                  'TxToRefRxDistance_m',12540,'nSegments',1,'nIterations',20, ...
                  'Fs',fs,'alpha',0,'initialAlpha',0);

    % ---- IMM tuning ----
    dt = 1;
    r_std_meas    = [4.9038, 0.9985];
    q_cruise      = [0.0048354, 0.0991];    % low-Q model (cruise)
    q_manoeuvre   = q_cruise .* 50;         % high-Q model (turns, ramps)

    [~, true_r, true_d] = computeGroundTruth(scenario_name);
    N_gt = numel(true_r);

    fprintf('Loading FERS baseband for %s...\n', scenario_name);
    [Ino,  Qno,  scale_no ] = loadfersHDF5(direct_h5);
    [Imov, Qmov, scale_mov] = loadfersHDF5(echo_h5);
    I_Qmov = (Imov + 1i*Qmov) .* scale_mov;
    I_Qno  = (Ino  + 1i*Qno ) .* scale_no;

    Ndelay = floor(delay*fs); Ndop = ceil(fs*dopp_bins/fs);
    range_axis_km = (0:Ndelay).' * (c/fs) / 1000;
    dopp_axis_hz  = (-Ndop:Ndop).';

    meas_r   = nan(1, N_scans); meas_d = nan(1, N_scans);
    est_r    = nan(1, N_scans); est_d  = nan(1, N_scans);
    mu_hist  = nan(2, N_scans);
    n_coast  = 0;   % scans where no cluster passed the gate
    gamma_assoc = 16;   % chi^2 2-dof, ~99.97% (loose, still rejects clutter)
    rmse_r_c = nan(1, N_scans); rmse_d_c = nan(1, N_scans);
    ll_scan  = nan(1, N_scans); nis_scan = nan(1, N_scans); nees_obs = nan(1, N_scans);
    loglik   = 0;  log2pi = log(2*pi);

    fig = figure('Name', sprintf('IMM live | %s', scenario_name), ...
                 'Position', [80 80 1100 700]);

    flt = [];
    initial = 1; current = fs;

    fprintf('\n==== IMM_CV2 | %s (N=%d) ====\n', scenario_name, N_scans);
    fprintf('%-4s | %-10s %-10s | %-10s %-10s | %-8s %-8s | %-12s\n', ...
            'scan','meas_r','meas_d','est_r','est_d','RMSE_r','RMSE_d','mu(low,high)');

    for k = 1:N_scans
        s1 = I_Qmov(initial:current);
        s2 = I_Qno (initial:current);
        s1 = procECA_Optimized(s2, s1, proc);
        [y, ~] = ardNoPlot(s1, s2, fs, dopp_bins, delay, k, []);

        [targetClusters, ~, ~] = ca_cfar(y.', 1e-7, fs, dopp_bins, delay, 20);
        [clusterCentroids, ~, ~, ~] = meanShift(targetClusters, 3, 5);

        if isempty(flt)
            % Bootstrap: no filter yet, take strongest centroid or truth
            if ~isempty(clusterCentroids)
                z_meas = clusterCentroids(:, 1);
            else
                z_meas = [true_r(min(1,N_gt)); true_d(min(1,N_gt))];
            end
            meas_r(k) = z_meas(1);  meas_d(k) = z_meas(2);
            init_state = [z_meas(1); 0; z_meas(2); 0];
            flt = IMM_CV2(dt, q_cruise, q_manoeuvre, ...
                          r_std_meas(1), r_std_meas(2), init_state);
            [x_pred, flt] = flt.predict();
            [xk, flt] = flt.update(z_meas);
            est_r(k) = xk(1); est_d(k) = xk(3);
            mu_hist(:, k) = flt.mu;
        else
            % Predict, then GNN gate cluster centroids against prediction
            [x_pred, flt] = flt.predict();
            Sk = flt.S;
            [z_gated, idx, ~] = gated_association(clusterCentroids, x_pred, ...
                                                  flt.H, Sk, gamma_assoc);
            if any(isnan(z_gated))
                % No cluster passed the gate — coast on prediction
                n_coast   = n_coast + 1;
                meas_r(k) = NaN;   meas_d(k) = NaN;
                est_r(k)  = x_pred(1);  est_d(k) = x_pred(3);
                mu_hist(:, k) = flt.mu;
            else
                meas_r(k) = z_gated(1);  meas_d(k) = z_gated(2);
                nu = z_gated - flt.H * x_pred;
                if det(Sk) > 0
                    ll_scan(k)  = -0.5 * (2*log2pi + log(det(Sk)) + nu.' * (Sk \ nu));
                    loglik      = loglik + ll_scan(k);
                    nis_scan(k) = nu.' * (Sk \ nu);
                end
                [xk, flt] = flt.update(z_gated);
                est_r(k)  = xk(1);  est_d(k) = xk(3);
                mu_hist(:, k) = flt.mu;
            end
        end

        if k <= N_gt
            er = est_r(1:k)' - true_r(1:k);
            ed = est_d(1:k)' - true_d(1:k);
            rmse_r_c(k) = sqrt(mean(er.^2, 'omitnan'));
            rmse_d_c(k) = sqrt(mean(ed.^2, 'omitnan'));
            HPH_post = flt.H * flt.P * flt.H.';
            ze = [true_r(k); true_d(k)] - flt.H * flt.X;
            if det(HPH_post) > 0
                nees_obs(k) = ze.' * (HPH_post \ ze);
            end
        end

        fprintf('%4d | %7.2fkm %7.2fHz | %7.2fkm %7.2fHz | %6.2fkm %5.1fHz | %.2f,%.2f\n', ...
                k, meas_r(k), meas_d(k), est_r(k), est_d(k), ...
                rmse_r_c(k), rmse_d_c(k), flt.mu(1), flt.mu(2));

        clf(fig);
        imagesc(range_axis_km, dopp_axis_hz, 10*log10(y.' + eps));
        axis xy; colormap(jet); hold on;
        xlabel('Bistatic range (km)','FontSize',12);
        ylabel('Bistatic doppler (Hz)','FontSize',12);
        title(sprintf('IMM_CV2 | %s | scan %d/%d  RMSE_r=%.2fkm  RMSE_d=%.1fHz  \\mu=[%.2f, %.2f]', ...
                      scenario_name, k, N_scans, rmse_r_c(k), rmse_d_c(k), flt.mu(1), flt.mu(2)), ...
              'FontSize', 12, 'Interpreter','tex');
        cb = colorbar; cb.Label.String = 'Level (dB)';

        if N_gt >= 1
            k_gt = min(k, N_gt);
            plot(true_r(1:k_gt), true_d(1:k_gt), 'g-', 'LineWidth', 1.5);
            plot(true_r(k_gt),   true_d(k_gt),   'g*', 'MarkerSize', 14, 'LineWidth', 2);
        end
        valid = ~isnan(meas_r);
        plot(meas_r(valid), meas_d(valid), 'k.', 'MarkerSize', 8);
        valid_est = ~isnan(est_r);
        plot(est_r(valid_est), est_d(valid_est), 'b-', 'LineWidth', 1.5);
        plot(est_r(k), est_d(k), 'b^', 'MarkerFaceColor','b','MarkerSize',8);
        if ~isempty(clusterCentroids)
            plot(clusterCentroids(1,1), clusterCentroids(2,1), 'wo', ...
                 'MarkerSize', 14, 'LineWidth', 2);
        end
        legend({'truth path','truth (*)','measurements','estimate path','estimate (\Delta)','strongest cluster (o)'}, ...
               'Location','northeast','TextColor','w','Color',[0 0 0 0.5],'FontSize',9);
        drawnow;

        initial = current + 1;
        current = current + fs;
    end

    N = min(N_scans, N_gt);
    er = est_r(1:N)' - true_r(1:N);
    ed = est_d(1:N)' - true_d(1:N);
    smooth_r = std(diff(est_r(1:N)), 'omitnan');
    smooth_d = std(diff(est_d(1:N)), 'omitnan');
    mean_nis = mean(nis_scan(1:N), 'omitnan');
    mean_nees = mean(nees_obs(1:N), 'omitnan');
    nis_ci_lo = 2 - 1.96*sqrt(4/N);  nis_ci_hi = 2 + 1.96*sqrt(4/N);

    out = struct( ...
        'scenario',    scenario_name, ...
        'N',           N, ...
        'rmse_r',      sqrt(mean(er.^2)), ...
        'rmse_d',      sqrt(mean(ed.^2)), ...
        'smooth_r',    smooth_r, ...
        'smooth_d',    smooth_d, ...
        'loglik',      loglik, ...
        'loglik_mean', loglik / N, ...
        'mean_nis',    mean_nis, ...
        'mean_nees',   mean_nees, ...
        'nis_ci',      [nis_ci_lo, nis_ci_hi], ...
        'mu_hist',     mu_hist(:,1:N), ...
        'est_r',       est_r(1:N), 'est_d', est_d(1:N), ...
        'll_scan',     ll_scan(1:N));

    fprintf('\n=============================================================\n');
    fprintf('  IMM_CV2 | %s | final stats\n', scenario_name);
    fprintf('=============================================================\n');
    fprintf('  N scans                  : %d\n', N);
    fprintf('  RMSE range               : %8.3f km\n',  out.rmse_r);
    fprintf('  RMSE doppler             : %8.3f Hz\n',  out.rmse_d);
    fprintf('  Smoothness range (std)   : %8.3f km\n',  out.smooth_r);
    fprintf('  Smoothness doppler (std) : %8.3f Hz\n',  out.smooth_d);
    fprintf('  Log-likelihood total     : %10.2f\n',   out.loglik);
    fprintf('  Log-likelihood per scan  : %10.4f\n',   out.loglik_mean);
    fprintf('  Mean NIS                 : %8.3f  (expected 2, 95%% CI [%.2f, %.2f])\n', ...
            out.mean_nis, out.nis_ci(1), out.nis_ci(2));
    fprintf('  Mean NEES (observable)   : %8.3f  (expected 2)\n', out.mean_nees);
    fprintf('  Final model probs        : mu_cruise=%.3f  mu_manoeuvre=%.3f\n', ...
            mu_hist(1, N), mu_hist(2, N));
    fprintf('  Coasted (no cluster)     : %d / %d scans (gate gamma=%.1f)\n', ...
            n_coast, N_scans, gamma_assoc);
    fprintf('=============================================================\n');

    % --- Summary figures ---
    t = 1:N;
    figure('Name','IMM | RMSE','Position',[90 90 800 420]);
    yyaxis left;  plot(t, rmse_r_c(1:N), 'b-', 'LineWidth', 1.8);
    grid on; ylabel('RMSE range (km)');
    yyaxis right; plot(t, rmse_d_c(1:N), 'r-', 'LineWidth', 1.8);
    ylabel('RMSE doppler (Hz)'); xlabel('scan');
    title(sprintf('IMM_CV2 | %s | Running RMSE', scenario_name), 'Interpreter','none');

    figure('Name','IMM | model probabilities','Position',[90 540 800 300]);
    plot(t, mu_hist(1,1:N), 'b-', 'LineWidth', 1.6); hold on;
    plot(t, mu_hist(2,1:N), 'r-', 'LineWidth', 1.6); grid on;
    xlabel('scan'); ylabel('model probability'); ylim([0 1]);
    legend('cruise (low Q)','manoeuvre (high Q)','Location','best');
    title(sprintf('IMM_CV2 | %s | Model probabilities', scenario_name), 'Interpreter','none');
end
