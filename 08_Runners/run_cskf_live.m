function out = run_cskf_live(scenario_name)
% RUN_CSKF_LIVE
%   Run CSKF (with adaptive R + per-axis adaptive Q, Scheme B2 tuning) on
%   a scenario, live-updating the ARD figure per scan with:
%     - ground truth path (green)
%     - measurement history (black dots)
%     - filter estimate path (blue triangles + line)
%     - CFAR mean-shift centroids (white o, current scan only)
%
%   Range on X-axis (km), doppler on Y-axis (Hz).
%   Final summary printed at end.
%
%   out = run_cskf_live('levelFlight')
%   out = run_cskf_live('landing')
%   out = run_cskf_live('orbit360')

    addpath('01_FERS/', '04_Detection/', '05_Clustering/', '03_DPI_Cancellation/', ...
            '07_Evaluation/GroundTruth/', '06_Tracking/Filters/KF/');

    scen_secs = struct('levelFlight',60,'landing',120,'takeoff',60,'orbit360',120);
    if ~isfield(scen_secs, scenario_name)
        error('Unknown scenario "%s". Use levelFlight | landing | takeoff | orbit360.', scenario_name);
    end
    N_scans = scen_secs.(scenario_name);
    direct_h5 = sprintf('direct_%s.h5', scenario_name);
    echo_h5   = sprintf('echo_%s.h5',   scenario_name);

    fs = 200000; dopp_bins = 200; delay = 233e-6; c = 299792458;
    range_delay = delay*c;
    proc = struct('cancellationMaxRange_m',range_delay,'cancellationMaxDoppler_Hz',4, ...
                  'TxToRefRxDistance_m',12540,'nSegments',1,'nIterations',20, ...
                  'Fs',fs,'alpha',0,'initialAlpha',0);

    % ---- CSKF Scheme B2 tuning ----
    dt = 1;
    kf_std_meas = [4.9038, 0.9985];
    % ---- Tuned from optimize_cskf_akhlaghi.m sweep ----
    Q_scale_r = 2.0;    % tuned
    Q_scale_d = 0.5;    % tuned
    kf_std_acc = [0.0048354 * Q_scale_r, 0.0991 * Q_scale_d];
    cs_alpha   = 0.98;  % Akhlaghi forgetting factor
    cs_beta    = 0.98;  % kept equal for verbatim Akhlaghi
    cs_q_amax  = 5.0;   % legacy field, unused

    % ---- Ground truth ----
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

    % ---- Preallocate ----
    meas_r   = nan(1, N_scans);
    meas_d   = nan(1, N_scans);
    est_r    = nan(1, N_scans);
    est_d    = nan(1, N_scans);
    trig_r   = false(1, N_scans);
    trig_d   = false(1, N_scans);
    rmse_r_c = nan(1, N_scans);
    rmse_d_c = nan(1, N_scans);
    ll_scan  = nan(1, N_scans);
    nis_scan = nan(1, N_scans);   % d' S^-1 d, expected 2 (dim z) under consistency
    nees_obs = nan(1, N_scans);   % (z_true - Hx)' (HPH')^-1 (z_true - Hx), expected 2
    loglik   = 0;
    log2pi   = log(2*pi);

    fig = figure('Name', sprintf('CSKF live | %s', scenario_name), ...
                 'Position', [80 80 1100 700]);

    flt = [];
    initial = 1; current = fs;

    fprintf('\n==== CSKF | %s (N=%d) ====\n', scenario_name, N_scans);
    fprintf('%-4s | %-10s %-10s | %-10s %-10s | %-8s %-8s | %-4s %-4s\n', ...
            'scan','meas_r','meas_d','est_r','est_d','RMSE_r','RMSE_d','trR','trD');

    for k = 1:N_scans
        s1 = I_Qmov(initial:current);
        s2 = I_Qno (initial:current);
        s1 = procECA_Optimized(s2, s1, proc);
        [y, ~] = ardNoPlot(s1, s2, fs, dopp_bins, delay, k, []);

        [targetClusters, ~, ~] = ca_cfar(y.', 1e-7, fs, dopp_bins, delay, 20);
        [clusterCentroids, ~, ~, ~] = meanShift(targetClusters, 3, 5);

        % ---- Bootstrap filter on first scan (strongest centroid) ----
        if isempty(flt)
            if ~isempty(clusterCentroids)
                z_boot = clusterCentroids(:, 1);
            else
                z_boot = [true_r(min(1,N_gt)); true_d(min(1,N_gt))];
            end
            meas_r(k)  = z_boot(1); meas_d(k) = z_boot(2);
            init_state = [z_boot(1); 0; z_boot(2); 0];
            flt = CSKF(dt, kf_std_acc, kf_std_meas(1), kf_std_meas(2), init_state);
            flt.cs_alpha       = cs_alpha;
            flt.cs_beta        = cs_beta;
            flt.cs_adapt_q     = true;
            flt.cs_q_alpha_max = cs_q_amax;
            [~, flt] = flt.predict();
            [xk, flt] = flt.update(z_boot);
            est_r(k) = xk(1); est_d(k) = xk(3);
            initial = current + 1; current = current + fs;
            continue;
        end

        % ---- Predict, GNN-gate cluster centroids against prediction, update ----
        [x_pred, flt] = flt.predict();
        Sk = flt.S;
        gamma_assoc = 25;   % chi^2 2-dof ~5-sigma (tuned)
        [z_gated, ~, ~] = gated_association(clusterCentroids, x_pred, ...
                                            flt.H, Sk, gamma_assoc);
        if any(isnan(z_gated))
            % No cluster passed the gate — coast on prediction
            meas_r(k) = NaN;  meas_d(k) = NaN;
            est_r(k)  = x_pred(1);  est_d(k) = x_pred(3);
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
        end

        % ---- Running RMSE + NEES (observable-space, post-update) ----
        if k <= N_gt
            er = est_r(1:k)' - true_r(1:k);
            ed = est_d(1:k)' - true_d(1:k);
            rmse_r_c(k) = sqrt(mean(er.^2, 'omitnan'));
            rmse_d_c(k) = sqrt(mean(ed.^2, 'omitnan'));

            HPH_post = flt.H * flt.P * flt.H.';
            ze       = [true_r(k); true_d(k)] - flt.H * flt.X;
            if det(HPH_post) > 0
                nees_obs(k) = ze.' * (HPH_post \ ze);
            end
        end

        fprintf('%4d | %7.2fkm %7.2fHz | %7.2fkm %7.2fHz | %6.2fkm %5.1fHz | %-4s %-4s\n', ...
                k, meas_r(k), meas_d(k), est_r(k), est_d(k), ...
                rmse_r_c(k), rmse_d_c(k), ...
                ternary(trig_r(k),'*',' '), ternary(trig_d(k),'*',' '));

        % ---- Live plot ----
        clf(fig);
        imagesc(range_axis_km, dopp_axis_hz, 10*log10(y.' + eps));
        axis xy; colormap(jet); hold on;
        xlabel('Bistatic range (km)','FontSize',12);
        ylabel('Bistatic doppler (Hz)','FontSize',12);
        title(sprintf('CSKF | %s | scan %d/%d  RMSE_r=%.2fkm  RMSE_d=%.1fHz', ...
                      scenario_name, k, N_scans, rmse_r_c(k), rmse_d_c(k)), ...
              'FontSize', 13);
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

    % ---- Final summary ----
    N = min(N_scans, N_gt);
    er = est_r(1:N)' - true_r(1:N);
    ed = est_d(1:N)' - true_d(1:N);
    smooth_r = std(diff(est_r(1:N)), 'omitnan');
    smooth_d = std(diff(est_d(1:N)), 'omitnan');

    mean_nis  = mean(nis_scan(1:N), 'omitnan');
    mean_nees = mean(nees_obs(1:N), 'omitnan');
    % Chi-squared 95% CI on mean for dim=2, N samples:
    % [2 - 1.96*sqrt(4/N), 2 + 1.96*sqrt(4/N)]
    nis_ci_lo = 2 - 1.96*sqrt(4/N);  nis_ci_hi = 2 + 1.96*sqrt(4/N);

    out = struct( ...
        'scenario',    scenario_name, ...
        'N',           N, ...
        'rmse_r',      sqrt(mean(er.^2)), ...
        'rmse_d',      sqrt(mean(ed.^2)), ...
        'mae_r',       mean(abs(er)), ...
        'mae_d',       mean(abs(ed)), ...
        'smooth_r',    smooth_r, ...
        'smooth_d',    smooth_d, ...
        'n_trig_r',    sum(trig_r), ...
        'n_trig_d',    sum(trig_d), ...
        'loglik',      loglik, ...
        'loglik_mean', loglik / N, ...
        'mean_nis',    mean_nis, ...
        'mean_nees',   mean_nees, ...
        'nis_ci',      [nis_ci_lo, nis_ci_hi]);

    fprintf('\n=============================================================\n');
    fprintf('  CSKF | %s | final stats\n', scenario_name);
    fprintf('=============================================================\n');
    fprintf('  N scans                  : %d\n', N);
    fprintf('  RMSE range               : %8.3f km\n',  out.rmse_r);
    fprintf('  RMSE doppler             : %8.3f Hz\n',  out.rmse_d);
    fprintf('  MAE  range               : %8.3f km\n',  out.mae_r);
    fprintf('  MAE  doppler             : %8.3f Hz\n',  out.mae_d);
    fprintf('  Smoothness range (std)   : %8.3f km\n',  out.smooth_r);
    fprintf('  Smoothness doppler (std) : %8.3f Hz\n',  out.smooth_d);
    fprintf('  Triggers  range/doppler  : %d / %d  of %d scans\n', ...
            out.n_trig_r, out.n_trig_d, N);
    fprintf('  Log-likelihood total     : %10.2f\n',   out.loglik);
    fprintf('  Log-likelihood per scan  : %10.4f\n',   out.loglik_mean);
    fprintf('  Mean NIS                 : %8.3f  (expected 2, 95%% CI [%.2f, %.2f])\n', ...
            out.mean_nis, out.nis_ci(1), out.nis_ci(2));
    fprintf('  Mean NEES (observable)   : %8.3f  (expected 2, same CI)\n', out.mean_nees);
    fprintf('=============================================================\n');

    % ---- End-of-run summary figures ----
    t = 1:N;

    % Figure: RMSE over time (dual y-axis)
    figure('Name', sprintf('CSKF | %s | RMSE', scenario_name), ...
           'Position', [90 90 800 420]);
    yyaxis left;
    plot(t, rmse_r_c(1:N), 'b-', 'LineWidth', 1.8); grid on;
    ylabel('RMSE range (km)'); ylim([0 max(rmse_r_c(1:N))*1.1 + eps]);
    yyaxis right;
    plot(t, rmse_d_c(1:N), 'r-', 'LineWidth', 1.8);
    ylabel('RMSE doppler (Hz)'); ylim([0 max(rmse_d_c(1:N))*1.1 + eps]);
    xlabel('scan');
    title(sprintf('CSKF | %s | Running RMSE   final: %.2f km / %.1f Hz', ...
                  scenario_name, out.rmse_r, out.rmse_d), 'FontSize', 13);

    % Figure: Log-likelihood over time
    figure('Name', sprintf('CSKF | %s | Log-likelihood', scenario_name), ...
           'Position', [90 540 800 420]);
    plot(t, ll_scan(1:N), 'm-', 'LineWidth', 1.4); grid on; hold on;
    % Mark trigger scans (dips flag outliers)
    trig_k = t(trig_r(1:N) | trig_d(1:N));
    if ~isempty(trig_k)
        yl = ylim;
        for i = 1:numel(trig_k), xline(trig_k(i), 'k:', 'Alpha', 0.35); end
        ylim(yl);
    end
    xlabel('scan'); ylabel('log-likelihood (per scan)');
    title(sprintf('CSKF | %s | Per-scan LL   sum=%.1f   mean=%.3f', ...
                  scenario_name, out.loglik, out.loglik_mean), 'FontSize', 13);

    % The ARD figure (fig) already shows the FULL truth path, measurement
    % history, and estimate history at the last scan — leave it visible.
    figure(fig);

    out.ll_scan  = ll_scan;
    out.rmse_r_c = rmse_r_c;
    out.rmse_d_c = rmse_d_c;
    out.trig_r   = trig_r;
    out.trig_d   = trig_d;
end

function s = ternary(cond, t, f)
    if cond, s = t; else, s = f; end
end
