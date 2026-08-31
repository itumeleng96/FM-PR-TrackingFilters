function out = run_filter_live(filter_name, scenario_name)
% RUN_FILTER_LIVE
%   Live ARD run for a chosen filter and scenario, using the tuned
%   parameters from optimize_all_filters.m and GNN gating.
%
%   filter_name  : 'KF' | 'UKF' | 'RGNF' | 'PF'
%   scenario_name: 'levelFlight' | 'landing' | 'orbit360' | 'takeoff'
%
%   Example:
%     out = run_filter_live('UKF', 'landing');

    addpath('01_FERS/', '04_Detection/', '05_Clustering/', '03_DPI_Cancellation/', ...
            '07_Evaluation/GroundTruth/', ...
            '06_Tracking/Filters/KF/', '06_Tracking/Filters/UKF/', ...
            '06_Tracking/Filters/RGNF/', '06_Tracking/Filters/PF/');

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

    % ---- Tuned params from optimize_all_filters.m ----
    dt = 1;
    switch upper(filter_name)
        case 'KF'
            std_meas_base = [4.9038, 0.9985];   std_acc_base = [0.0048354, 0.0991];
            alpha = 0.98; gamma_assoc = 25; Rs = 0.5; Qr = 2.0; Qd = 0.25;
        case 'UKF'
            std_meas_base = [0.9707, 0.79739];  std_acc_base = [0.0076533, 0.09938];
            alpha = 0.98; gamma_assoc = 25; Rs = 2.0; Qr = 2.0; Qd = 0.25;
            ukf_a = 0.5; ukf_k = 0; ukf_b = 2;
        case 'RGNF'
            std_meas_base = [2.046, 0.98];      std_acc_base = [0.057027, 0.047789];
            alpha = 0.98; gamma_assoc = 25; Rs = 1.0; Qr = 0.5; Qd = 0.25;
            rgnf_iter = 100; rgnf_lambda = 1;
        case 'PF'
            std_meas_base = [10, 2];            std_acc_base = [1.0, 1.0];
            alpha = 0.90; gamma_assoc = 16; Rs = 2.0; Qr = 0.5; Qd = 5.0;
            pf_N = 1000;
        otherwise
            error('Unknown filter "%s". Use KF | UKF | RGNF | PF.', filter_name);
    end
    sa  = [std_acc_base(1)*Qr, std_acc_base(2)*Qd];
    rst = std_meas_base * Rs;

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
    rmse_r_c = nan(1, N_scans); rmse_d_c = nan(1, N_scans);
    nis_scan = nan(1, N_scans); nees_obs = nan(1, N_scans);
    ll_scan  = nan(1, N_scans); loglik = 0; log2pi = log(2*pi);
    n_coast  = 0;

    fig = figure('Name', sprintf('%s live | %s', filter_name, scenario_name), ...
                 'Position', [80 80 1100 700]);

    flt = [];
    initial = 1; current = fs;

    fprintf('\n==== %s (tuned) | %s (N=%d) ====\n', filter_name, scenario_name, N_scans);
    fprintf('%-4s | %-10s %-10s | %-10s %-10s | %-8s %-8s\n', ...
            'scan','meas_r','meas_d','est_r','est_d','RMSE_r','RMSE_d');

    for k = 1:N_scans
        s1 = I_Qmov(initial:current);
        s2 = I_Qno (initial:current);
        s1 = procECA_Optimized(s2, s1, proc);
        [y, ~] = ardNoPlot(s1, s2, fs, dopp_bins, delay, k, []);
        [targetClusters, ~, ~]      = ca_cfar(y.', 1e-7, fs, dopp_bins, delay, 20);
        [clusterCentroids, ~, ~, ~] = meanShift(targetClusters, 3, 5);

        if isempty(flt)
            % Bootstrap: strongest centroid
            if ~isempty(clusterCentroids)
                z_boot = clusterCentroids(:, 1);
            else
                z_boot = [true_r(min(1,N_gt)); true_d(min(1,N_gt))];
            end
            meas_r(k) = z_boot(1); meas_d(k) = z_boot(2);
            init = [z_boot(1); 0; z_boot(2); 0];
            switch upper(filter_name)
                case 'KF';   flt = CSKF(dt, sa, rst(1), rst(2), init);
                case 'UKF';  flt = CSUKF(dt, sa, rst(1), rst(2), init', ukf_a, ukf_k, ukf_b);
                case 'RGNF'; flt = CSRGNF(dt, sa, rst(1), rst(2), init, rgnf_iter, rgnf_lambda);
                case 'PF';   flt = CSPF(dt, sa, rst, init, pf_N);
            end
            flt.cs_alpha = alpha;
            if isprop(flt, 'cs_adapt_q'), flt.cs_adapt_q = true; end
            [~, flt]  = flt.predict();
            [xk, flt] = flt.update(z_boot);
            est_r(k) = xk(1); est_d(k) = xk(3);
        else
            [x_pred, flt] = flt.predict();
            Sk = flt.S;
            [z_gated, ~, ~] = gated_association(clusterCentroids, x_pred, flt.H, Sk, gamma_assoc);
            if any(isnan(z_gated))
                n_coast   = n_coast + 1;
                meas_r(k) = NaN; meas_d(k) = NaN;
                est_r(k)  = x_pred(1); est_d(k) = x_pred(3);
            else
                meas_r(k) = z_gated(1); meas_d(k) = z_gated(2);
                nu = z_gated - flt.H * x_pred;
                if det(Sk) > 0
                    ll_scan(k)  = -0.5*(2*log2pi + log(det(Sk)) + nu.'*(Sk\nu));
                    loglik      = loglik + ll_scan(k);
                    nis_scan(k) = nu.'*(Sk\nu);
                end
                [xk, flt] = flt.update(z_gated);
                est_r(k)  = xk(1); est_d(k) = xk(3);
            end
        end

        if k <= N_gt
            er = est_r(1:k)' - true_r(1:k);
            ed = est_d(1:k)' - true_d(1:k);
            rmse_r_c(k) = sqrt(mean(er.^2, 'omitnan'));
            rmse_d_c(k) = sqrt(mean(ed.^2, 'omitnan'));
            HPH_post = flt.H * flt.P * flt.H.';
            x_col    = flt.X(:);
            ze       = [true_r(k); true_d(k)] - flt.H * x_col;
            if det(HPH_post) > 0
                nees_obs(k) = ze.' * (HPH_post \ ze);
            end
        end

        fprintf('%4d | %7.2fkm %7.2fHz | %7.2fkm %7.2fHz | %6.2fkm %5.1fHz\n', ...
                k, meas_r(k), meas_d(k), est_r(k), est_d(k), rmse_r_c(k), rmse_d_c(k));

        clf(fig);
        imagesc(range_axis_km, dopp_axis_hz, 10*log10(y.' + eps));
        axis xy; colormap(jet); hold on;
        xlabel('Bistatic range (km)','FontSize',12);
        ylabel('Bistatic doppler (Hz)','FontSize',12);
        title(sprintf('%s (tuned) | %s | scan %d/%d  RMSE_r=%.2fkm  RMSE_d=%.1fHz', ...
                      filter_name, scenario_name, k, N_scans, rmse_r_c(k), rmse_d_c(k)), ...
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

    out = struct( ...
        'filter',      filter_name, ...
        'scenario',    scenario_name, ...
        'N',           N, ...
        'rmse_r',      sqrt(mean(er.^2)), ...
        'rmse_d',      sqrt(mean(ed.^2)), ...
        'mae_r',       mean(abs(er)), ...
        'mae_d',       mean(abs(ed)), ...
        'smooth_r',    std(diff(est_r(1:N)), 'omitnan'), ...
        'smooth_d',    std(diff(est_d(1:N)), 'omitnan'), ...
        'n_coast',     n_coast, ...
        'loglik',      loglik, ...
        'loglik_mean', loglik / N, ...
        'mean_nis',    mean(nis_scan(1:N), 'omitnan'), ...
        'mean_nees',   mean(nees_obs(1:N), 'omitnan'));

    fprintf('\n=============================================================\n');
    fprintf('  %s (tuned) | %s | final stats\n', filter_name, scenario_name);
    fprintf('=============================================================\n');
    fprintf('  RMSE range               : %8.3f km\n', out.rmse_r);
    fprintf('  RMSE doppler             : %8.3f Hz\n', out.rmse_d);
    fprintf('  MAE  range               : %8.3f km\n', out.mae_r);
    fprintf('  MAE  doppler             : %8.3f Hz\n', out.mae_d);
    fprintf('  Smoothness range (std)   : %8.3f km\n', out.smooth_r);
    fprintf('  Smoothness doppler (std) : %8.3f Hz\n', out.smooth_d);
    fprintf('  Coasted scans            : %d / %d\n', out.n_coast, N_scans);
    fprintf('  Log-likelihood per scan  : %10.4f\n', out.loglik_mean);
    fprintf('  Mean NIS                 : %8.3f  (expected 2)\n', out.mean_nis);
    fprintf('  Mean NEES (observable)   : %8.3f  (expected 2)\n', out.mean_nees);
    fprintf('=============================================================\n');

    % ---- End-of-run summary figures ----
    t = 1:N;

    % Figure: RMSE over time (dual-axis)
    figure('Name', sprintf('%s | %s | RMSE', filter_name, scenario_name), ...
           'Position', [90 90 800 400]);
    yyaxis left;
    plot(t, rmse_r_c(1:N), 'b-', 'LineWidth', 1.8); grid on;
    ylabel('RMSE range (km)');
    yyaxis right;
    plot(t, rmse_d_c(1:N), 'r-', 'LineWidth', 1.8);
    ylabel('RMSE doppler (Hz)');
    xlabel('scan');
    title(sprintf('%s | %s | Running RMSE   final: %.2fkm / %.1fHz', ...
                  filter_name, scenario_name, out.rmse_r, out.rmse_d));

    % Figure: Log-likelihood over time
    figure('Name', sprintf('%s | %s | Log-likelihood', filter_name, scenario_name), ...
           'Position', [90 500 800 400]);
    plot(t, ll_scan(1:N), 'm-', 'LineWidth', 1.4); grid on;
    xlabel('scan'); ylabel('log-likelihood (per scan)');
    title(sprintf('%s | %s | Per-scan LL   sum=%.1f   mean=%.3f', ...
                  filter_name, scenario_name, out.loglik, out.loglik_mean));

    % Figure: NIS and NEES over time (dual-axis) with chi-squared band
    figure('Name', sprintf('%s | %s | NIS / NEES consistency', filter_name, scenario_name), ...
           'Position', [900 90 800 400]);
    yyaxis left;
    plot(t, nis_scan(1:N), 'b.-'); hold on;
    yline(2, 'b--', 'LineWidth', 1.2);
    ylabel('NIS (expected 2)');
    yyaxis right;
    plot(t, nees_obs(1:N), 'r.-');
    yline(2, 'r--', 'LineWidth', 1.2);
    ylabel('NEES obs (expected 2)');
    xlabel('scan'); grid on;
    title(sprintf('%s | %s | Consistency  mean NIS=%.2f  mean NEES=%.2f', ...
                  filter_name, scenario_name, out.mean_nis, out.mean_nees));

    out.ll_scan  = ll_scan(1:N);
    out.rmse_r_c = rmse_r_c(1:N);
    out.rmse_d_c = rmse_d_c(1:N);
    out.nis_scan = nis_scan(1:N);
    out.nees_obs = nees_obs(1:N);
end
