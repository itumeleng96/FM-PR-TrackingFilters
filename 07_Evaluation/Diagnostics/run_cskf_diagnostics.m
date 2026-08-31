% =============================================================
% run_cskf_diagnostics.m
%
% CSKF on all 4 scenarios with full per-scan diagnostics.
% For each scenario, produces a 2x3 figure:
%   (1,1) Range: truth / measurement / estimate  (+ trigger marks)
%   (1,2) Doppler: truth / measurement / estimate (+ trigger marks)
%   (1,3) Normalised innovation eta_r, eta_d vs trigger threshold
%   (2,1) Running RMSE (range, doppler dual-axis)
%   (2,2) Per-scan log-likelihood (Gaussian, uses adaptive S)
%   (2,3) Per-scan log-likelihood using NOMINAL R (pre-adaptation) —
%         this is the "did the filter FLAG the outlier" signal
%
% Saves figures/cskf_<scenario>.png and prints a summary.
% =============================================================

clc; clear; close all;

addpath('07_Evaluation/GroundTruth/', '06_Tracking/Filters/KF/');
if ~exist('figures', 'dir'); mkdir('figures'); end

scenarios = {'levelFlight', 'landing', 'takeoff', 'orbit360'};

dt = 1;
kf_std_meas = [4.9038, 0.9985];
% Scheme B2 winning tuning for CSKF (Q_scale_r=2, Q_scale_d=1).
Q_scale_r   = 2.0;
Q_scale_d   = 1.0;
kf_std_acc  = [0.0048354 * Q_scale_r, 0.0991 * Q_scale_d];
% Scheme B2 CS parameters
cs_alpha_tuned = 0.90;
cs_beta_tuned  = 0.90;
d_r = 1.5; d_d = 1.0;                       % FM PR bin sizes (km, Hz)

summary = struct();

for s = 1:numel(scenarios)
    name = scenarios{s};
    cache = sprintf('meas_%s.mat', name);
    S = load(cache); d = S.data;
    N = d.N_align;

    fprintf('\n==== CSKF | %s (N=%d) ====\n', name, N);

    initial_state = [d.meas(1,1); 0; d.meas(2,1); 0];
    flt = CSKF(dt, kf_std_acc, kf_std_meas(1), kf_std_meas(2), initial_state);
    flt.cs_alpha   = cs_alpha_tuned;
    flt.cs_beta    = cs_beta_tuned;
    flt.cs_adapt_q = true;                  % adaptive Q on (per-axis)

    est_r     = zeros(1, N);
    est_d     = zeros(1, N);
    eta_r     = zeros(1, N);      % d^2 / S(1,1)
    eta_d     = zeros(1, N);      % d^2 / S(2,2)
    trig_r    = false(1, N);
    trig_d    = false(1, N);
    ll_adapt  = zeros(1, N);      % log-likelihood with adapted S
    ll_nomR   = zeros(1, N);      % log-likelihood with nominal R (outlier flag)
    rmse_r_c  = zeros(1, N);
    rmse_d_c  = zeros(1, N);

    R_nom = flt.R_nominal;

    for k = 1:N
        z = d.meas(:, k);
        [~,  flt] = flt.predict();

        % Snapshot pre-update quantities for LL/trigger diagnostics.
        S_pre  = flt.H * flt.P * flt.H.' + flt.R;
        S_nom  = flt.H * flt.P * flt.H.' + R_nom;
        d_pre  = z - flt.H * flt.X;
        eta_r(k) = d_pre(1)^2 / S_pre(1,1);
        eta_d(k) = d_pre(2)^2 / S_pre(2,2);

        % Log-likelihood diagnostics
        if all(isfinite(d_pre)) && det(S_pre) > 0
            ll_adapt(k) = -0.5*(2*log(2*pi) + log(det(S_pre)) + d_pre.' * (S_pre \ d_pre));
        end
        if all(isfinite(d_pre)) && det(S_nom) > 0
            ll_nomR(k)  = -0.5*(2*log(2*pi) + log(det(S_nom)) + d_pre.' * (S_nom \ d_pre));
        end

        % Detect trigger by watching residual-buffer trimming (CSKF trims on fire).
        n_before_r = numel(flt.residuals_x);
        n_before_d = numel(flt.residuals_y);
        [xk, flt] = flt.update(z);
        n_after_r = numel(flt.residuals_x);
        n_after_d = numel(flt.residuals_y);
        trig_r(k) = n_after_r < n_before_r + 1;   % buffer trimmed => fired
        trig_d(k) = n_after_d < n_before_d + 1;

        est_r(k) = xk(1);
        est_d(k) = xk(3);

        er = est_r(1:k)' - d.true_r(1:k);
        ed = est_d(1:k)' - d.true_d(1:k);
        rmse_r_c(k) = sqrt(mean(er.^2));
        rmse_d_c(k) = sqrt(mean(ed.^2));
    end

    rmse_r = rmse_r_c(end);
    rmse_d = rmse_d_c(end);
    score  = rmse_r/d_r + rmse_d/d_d;
    n_trig_r = sum(trig_r);
    n_trig_d = sum(trig_d);

    fprintf('  RMSE     range=%.2f km   doppler=%.2f Hz   score=%.2f (bins)\n', ...
            rmse_r, rmse_d, score);
    fprintf('  Triggers range=%d   doppler=%d (of %d scans)\n', ...
            n_trig_r, n_trig_d, N);
    fprintf('  LL(adapt) sum=%.1f   mean/scan=%.3f\n', sum(ll_adapt), mean(ll_adapt));
    fprintf('  LL(nomR)  sum=%.1f   mean/scan=%.3f (spikes flag outliers)\n', ...
            sum(ll_nomR), mean(ll_nomR));

    summary.(name) = struct('rmse_r',rmse_r,'rmse_d',rmse_d,'score',score, ...
                            'n_trig_r',n_trig_r,'n_trig_d',n_trig_d, ...
                            'll_adapt',sum(ll_adapt),'ll_nomR',sum(ll_nomR));

    % ---- Plot ----
    t = 1:N;
    trig_r_k = t(trig_r);
    trig_d_k = t(trig_d);

    fig = figure('Position',[100 100 1500 750],'Visible','off');

    % Range panel
    subplot(2,3,1); hold on; grid on;
    plot(t, d.meas(1,1:N),'k.','MarkerSize',6,'DisplayName','measurement');
    plot(t, d.true_r,     'g-','LineWidth',2,'DisplayName','truth');
    plot(t, est_r,        'b-','LineWidth',1.5,'DisplayName','estimate');
    if ~isempty(trig_r_k)
        yl = ylim;
        for i = 1:numel(trig_r_k), xline(trig_r_k(i),'r:','Alpha',0.4); end
        ylim(yl);
    end
    xlabel('scan'); ylabel('bistatic range (km)');
    title(sprintf('%s | range  RMSE=%.2f km  triggers=%d', name, rmse_r, n_trig_r));
    legend('Location','best');

    % Doppler panel
    subplot(2,3,2); hold on; grid on;
    plot(t, d.meas(2,1:N),'k.','MarkerSize',6);
    plot(t, d.true_d,     'g-','LineWidth',2);
    plot(t, est_d,        'b-','LineWidth',1.5);
    if ~isempty(trig_d_k)
        yl = ylim;
        for i = 1:numel(trig_d_k), xline(trig_d_k(i),'r:','Alpha',0.4); end
        ylim(yl);
    end
    xlabel('scan'); ylabel('bistatic doppler (Hz)');
    title(sprintf('doppler  RMSE=%.2f Hz  triggers=%d', rmse_d, n_trig_d));

    % Normalised innovation (eta) with trigger thresholds
    subplot(2,3,3); hold on; grid on;
    semilogy(t, max(eta_r,1e-3), 'r-', 'LineWidth',1.2, 'DisplayName','\eta_r');
    semilogy(t, max(eta_d,1e-3), 'b-', 'LineWidth',1.2, 'DisplayName','\eta_d');
    yline(1, 'k--', 'threshold \tau=1', 'Alpha',0.6);
    set(gca,'YScale','log');
    xlabel('scan'); ylabel('\eta = d^2 / S (log)');
    title('Normalised innovation (outlier statistic)');
    legend('Location','best');

    % Running RMSE
    subplot(2,3,4); grid on;
    yyaxis left;  plot(t, rmse_r_c,'b-','LineWidth',1.5); ylabel('RMSE range (km)');
    yyaxis right; plot(t, rmse_d_c,'r-','LineWidth',1.5); ylabel('RMSE doppler (Hz)');
    xlabel('scan'); title('Running RMSE');

    % Log-likelihood (adapted S)
    subplot(2,3,5); hold on; grid on;
    plot(t, ll_adapt, 'b-','LineWidth',1.2);
    xlabel('scan'); ylabel('log-likelihood (adapted R)');
    title(sprintf('LL with adapted R  sum=%.1f', sum(ll_adapt)));

    % Log-likelihood (nominal R) — spikes flag outliers
    subplot(2,3,6); hold on; grid on;
    plot(t, ll_nomR, 'm-','LineWidth',1.2);
    if ~isempty(trig_r_k)
        for i = 1:numel(trig_r_k), xline(trig_r_k(i),'r:','Alpha',0.4); end
    end
    if ~isempty(trig_d_k)
        for i = 1:numel(trig_d_k), xline(trig_d_k(i),'b:','Alpha',0.4); end
    end
    xlabel('scan'); ylabel('log-likelihood (nominal R)');
    title('LL with fixed nominal R (dips flag outliers)');

    sgtitle(sprintf('CSKF diagnostics | %s', name), 'FontWeight','bold','FontSize',13);

    out_png = sprintf('figures/cskf_%s.png', name);
    exportgraphics(fig, out_png, 'Resolution', 130);
    close(fig);
    fprintf('  Saved %s\n', out_png);
end

% ---- Summary table ----
fprintf('\n\n=============================================================\n');
fprintf('  CSKF summary across scenarios (adaptive Q on, per-axis)\n');
fprintf('=============================================================\n');
fprintf('%-14s | %8s %8s %8s | %6s %6s | %10s %10s\n', ...
        'Scenario','RMSE_r','RMSE_d','score','tr_r','tr_d','LL_adapt','LL_nomR');
for s = 1:numel(scenarios)
    n = scenarios{s}; m = summary.(n);
    fprintf('%-14s | %6.2fkm %6.1fHz %8.2f | %6d %6d | %10.1f %10.1f\n', ...
            n, m.rmse_r, m.rmse_d, m.score, m.n_trig_r, m.n_trig_d, m.ll_adapt, m.ll_nomR);
end
save('cskf_diagnostics.mat', 'summary');
fprintf('\nSaved figures/ and cskf_diagnostics.mat\n');
