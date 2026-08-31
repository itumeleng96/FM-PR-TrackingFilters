% =============================================================
% single_run_cspf.m
%
% Purpose:
%   One-shot single-target run of the refactored CSPF (no Monte
%   Carlo, no timing). Uses the same measurement pipeline as
%   runComputationalLoadRawFilter.m (FERS -> ECA -> ARD -> CFAR ->
%   mean-shift centroids) then feeds each per-scan centroid into
%   CSPF and plots estimate vs. measurement.
% =============================================================

clc; clear; close all;

addpath('01_FERS/', ...
        '04_Detection/', ...
        '05_Clustering/', ...
        '03_DPI_Cancellation/', ...
        '06_Tracking/Filters/PF/');

% -------- FERS data (skip regen if h5 files exist) --------
fers_bin = '/home/itumeleng/Documents/Academia/MscEng/FERS/build/src/fers';
sys_libs = '/usr/lib/gcc/x86_64-pc-linux-gnu/15:/usr/lib64';
if ~exist('direct.h5', 'file') || ~exist('echo.h5', 'file')
    fers_cmd = ['env LD_LIBRARY_PATH=' sys_libs ':$LD_LIBRARY_PATH ' ...
                fers_bin ' FERS/BackupScenarios/scenario_1_singleFile.fersxml'];
    system(fers_cmd);
end
[Ino,  Qno,  scale_no ] = loadfersHDF5('direct.h5');
[Imov, Qmov, scale_mov] = loadfersHDF5('echo.h5');

I_Qmov = (Imov + 1i*Qmov) .* scale_mov;
I_Qno  = (Ino  + 1i*Qno ) .* scale_no;

fs              = 200000;
simulation_time = 60;

dopp_bins   = 200;
delay       = 233e-6;
c           = 299792458;
range_delay = delay * c;

proc = struct('cancellationMaxRange_m',   range_delay, ...
              'cancellationMaxDoppler_Hz', 4, ...
              'TxToRefRxDistance_m',       12540, ...
              'nSegments',                 1, ...
              'nIterations',               20, ...
              'Fs',                        fs, ...
              'alpha',                     0, ...
              'initialAlpha',              0);

% -------- Per-scan measurement pipeline --------
measurements = zeros(2, simulation_time);
initial = 1; current = fs;
fprintf('Pre-computing measurement pipeline for %d scans...\n', simulation_time);
for i = 1:simulation_time
    s1 = I_Qmov(initial:current);
    s2 = I_Qno(initial:current);
    s1 = procECA_Optimized(s2, s1, proc);
    [y, ~] = ardNoPlot(s1, s2, fs, dopp_bins, delay, i, []);
    [targetClusters, ~, ~] = ca_cfar(y.', 1e-7, fs, dopp_bins, delay, 20);
    [clusterCentroids, ~, ~, ~] = meanShift(targetClusters, 10, 8);
    if ~isempty(clusterCentroids)
        measurements(:, i) = clusterCentroids(1:2, 1);
    elseif i > 1
        measurements(:, i) = measurements(:, i-1);
    end
    initial = current + 1;
    current = current + fs;
end

initial_state = [measurements(1,1); 0; measurements(2,1); 0];
fprintf('Init from scan 1 centroid: range=%.2f  doppler=%.2f\n', ...
        initial_state(1), initial_state(3));

% -------- CSPF setup (mirrors track.m case 3) --------
dt          = 1;
pf_std_acc  = [1.429, 1.9452];
pf_std_meas = [10, 2];
N_particles = 10000;

pf_obj = CSPF(dt, pf_std_acc, pf_std_meas, initial_state, N_particles);

% -------- Filter loop with diagnostics --------
est          = zeros(4, simulation_time);
neff_trace   = zeros(1, simulation_time);
eps_r_trace  = zeros(1, simulation_time);
eps_d_trace  = zeros(1, simulation_time);
outlier_r    = false(1, simulation_time);
outlier_d    = false(1, simulation_time);
fallback     = false(1, simulation_time);

M_window = 6;   % must match CSPF.m
for k = 1:simulation_time
    if k > 1 && all(measurements(:, k) == measurements(:, k-1))
        fallback(k) = true;
    end

    [~,  pf_obj] = pf_obj.predict();
    [xk, pf_obj] = pf_obj.update(measurements(:, k));
    est(:, k) = xk(:);

    neff_trace(k)  = 1 / sum(pf_obj.weights .^ 2);
    eps_r_trace(k) = pf_obj.epsRange(end);
    eps_d_trace(k) = pf_obj.epsDoppler(end);

    % Recompute the outlier trigger externally (matches CSPF.m logic).
    if numel(pf_obj.epsRange) > M_window
        w = pf_obj.epsRange(end-M_window : end-1);
        outlier_r(k) = eps_r_trace(k) > 2 && eps_r_trace(k) > (mean(w) + std(w));
    end
    if numel(pf_obj.epsDoppler) > M_window
        w = pf_obj.epsDoppler(end-M_window : end-1);
        outlier_d(k) = eps_d_trace(k) > 1 && eps_d_trace(k) > (mean(w) + std(w));
    end
end

% -------- Plot with diagnostics --------
t = 1:simulation_time;
figure('Name', 'CSPF single-run track', 'Position', [100 100 900 800]);

subplot(3,1,1);
plot(t, measurements(1,:), 'kx', 'DisplayName', 'Measurement');
hold on; grid on;
plot(t, est(1,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'CSPF est');
plot(t(fallback),  measurements(1, fallback),  'ms', 'MarkerSize', 10, ...
     'MarkerFaceColor', 'm', 'DisplayName', 'Missed-det fallback');
plot(t(outlier_r), est(1, outlier_r), 'ro', 'MarkerSize', 10, ...
     'LineWidth', 1.5, 'DisplayName', 'Range outlier trigger');
xlabel('scan'); ylabel('Range (m)'); legend('Location','best');
title('Range');

subplot(3,1,2);
plot(t, measurements(2,:), 'kx', 'DisplayName', 'Measurement');
hold on; grid on;
plot(t, est(3,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'CSPF est');
plot(t(fallback),  measurements(2, fallback),  'ms', 'MarkerSize', 10, ...
     'MarkerFaceColor', 'm', 'DisplayName', 'Missed-det fallback');
plot(t(outlier_d), est(3, outlier_d), 'ro', 'MarkerSize', 10, ...
     'LineWidth', 1.5, 'DisplayName', 'Doppler outlier trigger');
xlabel('scan'); ylabel('Doppler (Hz)'); legend('Location','best');
title('Doppler');

subplot(3,1,3);
semilogy(t, neff_trace, 'k-', 'LineWidth', 1.2); grid on;
yline(N_particles/2, 'r--', 'DisplayName', 'Resample threshold N/2');
xlabel('scan'); ylabel('N_{eff}');
title('Effective sample size');
legend('Location','best');

sgtitle('CSPF single-target run (post-refactor + diagnostics)');

% -------- Summary --------
res_range = est(1,:) - measurements(1,:);
res_dopp  = est(3,:) - measurements(2,:);
fprintf('\nEstimate-vs-measurement residual stats:\n');
fprintf('  Range   : mean=%+8.3f  std=%8.3f\n', mean(res_range), std(res_range));
fprintf('  Doppler : mean=%+8.3f  std=%8.3f\n', mean(res_dopp),  std(res_dopp));

fprintf('\nDiagnostic summary:\n');
fprintf('  Missed-detection fallbacks : %d / %d scans\n', sum(fallback), simulation_time);
fprintf('  Range outlier triggers     : %d\n', sum(outlier_r));
fprintf('  Doppler outlier triggers   : %d\n', sum(outlier_d));
fprintf('  NEFF min / median / max    : %.0f / %.0f / %.0f  (threshold %.0f)\n', ...
        min(neff_trace), median(neff_trace), max(neff_trace), N_particles/2);

% -------- Per-scan dump (for offline analysis) --------
fprintf('\n%3s %10s %10s %10s %10s %10s %10s %4s %4s %4s\n', ...
        'k','meas_r','meas_d','est_r','est_d','eps_r','eps_d','out_r','out_d','fb');
for k = 1:simulation_time
    fprintf('%3d %10.2f %10.2f %10.2f %10.2f %10.3f %10.3f %4d %4d %4d\n', ...
            k, measurements(1,k), measurements(2,k), est(1,k), est(3,k), ...
            eps_r_trace(k), eps_d_trace(k), ...
            outlier_r(k), outlier_d(k), fallback(k));
end
