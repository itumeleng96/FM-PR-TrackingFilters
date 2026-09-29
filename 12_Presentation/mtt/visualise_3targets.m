% visualise_3targets.m
% -----------------------------------------------------------------
% Live-updating phosphor diagnostic for the 3-target demo:
%   Left : ARD (range-Doppler) of current frame
%   Right: CFAR detections + centroids (bright) + phosphor history
%
% Centroids from past frames linger as a phosphor trail —
% older ones fade with exponential decay.
%
% Run from repo root:
%   >> visualise_3targets
%
% Author: Itumeleng Malemela
% -----------------------------------------------------------------
clc; clear; close all;

addpath('01_FERS/', '04_Detection/', '05_Clustering/', ...
        '03_DPI_Cancellation/', '07_Evaluation/GroundTruth/', ...
        '06_Tracking/MTT/', ...
        '06_Tracking/Filters/KF/', '06_Tracking/Filters/UKF/', ...
        '06_Tracking/Filters/RGNF/', '06_Tracking/Filters/PF/');

% ---- Parameters ----
fs        = 200000;
dopp_bins = 200;
delay     = 233e-6;
c         = 299792458;
fc        = 94e6;
lambda    = c / fc;
Baseline_m = 12540;
Tx_Pos    = [6440; 10760; 1056];
RefRxPos  = [0; 0; 1000];
range_delay = delay * c;

proc = struct('cancellationMaxRange_m',   range_delay, ...
              'cancellationMaxDoppler_Hz', 4, ...
              'TxToRefRxDistance_m',       12540, ...
              'nSegments', 1, 'nIterations', 20, ...
              'Fs', fs, 'alpha', 0, 'initialAlpha', 0);

% ---- Ground truth ----
multi_tgts = struct( ...
    'T_slow', struct('pos', [15000, 12000, 2800; 15000, 18000, 2800], 'times', [0, 60]), ...
    'T_mid',  struct('pos', [20000, 14000, 3200; 20000, 23000, 3200], 'times', [0, 60]), ...
    'T_fast', struct('pos', [25000, 12000, 3600; 25000, 24000, 3600], 'times', [0, 60]));
t = (0:1:60)';
names = fieldnames(multi_tgts);
NT = numel(names);
gt_r = zeros(numel(t), NT);
gt_d = zeros(numel(t), NT);
for ti = 1:NT
    mt = multi_tgts.(names{ti});
    pos_t = interp1(mt.times, mt.pos, t, 'linear').';
    r_tx  = vecnorm(pos_t - Tx_Pos);
    r_ref = vecnorm(pos_t - RefRxPos);
    total_m = r_tx + r_ref;
    gt_r(:, ti) = (total_m - Baseline_m).' / 1000;
    d_path = zeros(1, numel(t));
    d_path(1) = total_m(2) - total_m(1);
    d_path(end) = total_m(end) - total_m(end-1);
    if numel(t) > 2
        d_path(2:end-1) = (total_m(3:end) - total_m(1:end-2)) / 2;
    end
    gt_d(:, ti) = (-d_path ./ (1 * lambda)).';
end

% ---- Load FERS ----
fprintf('Loading FERS baseband...\n');
[Ino,  Qno,  scale_no ] = loadfersHDF5('01_FERS/direct_3demo.h5');
[Imov, Qmov, scale_mov] = loadfersHDF5('01_FERS/echo_3demo.h5');
I_Qno  = (Ino  + 1i*Qno ) .* scale_no;
I_Qmov = (Imov + 1i*Qmov) .* scale_mov;

% ---- Axes ----
Ndelay = floor(delay * fs);
range_axis = (0:Ndelay/fs) * c / 1000;
dopp_axis  = -dopp_bins:1:dopp_bins;
gt_colors  = lines(NT);

% ---- MTT ----
% filterType 2 = CSKF (matches Scheme B2 tuning)
mtt = multiTargetTracker(4, 6, 15, 2);

% ---- Figure ----
fig = figure('Color','w','Position',[50 50 1300 550]);
ax_left  = axes('Position', [0.03 0.05 0.43 0.92]);
ax_right = axes('Position', [0.50 0.05 0.47 0.92]);

% Left: ARD (labels set each frame inside loop)
y_db_full = NaN(size(range_axis,2), 2*dopp_bins+1);
imagesc(ax_left, range_axis, dopp_axis, y_db_full);
set(ax_left, 'YDir','normal');
colormap(ax_left, turbo); caxis(ax_left, [-50 0]);

% Right: phosphor canvas
imagesc(ax_right, range_axis, dopp_axis, NaN(size(range_axis,2), 2*dopp_bins+1));
set(ax_right, 'YDir','normal');
colormap(ax_right, turbo); caxis(ax_right, [-50 0]);
xlabel(ax_right, 'Bistatic Range [km]');
ylabel(ax_right, 'Bistatic Doppler [Hz]');
hold(ax_right, 'on');  % must be on before creating persistent handles below

% History scatter handles
h_history = scatter(ax_right, [], [], 36, 'filled', ...
                    'Marker', 'h', ...
                    'MarkerFaceColor', [0.8 0.8 0.1], ...
                    'MarkerEdgeColor', [0.8 0.8 0.1], ...
                    'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.25, ...
                    'DisplayName', 'Centroid history', ...
                    'HandleVisibility','off');
h_current = scatter(ax_right, [], [], 100, 'filled', ...
                    'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'w', ...
                    'DisplayName', 'CFAR centroid');
h_cfars   = scatter(ax_right, [], [], 4, 'r', 'filled', ...
                    'HandleVisibility','off');
h_gt      = gobjects(NT, 1);
for j = 1:NT
    h_gt(j) = plot(ax_right, gt_r(1,j), gt_d(1,j), 's', ...
                   'MarkerSize', 8, 'MarkerFaceColor', gt_colors(j,:), ...
                   'MarkerEdgeColor', 'k', ...
                   'DisplayName', sprintf('GT Target %d', j));
end

% Legend dummy handles for track states (real handles created per-frame;
% these dummies keep a stable legend entry)
h_leg_conf = plot(ax_right, NaN, NaN, '^', 'MarkerSize', 9, ...
                  'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k', ...
                  'LineStyle', 'none', 'DisplayName', 'Confirmed track');
h_leg_tent = plot(ax_right, NaN, NaN, '^', 'MarkerSize', 7, ...
                  'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'b', ...
                  'LineStyle', 'none', 'LineWidth', 1.4, ...
                  'DisplayName', 'Tentative track');

legend(ax_right, 'Location', 'northwest', 'FontSize', 8, 'Box', 'on', ...
       'TextColor', 'w', 'Color', [0.1 0.1 0.1], 'EdgeColor', [0.4 0.4 0.4]);

h_timer_text = text(0.05, 0.95, '', 'Parent', ax_right, ...
                    'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold', ...
                    'Color', 'w', 'BackgroundColor', [0 0 0 0.6]);

% Track overlays (redrawn each frame via delete/plot)
track_plot_handles = gobjects(0);

% ---- History buffer ----
max_hist = 10000;
hist_r = zeros(1, max_hist);
hist_d = zeros(1, max_hist);
n_hist = 0;

% ---- Animate ----
N = 60;
for k = 1:N
    initial = (k-1)*fs + 1;
    current = k * fs;

    s1 = I_Qmov(initial:current);
    s2 = I_Qno (initial:current);
    s1 = procECA_Optimized(s2, s1, proc);
    [y, ~] = ardNoPlot(s1, s2, fs, dopp_bins, delay, k, []);
    y_mag = abs(y);
    y_db  = 20*log10(y_mag / max(y_mag(:)) + eps);

    % ---- Left: ARD ----
    imagesc(ax_left, range_axis, dopp_axis, y_db);   % hold off → clears previous frame
    hold(ax_left, 'on');
    set(ax_left, 'YDir','normal');
    colormap(ax_left, turbo); caxis(ax_left, [-50 0]);
    xlabel(ax_left, 'Bistatic Range [km]');
    ylabel(ax_left, 'Bistatic Doppler [Hz]');
    % GT markers for all three targets
    for j = 1:NT
        plot(ax_left, gt_r(k,j), gt_d(k,j), 'w*', 'MarkerSize', 8, 'LineWidth', 1.2);
    end
    hold(ax_left, 'off');

    % ---- Right: CFAR + centroids + phosphor history ----
    [targetClusters, ~, ~] = ca_cfar(y.', 1e-7, fs, dopp_bins, delay, 20);
    [clusterCentroids, ~, ~, ~] = meanShift(targetClusters, 3, 5);

    % CFAR faint
    if ~isempty(targetClusters) && size(targetClusters,2) > 0
        h_cfars.XData = targetClusters(1, :);
        h_cfars.YData = targetClusters(2, :);
        h_cfars.Visible = 'on';
    else
        h_cfars.Visible = 'off';
    end

    % Current centroids (bright red hexagons)
    if ~isempty(clusterCentroids)
        cents = clusterCentroids;                  % 2×N: row1=range, row2=doppler
        if size(cents,1) == 2, cents = cents.'; end % → N×2
        if size(cents,1) > 0
            h_current.XData = cents(:,1);
            h_current.YData = cents(:,2);
            h_current.Visible = 'on';
        else
            h_current.Visible = 'off';
        end
    else
        h_current.Visible = 'off';
    end

    % ---- MTT pipeline ----
    if ~isempty(clusterCentroids)
        mtt = mtt.createNewTracks(clusterCentroids, k);
    end
    mtt = mtt.maintainTracks();
    mtt = mtt.predictionStage();
    if ~isempty(clusterCentroids)
        mtt = mtt.updateStage(clusterCentroids, k);
    end

    % Overlay track estimates on right panel
    delete(track_plot_handles);
    track_plot_handles = gobjects(0);
    for ti = 1:numel(mtt.tracks)
        tr = mtt.tracks(ti);
        if tr.deleted || isempty(tr.predictedTrack), continue; end
        r_est = tr.predictedTrack(1, end);
        d_est = tr.predictedTrack(2, end);
        if tr.confirmed
            h_tr = plot(ax_right, r_est, d_est, 'c^', ...
                        'MarkerSize', 9, 'LineWidth', 1.5, ...
                        'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k');
            h_lbl = text(r_est + 0.3, d_est + 3, sprintf('T%d', tr.trackId), ...
                         'Parent', ax_right, 'Color', 'c', 'FontSize', 8, 'FontWeight', 'bold');
            track_plot_handles(end+1) = h_lbl;
        else
            h_tr = plot(ax_right, r_est, d_est, 'b^', ...
                        'MarkerSize', 7, 'LineWidth', 1.2, ...
                        'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'b');
        end
        track_plot_handles(end+1) = h_tr;
    end

    % Phosphor history (yellow hexagons, fixed alpha)
    if n_hist > 0
        h_history.XData = hist_r(1:n_hist);
        h_history.YData = hist_d(1:n_hist);
        h_history.Visible = 'on';
    else
        h_history.Visible = 'off';
    end

    % GT squares
    for j = 1:NT
        h_gt(j).XData = gt_r(k, j);
        h_gt(j).YData = gt_d(k, j);
    end

    % Text overlay
    h_timer_text.String = sprintf('Scan %d / %d   t = %.1f s', k, N, k);

    % Store current centroids in history
    if ~isempty(clusterCentroids)
        cents = clusterCentroids;
        if size(cents,1) == 2, cents = cents.'; end
        n_new = size(cents,1);
        if n_hist + n_new <= max_hist
            hist_r(1, n_hist+1:n_hist+n_new) = cents(:,1).';
            hist_d(1, n_hist+1:n_hist+n_new) = cents(:,2).';
            n_hist = n_hist + n_new;
        end
    end

    drawnow;
end

fprintf('Done. %d frames, %d history centroids.\n', N, n_hist);
fprintf('Right panel shows phosphor trail of all past centroids.\n');
