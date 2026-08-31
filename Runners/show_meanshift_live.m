function show_meanshift_live(scenario_name)
% SHOW_MEANSHIFT_LIVE
%   Live-updating ARD figure with mean-shift centroids and ground truth
%   overlaid, scan-by-scan. Range on x-axis, doppler on y-axis (standard
%   FM PR presentation).
%
%   show_meanshift_live('landing')
%
% Requires direct_<scenario>.h5 and echo_<scenario>.h5.

    addpath('01_FERS/', '04_Detection/', '05_Clustering/', '03_DPI_Cancellation/', ...
            '07_Evaluation/GroundTruth/');

    scen_secs = struct('levelFlight',60,'landing',120,'takeoff',60,'orbit360',120);
    if ~isfield(scen_secs, scenario_name)
        error('Unknown scenario "%s".', scenario_name);
    end
    N_scans = scen_secs.(scenario_name);
    direct_h5 = sprintf('direct_%s.h5', scenario_name);
    echo_h5   = sprintf('echo_%s.h5',   scenario_name);

    fs = 200000; dopp_bins = 200; delay = 233e-6; c = 299792458;
    range_delay = delay*c;
    proc = struct('cancellationMaxRange_m',range_delay,'cancellationMaxDoppler_Hz',4, ...
                  'TxToRefRxDistance_m',12540,'nSegments',1,'nIterations',20, ...
                  'Fs',fs,'alpha',0,'initialAlpha',0);

    [~, true_r, true_d] = computeGroundTruth(scenario_name);
    N_gt = numel(true_r);

    fprintf('Loading FERS data...\n');
    [Ino,  Qno,  scale_no ] = loadfersHDF5(direct_h5);
    [Imov, Qmov, scale_mov] = loadfersHDF5(echo_h5);
    I_Qmov = (Imov + 1i*Qmov) .* scale_mov;
    I_Qno  = (Ino  + 1i*Qno ) .* scale_no;

    Ndelay = floor(delay*fs); Ndop = ceil(fs*dopp_bins/fs);
    range_axis_km = (0:Ndelay).' * (c/fs) / 1000;   % x-axis
    dopp_axis_hz  = (-Ndop:Ndop).';                 % y-axis

    % Live figure (visible)
    fig = figure('Name', sprintf('Meanshift live | %s', scenario_name), ...
                 'Position', [80 80 1000 620]);
    initial = 1; current = fs;

    meas_r_hist = nan(1, N_scans);
    meas_d_hist = nan(1, N_scans);

    for k = 1:N_scans
        s1 = I_Qmov(initial:current);
        s2 = I_Qno (initial:current);
        s1 = procECA_Optimized(s2, s1, proc);
        [y, ~] = ardNoPlot(s1, s2, fs, dopp_bins, delay, k, []);

        [targetClusters, ~, ~] = ca_cfar(y.', 1e-7, fs, dopp_bins, delay, 20);
        [clusterCentroids, ~, ~, ~] = meanShift(targetClusters, 3, 5);

        if ~isempty(clusterCentroids)
            meas_r_hist(k) = clusterCentroids(1, 1);
            meas_d_hist(k) = clusterCentroids(2, 1);
        end

        % --- Plot (range on x, doppler on y — TRANSPOSE the ARD) ---
        clf(fig);
        imagesc(range_axis_km, dopp_axis_hz, 10*log10(y.' + eps));
        axis xy; colormap(jet); hold on;
        xlabel('Bistatic range (km)', 'FontSize', 12);
        ylabel('Bistatic doppler (Hz)', 'FontSize', 12);
        title(sprintf('%s | scan %d/%d  (range on x, doppler on y)', ...
                      scenario_name, k, N_scans), 'FontSize', 13);
        cb = colorbar; cb.Label.String = 'Level (dB)';

        % Ground-truth path so far (green line + current *)
        if N_gt >= 1
            k_gt = min(k, N_gt);
            plot(true_r(1:k_gt), true_d(1:k_gt), 'g-', 'LineWidth', 1.5);
            plot(true_r(k_gt),   true_d(k_gt),   'g*', 'MarkerSize', 16, 'LineWidth', 2);
        end

        % Measurement history (black trail) and current strongest (white o)
        valid = ~isnan(meas_r_hist(1:k));
        plot(meas_r_hist(valid), meas_d_hist(valid), 'k.-', 'MarkerSize', 6);
        if ~isempty(clusterCentroids)
            plot(clusterCentroids(1,1), clusterCentroids(2,1), 'wo', ...
                 'MarkerSize', 14, 'LineWidth', 2);
            if size(clusterCentroids, 2) > 1
                plot(clusterCentroids(1,2:end), clusterCentroids(2,2:end), ...
                     'w.', 'MarkerSize', 10);
            end
        end

        legend({'truth path','truth (*)','meas trail','strongest (o)','other clusters'}, ...
               'Location','northeast','TextColor','w','Color',[0 0 0 0.5], 'FontSize', 9);

        drawnow;
        fprintf('  scan %3d/%d   meas=[%.2f km, %.2f Hz]  truth=[%.2f km, %.2f Hz]\n', ...
                k, N_scans, meas_r_hist(k), meas_d_hist(k), ...
                (k<=N_gt)*true_r(min(k,N_gt)), (k<=N_gt)*true_d(min(k,N_gt)));

        initial = current + 1;
        current = current + fs;
    end
end
