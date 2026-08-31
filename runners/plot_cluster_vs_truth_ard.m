function plot_cluster_vs_truth_ard(scenario_name, snapshot_scans)
% PLOT_CLUSTER_VS_TRUTH_ARD
%   For a scenario, produce a montage of ARD snapshots with the mean-shift
%   centroid (measurement) and ground truth overlaid on each panel.
%   Lets you visually check measurement quality vs truth.
%
%   plot_cluster_vs_truth_ard('landing')
%   plot_cluster_vs_truth_ard('orbit360', [10 40 70 100])
%
% Requires direct_<scenario>.h5 and echo_<scenario>.h5.
% Saves figures/ard_cluster_vs_truth_<scenario>.png.

    addpath('FERS/', 'cfar/', 'meanShiftCluster/', 'DPI_Suppression/', ...
            'groundTruthCalculations/');
    if ~exist('figures','dir'); mkdir('figures'); end

    scen_secs = struct('levelFlight',60,'landing',120,'takeoff',60,'orbit360',120);
    N_scans   = scen_secs.(scenario_name);
    direct_h5 = sprintf('direct_%s.h5', scenario_name);
    echo_h5   = sprintf('echo_%s.h5',   scenario_name);

    fs = 200000; dopp_bins = 200; delay = 233e-6; c = 299792458;
    range_delay = delay*c;
    proc = struct('cancellationMaxRange_m',range_delay,'cancellationMaxDoppler_Hz',4, ...
                  'TxToRefRxDistance_m',12540,'nSegments',1,'nIterations',20, ...
                  'Fs',fs,'alpha',0,'initialAlpha',0);

    if nargin < 2 || isempty(snapshot_scans)
        snapshot_scans = round(linspace(5, N_scans-5, 4));
    end
    n_snap = numel(snapshot_scans);

    [~, true_r, true_d] = computeGroundTruth(scenario_name);
    N_gt = numel(true_r);

    [Ino,  Qno,  scale_no ] = loadfersHDF5(direct_h5);
    [Imov, Qmov, scale_mov] = loadfersHDF5(echo_h5);
    I_Qmov = (Imov + 1i*Qmov) .* scale_mov;
    I_Qno  = (Ino  + 1i*Qno ) .* scale_no;

    Ndelay = floor(delay*fs); Ndop = ceil(fs*dopp_bins/fs);
    range_axis = (0:Ndelay).' * (c/fs) / 1000;
    dopp_axis  = (-Ndop:Ndop).';

    fig = figure('Position',[80 80 400*n_snap 460],'Visible','off');

    for ii = 1:n_snap
        k = snapshot_scans(ii);
        seg_start = (k-1)*fs + 1;
        seg_end   = k*fs;
        s1 = I_Qmov(seg_start:seg_end);
        s2 = I_Qno (seg_start:seg_end);
        s1 = procECA_Optimized(s2, s1, proc);
        [y, ~] = ardNoPlot(s1, s2, fs, dopp_bins, delay, k, []);

        [targetClusters, ~, ~] = ca_cfar(y.', 1e-7, fs, dopp_bins, delay, 20);
        [clusterCentroids, ~, ~, ~] = meanShift(targetClusters, 3, 5);

        subplot(1, n_snap, ii);
        imagesc(dopp_axis, range_axis, 10*log10(y + eps));
        axis xy; colormap(jet); hold on;
        xlabel('Bistatic doppler (Hz)'); ylabel('Bistatic range (km)');
        title(sprintf('scan %d', k));

        % Overlay ground truth (green star)
        if k <= N_gt
            plot(true_d(k), true_r(k), 'g*', 'MarkerSize', 14, 'LineWidth', 2);
        end
        % Overlay mean-shift centroid (top-vote first after sort A)
        if ~isempty(clusterCentroids)
            % Highlight strongest centroid (first after sort) in white 'o'
            plot(clusterCentroids(2,1), clusterCentroids(1,1), 'wo', ...
                 'MarkerSize', 12, 'LineWidth', 2);
            % Show any additional centroids as smaller white circles
            if size(clusterCentroids, 2) > 1
                plot(clusterCentroids(2,2:end), clusterCentroids(1,2:end), 'w.', ...
                     'MarkerSize', 12);
            end
        end
        colorbar;
        if ii == 1
            legend('truth (*)', 'meas #1 (o)', 'other clusters (.)', ...
                   'Location','southoutside','TextColor','w','Color',[0 0 0 0.5]);
        end
    end

    sgtitle(sprintf('Mean-shift centroids vs Ground Truth on ARD | %s', scenario_name), ...
            'FontWeight','bold','FontSize',13);
    out_png = sprintf('figures/ard_cluster_vs_truth_%s.png', scenario_name);
    exportgraphics(fig, out_png, 'Resolution', 130);
    close(fig);
    fprintf('Saved %s\n', out_png);
end
