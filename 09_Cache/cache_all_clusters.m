% =============================================================
% cache_all_clusters.m
%
% One-shot cache of ALL mean-shift cluster centroids per scan for
% each scenario, so downstream tuning can run GNN gating without
% re-running ECA/CFAR/meanShift.
% Saves clusters_<scenario>.mat with fields:
%   scans          : 1xN cell array, scans{k} = 2xM matrix of centroids
%   true_r, true_d : ground-truth arrays
% =============================================================

clc;

addpath('01_FERS/', '04_Detection/', '05_Clustering/', '03_DPI_Cancellation/', ...
        '07_Evaluation/GroundTruth/');

scenarios = {'levelFlight','landing','orbit360'};
scen_secs = struct('levelFlight',60,'landing',120,'orbit360',120);

fs = 200000; dopp_bins = 200; delay = 233e-6; c = 299792458;
range_delay = delay*c;
proc = struct('cancellationMaxRange_m',range_delay,'cancellationMaxDoppler_Hz',4, ...
              'TxToRefRxDistance_m',12540,'nSegments',1,'nIterations',20, ...
              'Fs',fs,'alpha',0,'initialAlpha',0);

for s = 1:numel(scenarios)
    name = scenarios{s};
    N_scans = scen_secs.(name);
    direct_h5 = sprintf('direct_%s.h5', name);
    echo_h5   = sprintf('echo_%s.h5',   name);

    [~, true_r, true_d] = computeGroundTruth(name);

    fprintf('\n[%s] loading FERS baseband...\n', name);
    [Ino,  Qno,  scale_no ] = loadfersHDF5(direct_h5);
    [Imov, Qmov, scale_mov] = loadfersHDF5(echo_h5);
    I_Qmov = (Imov + 1i*Qmov) .* scale_mov;
    I_Qno  = (Ino  + 1i*Qno ) .* scale_no;

    scans = cell(1, N_scans);
    initial = 1; current = fs;
    t0 = tic;
    for k = 1:N_scans
        s1 = I_Qmov(initial:current);
        s2 = I_Qno (initial:current);
        s1 = procECA_Optimized(s2, s1, proc);
        [y, ~] = ardNoPlot(s1, s2, fs, dopp_bins, delay, k, []);
        [targetClusters, ~, ~]      = ca_cfar(y.', 1e-7, fs, dopp_bins, delay, 20);
        [clusterCentroids, ~, ~, ~] = meanShift(targetClusters, 3, 5);
        scans{k} = clusterCentroids;   % 2xM (may be empty)
        if mod(k, 20) == 0
            fprintf('  scan %3d/%d   t=%.1fs\n', k, N_scans, toc(t0));
        end
        initial = current + 1;
        current = current + fs;
    end

    outname = sprintf('clusters_%s.mat', name);
    save(outname, 'scans', 'true_r', 'true_d', 'N_scans');
    fprintf('[%s] saved %s   (%d scans, %.1fs total)\n', name, outname, N_scans, toc(t0));
end

fprintf('\nDone. Caches: clusters_levelFlight.mat, clusters_landing.mat, clusters_orbit360.mat\n');
