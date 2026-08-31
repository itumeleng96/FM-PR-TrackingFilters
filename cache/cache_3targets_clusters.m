% =============================================================
% cache_3targets_clusters.m
%
% One-shot cache of mean-shift cluster centroids per scan for the
% 3-target scenario (FERS/BackupScenarios/scenario_N3_targets.fersxml,
% h5 pair: direct_N3.h5 / echo_N3.h5).
%
% Saves clusters_3targets.mat with fields:
%   scans          : 1xN cell array, scans{k} = 2xM matrix of centroids
%   true_r, true_d : N x 3 ground-truth arrays (columns = targets)
%   N_scans        : number of scans
%   n_range_bins, n_doppler_bins : ARD dimensions (for Poisson FA rate)
% =============================================================
clc;

addpath('FERS/', 'cfar/', 'meanShiftCluster/', 'DPI_Suppression/', ...
        'groundTruthCalculations/');

name       = '3targets';
N_scans    = 60;
direct_h5  = 'direct_N3.h5';
echo_h5    = 'echo_N3.h5';

fs = 200000; dopp_bins = 200; delay = 233e-6; c = 299792458;
range_delay = delay * c;
proc = struct('cancellationMaxRange_m', range_delay, ...
              'cancellationMaxDoppler_Hz', 4, ...
              'TxToRefRxDistance_m', 12540, ...
              'nSegments', 1, 'nIterations', 20, ...
              'Fs', fs, 'alpha', 0, 'initialAlpha', 0);

[~, true_r, true_d] = computeGroundTruth(name);   % (61 x 3) each

fprintf('[%s] loading FERS baseband...\n', name);
[Ino,  Qno,  scale_no ] = loadfersHDF5(direct_h5);
[Imov, Qmov, scale_mov] = loadfersHDF5(echo_h5);
I_Qmov = (Imov + 1i*Qmov) .* scale_mov;
I_Qno  = (Ino  + 1i*Qno ) .* scale_no;

scans = cell(1, N_scans);
n_range_bins = 0; n_doppler_bins = 0;
initial = 1; current = fs;
t0 = tic;
for k = 1:N_scans
    s1 = I_Qmov(initial:current);
    s2 = I_Qno (initial:current);
    s1 = procECA_Optimized(s2, s1, proc);
    [y, ~] = ardNoPlot(s1, s2, fs, dopp_bins, delay, k, []);
    if k == 1
        n_range_bins   = size(y, 2);   % after transpose in ca_cfar, this is what CFAR sees
        n_doppler_bins = size(y, 1);
    end
    [targetClusters, ~, ~]      = ca_cfar(y.', 1e-7, fs, dopp_bins, delay, 20);
    [clusterCentroids, ~, ~, ~] = meanShift(targetClusters, 3, 5);
    scans{k} = clusterCentroids;
    if mod(k, 10) == 0
        fprintf('  scan %3d/%d   t=%.1fs\n', k, N_scans, toc(t0));
    end
    initial = current + 1;
    current = current + fs;
end

outname = 'clusters_3targets.mat';
save(outname, 'scans', 'true_r', 'true_d', 'N_scans', ...
     'n_range_bins', 'n_doppler_bins');
fprintf('[%s] saved %s (%d scans, %.1fs total, ARD %dx%d)\n', ...
        name, outname, N_scans, toc(t0), n_range_bins, n_doppler_bins);
