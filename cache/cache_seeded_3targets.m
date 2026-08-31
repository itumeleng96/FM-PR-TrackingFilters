function cache_seeded_3targets(N_MC, use_parfor, pfa)
% CACHE_SEEDED_3TARGETS  Extract mean-shift centroids for the 3-target
% scenario across all seeded FERS runs at a given CFAR P_fa.
%
%   cache_seeded_3targets(N_MC, use_parfor, pfa)
%
% For each seed s = 1..N_MC, loads seeds/seed_SSS/{direct,echo}_3targets.h5,
% runs ECA + CFAR (at pfa) + mean-shift, and saves
% seeds/seed_SSS/clusters_3targets_pfa<E>.mat  (E = |log10(pfa)|)
% with fields:
%   scans          : 1xN cell, scans{k} = 2xM centroid matrix
%   true_r, true_d : ground-truth (N_scans x 3), one column per target
%   N_scans, pfa
%   n_range_bins, n_doppler_bins : ARD dimensions
%
% Skips (seed) pairs already cached (idempotent).
%
% use_parfor : true to parfor over seeds (requires PCT). Default false.
% pfa        : CFAR probability of false alarm. Default 1e-5.

    if nargin < 1, N_MC = 100;     end
    if nargin < 2, use_parfor = false; end
    if nargin < 3, pfa = 1e-5;    end

    addpath('FERS/', 'cfar/', 'meanShiftCluster/', 'DPI_Suppression/', ...
            'groundTruthCalculations/');

    fs = 200000; dopp_bins = 200; delay = 233e-6; c = 299792458;
    range_delay = delay * c;
    proc = struct('cancellationMaxRange_m',range_delay, ...
                  'cancellationMaxDoppler_Hz',4, ...
                  'TxToRefRxDistance_m',12540, ...
                  'nSegments',1,'nIterations',20, ...
                  'Fs',fs,'alpha',0,'initialAlpha',0);

    [~, true_r, true_d] = computeGroundTruth('3targets');   % 61 x 3
    pfa_tag = abs(round(log10(pfa)));

    if use_parfor
        parfor s = 1:N_MC
            process_seed(s, fs, dopp_bins, delay, proc, true_r, true_d, pfa, pfa_tag);
        end
    else
        for s = 1:N_MC
            process_seed(s, fs, dopp_bins, delay, proc, true_r, true_d, pfa, pfa_tag);
        end
    end
end

function process_seed(s, fs, dopp_bins, delay, proc, true_r, true_d, pfa, pfa_tag)
    seed_pad  = sprintf('%03d', s);
    seed_dir  = fullfile('seeds', ['seed_' seed_pad]);
    outmat    = fullfile(seed_dir, sprintf('clusters_3targets_pfa%d.mat', pfa_tag));
    direct_h5 = fullfile(seed_dir, 'direct_3targets.h5');
    echo_h5   = fullfile(seed_dir, 'echo_3targets.h5');

    if isfile(outmat)
        return;
    end
    if ~isfile(direct_h5) || ~isfile(echo_h5)
        warning('Missing H5 for seed=%s', seed_pad);
        return;
    end

    t_seed = tic;
    [Ino,  Qno,  scale_no ] = loadfersHDF5(direct_h5);
    [Imov, Qmov, scale_mov] = loadfersHDF5(echo_h5);
    I_Qmov = (Imov + 1i*Qmov) .* scale_mov;
    I_Qno  = (Ino  + 1i*Qno ) .* scale_no;

    N_scans = 60;
    scans = cell(1, N_scans);
    n_range_bins = 0; n_doppler_bins = 0;
    initial = 1; current = fs;
    for k = 1:N_scans
        s1 = I_Qmov(initial:current);
        s2 = I_Qno (initial:current);
        s1 = procECA_Optimized(s2, s1, proc);
        [y, ~] = ardNoPlot(s1, s2, fs, dopp_bins, delay, k, []);
        if k == 1
            n_range_bins   = size(y, 2);
            n_doppler_bins = size(y, 1);
        end
        [tc, ~, ~]  = ca_cfar(y.', pfa, fs, dopp_bins, delay, 20);
        [cc, ~, ~, ~] = meanShift(tc, 3, 5);
        scans{k} = cc;
        initial = current + 1;
        current = current + fs;
    end

    save(outmat, 'scans', 'true_r', 'true_d', 'N_scans', ...
         'n_range_bins', 'n_doppler_bins', 'pfa');
    fprintf('[seed %s] 3targets pfa=%.0e done (%d scans, %.1fs)\n', ...
            seed_pad, pfa, N_scans, toc(t_seed));
end
