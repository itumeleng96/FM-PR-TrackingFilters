function cache_seeded_clusters(N_MC, use_parfor)
% CACHE_SEEDED_CLUSTERS  Extract mean-shift cluster centroids per (seed, scenario).
%
%   cache_seeded_clusters(N_MC, use_parfor)
%
% For each seed s = 1..N_MC and scenario in {landing, takeoff, orbit360},
% loads seeds/seed_SSS/{direct,echo}_<scen>.h5, runs ECA + CFAR + mean-shift,
% and saves seeds/seed_SSS/clusters_<scen>.mat with fields:
%   scans          : 1xN cell, scans{k} = 2xM centroid matrix
%   true_r, true_d : ground-truth arrays (scenario-fixed)
%   N_scans
%
% Skips (seed, scenario) pairs already cached (idempotent).
%
% use_parfor : true to parfor over seeds (requires PCT). Default false.

    if nargin < 1, N_MC = 3; end
    if nargin < 2, use_parfor = false; end

    addpath('FERS/', 'cfar/', 'meanShiftCluster/', 'DPI_Suppression/', ...
            'groundTruthCalculations/');

    scenarios = {'landing','takeoff','orbit360','3targets'};
    scen_secs = struct('landing',120,'takeoff',60,'orbit360',120,'x3targets',60);
    % Note: struct field can't start with a digit; use x3targets for lookup and remap.
    scen_secs_lookup = @(s) scen_secs.(regexprep(s, '^(\d)', 'x$1'));

    fs = 200000; dopp_bins = 200; delay = 233e-6; c = 299792458;
    range_delay = delay * c;
    proc = struct('cancellationMaxRange_m',range_delay, ...
                  'cancellationMaxDoppler_Hz',4, ...
                  'TxToRefRxDistance_m',12540, ...
                  'nSegments',1,'nIterations',20, ...
                  'Fs',fs,'alpha',0,'initialAlpha',0);

    % Pre-compute ground truth per scenario (same for all seeds).
    gt = struct();
    for si = 1:numel(scenarios)
        scen_key = regexprep(scenarios{si}, '^(\d)', 'x$1');
        [~, tr, td] = computeGroundTruth(scenarios{si});
        gt.(scen_key) = struct('true_r',tr,'true_d',td);
    end

    if use_parfor
        parfor s = 1:N_MC
            process_one_seed(s, scenarios, scen_secs_lookup, fs, dopp_bins, delay, proc, gt);
        end
    else
        for s = 1:N_MC
            process_one_seed(s, scenarios, scen_secs_lookup, fs, dopp_bins, delay, proc, gt);
        end
    end
end

function process_one_seed(s, scenarios, scen_secs_lookup, fs, dopp_bins, delay, proc, gt)
    seed_pad = sprintf('%03d', s);
    seed_dir = fullfile('seeds', ['seed_' seed_pad]);
    if ~isfolder(seed_dir)
        warning('Missing %s; skipping seed %d', seed_dir, s);
        return;
    end
    fprintf('[seed %s] extracting centroids...\n', seed_pad);
    t_seed = tic;

    for si = 1:numel(scenarios)
        scen  = scenarios{si};
        outmat = fullfile(seed_dir, sprintf('clusters_%s.mat', scen));
        if isfile(outmat)
            continue;    % idempotent
        end
        direct_h5 = fullfile(seed_dir, sprintf('direct_%s.h5', scen));
        echo_h5   = fullfile(seed_dir, sprintf('echo_%s.h5',   scen));
        if ~isfile(direct_h5) || ~isfile(echo_h5)
            warning('Missing H5 for seed=%s scen=%s', seed_pad, scen);
            continue;
        end
        N_scans = scen_secs_lookup(scen);
        scen_key = regexprep(scen, '^(\d)', 'x$1');
        true_r  = gt.(scen_key).true_r;
        true_d  = gt.(scen_key).true_d;

        [Ino,  Qno,  scale_no ] = loadfersHDF5(direct_h5);
        [Imov, Qmov, scale_mov] = loadfersHDF5(echo_h5);
        I_Qmov = (Imov + 1i*Qmov) .* scale_mov;
        I_Qno  = (Ino  + 1i*Qno ) .* scale_no;

        scans = cell(1, N_scans);
        initial = 1; current = fs;
        for k = 1:N_scans
            s1 = I_Qmov(initial:current);
            s2 = I_Qno (initial:current);
            s1 = procECA_Optimized(s2, s1, proc);
            [y, ~] = ardNoPlot(s1, s2, fs, dopp_bins, delay, k, []);
            [tc, ~, ~]  = ca_cfar(y.', 1e-7, fs, dopp_bins, delay, 20);
            [cc, ~, ~, ~] = meanShift(tc, 3, 5);
            scans{k} = cc;
            initial = current + 1;
            current = current + fs;
        end
        save(outmat, 'scans', 'true_r', 'true_d', 'N_scans');
        fprintf('  [seed %s] %-9s done (%d scans)\n', seed_pad, scen, N_scans);
    end
    fprintf('[seed %s] total %.1fs\n', seed_pad, toc(t_seed));
end
