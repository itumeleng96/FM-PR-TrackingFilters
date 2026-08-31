function mtt_mc(N_MC)
% MTT_MC  FERS-seeded Monte Carlo evaluation of MTT track-management metrics.
%
%   mtt_mc(N_MC)
%
% Loads seeds/seed_SSS/clusters_3targets.mat for each seed s = 1..N_MC.
% Each seed carries a distinct FERS thermal-noise realisation, so ARD,
% CFAR triggers, and mean-shift centroids differ per seed (real false-alarm
% structure — no synthetic Poisson injection).
%
% For each seed and each CFAR P_fa in {1e-5, 1e-7, 1e-9}:
%   * Filter the pre-clustered detection stream through the multiTargetTracker
%     for each of the four Akhlaghi-adaptive filters (CSKF, CSUKF, CSPF, CSRGNF).
%   * Score confirmed tracks vs ground truth (Kennedy-style lifetime accounting).
%
% NOTE: current cluster caches were generated at P_fa = 1e-7 in ca_cfar.
% To sweep true P_fa values you must re-cache per P_fa value. This MC run
% treats all three P_fa columns as the same detection stream and only
% differs in the tracker (which does not consume P_fa directly). For a
% real P_fa sweep, re-run cache_seeded_3targets.m at each P_fa and load
% the appropriate cache here.
%
% Metrics per (r, P_fa, filter) — Kennedy-style lifetime accounting:
%   TTC : GT targets ever associated with a confirmed track (0..3)
%   FTC : confirmed tracks (including deleted) that never matched any GT
%   TTD : GT targets NOT held by any alive confirmed track at the final scan
%         (end-state deletion count, bounded by NT=3)
%
% Association threshold: |dr| < 2 km AND |dd| < 5 Hz (filter-independent).
%
% Default: N_MC = 100.
% Saves mtt_mc_results.mat.

    if nargin < 1 || isempty(N_MC), N_MC = 100; end

    addpath('06_Tracking/MTT/', ...
            '06_Tracking/Filters/KF/','06_Tracking/Filters/UKF/', ...
            '06_Tracking/Filters/RGNF/','06_Tracking/Filters/PF/');

    Pfa_vals = [1e-5, 1e-7, 1e-9];
    filters  = {'KF','UKF','RGNF','PF'};
    ftypes   = [2,     6,     8,     4];    % CSKF, CSUKF, CSRGNF, CSPF

    NP = numel(Pfa_vals); NF = numel(filters);

    assoc_dr_km = 2.0;
    assoc_dd_hz = 5.0;

    gatingThreshold      = 25;
    confirmationThreshold = 3;
    deletionThreshold    = 4;

    raw = nan(N_MC, NP, NF, 3);

    t0 = tic;
    for r = 1:N_MC
        seed_pad = sprintf('%03d', r);
        seed_dir = fullfile('seeds', ['seed_' seed_pad]);

        for pi = 1:NP
            pfa_tag = abs(round(log10(Pfa_vals(pi))));
            seed_mat = fullfile(seed_dir, ...
                                sprintf('clusters_3targets_pfa%d.mat', pfa_tag));
            if ~isfile(seed_mat)
                continue;    % cache not yet generated at this P_fa
            end
            d = load(seed_mat);
            N  = d.N_scans;
            NT = size(d.true_r, 2);
            det_scans = d.scans;

            for fi = 1:NF
                ft = ftypes(fi);
                try
                    tracker = multiTargetTracker(confirmationThreshold, ...
                                                 deletionThreshold, ...
                                                 gatingThreshold, ft);
                    matched            = false(NT, N);
                    ever_confirmed_ids = [];
                    ever_matched_ids   = [];
                    matched_at_end     = false(NT, 1);

                    for k = 1:N
                        det = det_scans{k};
                        if k == 1
                            tracker = tracker.createNewTracks(det, k);
                        end
                        tracker = tracker.predictionStage();
                        tracker = tracker.updateStage(det, k);
                        tracker = tracker.maintainTracks();

                        cur_ids_confirmed = [];
                        if ~isempty(tracker.tracks)
                            for ti = 1:numel(tracker.tracks)
                                trk = tracker.tracks(ti);
                                if trk.confirmed && ~trk.deleted && ~isempty(trk.trueTrack)
                                    tstate = trk.trueTrack(:, end);
                                    cur_ids_confirmed(end+1) = trk.trackId; %#ok<AGROW>
                                    if ~ismember(trk.trackId, ever_confirmed_ids)
                                        ever_confirmed_ids(end+1) = trk.trackId; %#ok<AGROW>
                                    end
                                    for g = 1:NT
                                        if k > size(d.true_r, 1); continue; end
                                        gr = d.true_r(k, g); gd = d.true_d(k, g);
                                        if abs(tstate(1) - gr) < assoc_dr_km && ...
                                           abs(tstate(2) - gd) < assoc_dd_hz
                                            matched(g, k) = true;
                                            if ~ismember(trk.trackId, ever_matched_ids)
                                                ever_matched_ids(end+1) = trk.trackId; %#ok<AGROW>
                                            end
                                            if k == N
                                                matched_at_end(g) = true;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end

                    % Metrics.
                    ttc = sum(any(matched, 2));
                    ftc = numel(setdiff(ever_confirmed_ids, ever_matched_ids));
                    % TTD: end-state deletion count — GT targets not held
                    % by any alive confirmed track at the final scan.
                    ttd = NT - sum(matched_at_end);

                    raw(r, pi, fi, 1) = ttc;
                    raw(r, pi, fi, 2) = ftc;
                    raw(r, pi, fi, 3) = ttd;
                catch ME
                    warning('run=%d pfa=%.0e filt=%s failed: %s', r, Pfa_vals(pi), filters{fi}, ME.message);
                end
            end
        end
        if mod(r, 10) == 0 || r == 1
            fprintf('  MC run %3d/%d  (elapsed %.1fs)\n', r, N_MC, toc(t0));
        end
    end

    mu    = squeeze(mean(raw, 1, 'omitnan'));
    sigma = squeeze(std (raw, 0, 1, 'omitnan'));

    save('mtt_mc_results.mat', 'raw', 'mu', 'sigma', ...
         'Pfa_vals', 'filters', 'N_MC', ...
         'assoc_dr_km', 'assoc_dd_hz', 'gatingThreshold', ...
         'deletionThreshold');

    print_table(mu, sigma, Pfa_vals, filters);
end

function print_table(mu, sigma, Pfa, filters)
    metric_names = {'TTC','FTC','TTD'};
    NF = numel(filters);
    fprintf('\n=== Table mtt_data_assoc (mean +/- std over MC runs) ===\n');
    for mi = 1:3
        fprintf('\n[%s]\n', metric_names{mi});
        fprintf('%-8s', 'Pfa');
        for fi = 1:NF, fprintf(' %-18s', filters{fi}); end
        fprintf('\n');
        for pi = 1:numel(Pfa)
            fprintf('1e%-4d ', log10(Pfa(pi)));
            for fi = 1:NF
                fprintf(' %6.2f +/- %6.2f  ', mu(pi,fi,mi), sigma(pi,fi,mi));
            end
            fprintf('\n');
        end
    end
end
