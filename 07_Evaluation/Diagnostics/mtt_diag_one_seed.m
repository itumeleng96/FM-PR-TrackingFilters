function mtt_diag_one_seed(seed, filter_name)
% MTT_DIAG_ONE_SEED  Verify what TTD is actually counting.
%
%   mtt_diag_one_seed(seed, filter_name)
%
% For a single seed and filter, prints per-scan:
%   scan | confirmed_ids | ids_matched_to_any_GT | new_confirmations | new_deletions
% and finally lists which ever-matched IDs are missing at end.

    if nargin < 1, seed = 1; end
    if nargin < 2, filter_name = 'PF'; end

    addpath('06_Tracking/MTT/', ...
            '06_Tracking/Filters/KF/','06_Tracking/Filters/UKF/', ...
            '06_Tracking/Filters/RGNF/','06_Tracking/Filters/PF/');

    ftype_map = struct('KF',2,'UKF',6,'RGNF',8,'PF',4);
    ft = ftype_map.(filter_name);

    seed_dir = fullfile('seeds', sprintf('seed_%03d', seed));
    d = load(fullfile(seed_dir, 'clusters_3targets_pfa7.mat'));

    N  = d.N_scans; NT = size(d.true_r, 2);
    assoc_dr_km = 2.0; assoc_dd_hz = 5.0;

    tracker = multiTargetTracker(3, 4, 25, ft);
    prev_confirmed_ids = [];
    ever_confirmed_ids = [];
    ever_matched_ids   = [];
    alive_at_end       = [];

    fprintf('seed=%d filter=%s\n', seed, filter_name);
    fprintf('scan | confirmed_ids               | matched_ids     | new_conf | deleted_since_prev\n');

    for k = 1:N
        det = d.scans{k};
        if k == 1
            tracker = tracker.createNewTracks(det, k);
        end
        tracker = tracker.predictionStage();
        tracker = tracker.updateStage(det, k);
        tracker = tracker.maintainTracks();

        cur_conf = [];
        matched_now = [];
        if ~isempty(tracker.tracks)
            for ti = 1:numel(tracker.tracks)
                trk = tracker.tracks(ti);
                if trk.confirmed && ~trk.deleted && ~isempty(trk.trueTrack)
                    tid = trk.trackId;
                    cur_conf(end+1) = tid; %#ok<AGROW>
                    if ~ismember(tid, ever_confirmed_ids)
                        ever_confirmed_ids(end+1) = tid; %#ok<AGROW>
                    end
                    tstate = trk.trueTrack(:, end);
                    for g = 1:NT
                        if k > size(d.true_r, 1); continue; end
                        gr = d.true_r(k, g); gd = d.true_d(k, g);
                        if abs(tstate(1) - gr) < assoc_dr_km && ...
                           abs(tstate(2) - gd) < assoc_dd_hz
                            matched_now(end+1) = tid; %#ok<AGROW>
                            if ~ismember(tid, ever_matched_ids)
                                ever_matched_ids(end+1) = tid; %#ok<AGROW>
                            end
                        end
                    end
                end
            end
        end
        matched_now = unique(matched_now);

        new_conf   = setdiff(cur_conf, prev_confirmed_ids);
        deleted    = setdiff(prev_confirmed_ids, cur_conf);
        if ~isempty(new_conf) || ~isempty(deleted) || mod(k,10)==0 || k==N
            fprintf('%3d  | %-28s| %-15s | %-8s | %s\n', k, ...
                    mat2str(sort(cur_conf)), mat2str(sort(matched_now)), ...
                    mat2str(sort(new_conf)), mat2str(sort(deleted)));
        end
        prev_confirmed_ids = cur_conf;
        if k == N, alive_at_end = cur_conf; end
    end

    fprintf('\nSummary for seed %d filter %s:\n', seed, filter_name);
    fprintf('  ever_confirmed_ids : %s\n', mat2str(sort(ever_confirmed_ids)));
    fprintf('  ever_matched_ids   : %s\n', mat2str(sort(ever_matched_ids)));
    fprintf('  alive_at_end       : %s\n', mat2str(sort(alive_at_end)));
    fprintf('  TTC = # GT matched : %d\n', numel(unique(ever_matched_ids))); % (upper bound)
    fprintf('  FTC = confirmed but never matched : %d\n', ...
            numel(setdiff(ever_confirmed_ids, ever_matched_ids)));
    fprintf('  TTD = matched IDs not alive at end : %d\n', ...
            numel(setdiff(ever_matched_ids, alive_at_end)));
end
