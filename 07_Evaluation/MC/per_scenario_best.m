% =============================================================
% per_scenario_best.m
%
% For each (filter, scenario) pair, pick the tuning combo that
% minimises that pair's own normalised RMSE score from the Scheme B2
% grid. Report the per-scenario tuning table and the resulting
% per-scenario RMSE.
%
% Rationale: The globally-best tuning trades performance across
% scenarios and can leave hard scenarios (landing, orbit360) with
% divergent doppler for the KF-family. Per-scenario tuning shows the
% best each filter can achieve when allowed to specialise.
% =============================================================

clc; close all;

load('scheme_b2_results.mat', 'results', 'combos', 'scenarios', 'filters');
n_filt = numel(filters);
n_scen = size(scenarios, 1);

fprintf('\nPer-scenario best tuning (Scheme B2 grid):\n');
fprintf('%s\n', repmat('=', 1, 88));
fprintf('%-8s %-12s | %8s %8s %10s %10s | %10s %10s\n', ...
        'Filter','Scenario','cs_alpha','cs_beta','Q_scale_r','Q_scale_d','RMSE_r','RMSE_d');
fprintf('%s\n', repmat('-', 1, 88));

per_scen_best = struct();

for f = 1:n_filt
    R = results{f};                     % n_combo x n_scen x 2  (rmse_r, rmse_d)
    for s = 1:n_scen
        r_col = R(:, s, 1);
        d_col = R(:, s, 2);
        r_lo  = min(r_col, [], 'omitnan');   r_hi = max(r_col, [], 'omitnan');
        d_lo  = min(d_col, [], 'omitnan');   d_hi = max(d_col, [], 'omitnan');
        r_norm = (r_col - r_lo) ./ max(r_hi - r_lo, eps);
        d_norm = (d_col - d_lo) ./ max(d_hi - d_lo, eps);
        score  = (r_norm + d_norm) / 2;
        [~, cbest] = min(score);

        per_scen_best.(filters{f}).(scenarios{s,1}) = struct( ...
            'cs_alpha',  combos(cbest, 1), ...
            'cs_beta',   combos(cbest, 2), ...
            'Q_scale_r', combos(cbest, 3), ...
            'Q_scale_d', combos(cbest, 4), ...
            'rmse_r',    R(cbest, s, 1), ...
            'rmse_d',    R(cbest, s, 2));

        fprintf('%-8s %-12s | %8.2f %8.2f %10.2f %10.2f | %8.2f km  %6.1f Hz\n', ...
                filters{f}, scenarios{s,1}, ...
                combos(cbest,1), combos(cbest,2), combos(cbest,3), combos(cbest,4), ...
                R(cbest, s, 1), R(cbest, s, 2));
    end
    if f < n_filt; fprintf('%s\n', repmat('-', 1, 88)); end
end
fprintf('%s\n', repmat('=', 1, 88));

save('per_scenario_best.mat', 'per_scen_best');
fprintf('\nSaved per-scenario best tunings to per_scenario_best.mat\n');
